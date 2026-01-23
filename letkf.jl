function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, ρ,
        localization, statetypes = parameters()
    @unpack modelproj = common_simulation()

    shuffle_xy = get(parameters(),:shuffle_xy, false)

    lengthscale = only(modelsteps).lengthscale

    paths = buildtruepaths()
    npaths = length(paths)
    
    gridshape = (length(y_grid), length(x_grid))  # useful later
    
    CI = CartesianIndices(gridshape)

    xy_grid = densify(x_grid, y_grid)
    trans = Proj.Transformation(modelproj, wgs84())
    lola = trans.(parent(parent(xy_grid)))

    distarr = lonlatgrid_dists(lola)
    gc = gaspari1999_410(distarr, compactlengthscale(lengthscale))

    # We wrap in Symmetric because we know it is. Without it, eigvals of β can be small
    # (but positive) and still fail the `isposdef` check
    hdistribution = Distributions.MvNormal(h0, Symmetric(LinearAlgebra.Diagonal(hB)*gc*LinearAlgebra.Diagonal(hB)))  # yes, w/ matrix argument we need variance
    bdistribution = Distributions.MvNormal(b0, Symmetric(LinearAlgebra.Diagonal(bB)*gc*LinearAlgebra.Diagonal(bB)))

    # Initial ensemble
    h_init = reshape(rand(rng, hdistribution, ens_size), gridshape..., ens_size)
    b_init = reshape(rand(rng, bdistribution, ens_size), gridshape..., ens_size)
     
    replace!(x->x < MIN_BETA ? MIN_BETA : x, b_init) #TODO Check structure here; GPT thinks this is not an appropriate way to do this
    replace!(x->x > MAX_BETA ? MAX_BETA : x, b_init) # possible correction: b_init[b_init .< MIN_BETA] .= MIN_BETA
    replace!(x->x < MIN_H ? MIN_H : x, h_init) 
    replace!(x->x > MAX_H ? MAX_H : x, h_init)

    locmask = anylocal(localization)
    h_init[CI[.!locmask],:] .= NaN
    b_init[CI[.!locmask],:] .= NaN

    xy_state = KeyedArray(fill(NaN, 2, length(y_grid), length(x_grid), ens_size, ntimes+1),
    field=[:h, :b], y=y_grid, x=x_grid, ens=1:ens_size, t=0:ntimes)

    xy_state(:h)(t=0) .= h_init
    xy_state(:b)(t=0) .= b_init

    if :tx in statetypes
        @unpack σNLK, σNML, shuffle_tx, NLKb, NMLb = parameters()

        NMLdistribution = Distributions.Normal(NMLb, σNML) #Nominal values from FWM dictionary
        NLKdistribution = Distributions.Normal(NLKb, σNLK)

        NML_init = rand(rng, NMLdistribution, ens_size)
        NLK_init = rand(rng, NLKdistribution, ens_size)

        tx_pwrs = KeyedArray(fill(NaN,2, ens_size, ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, t=0:ntimes)

        tx_pwrs(:NML)(t=0) .= NML_init
        tx_pwrs(:NLK)(t=0) .= NLK_init
    end

    if :rx in statetypes
        @unpack precondition_rx, shuffle_rx = parameters()
        rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size,ntimes+1), path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
        rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, ens_size)) #with only 4 possible values, we initialize with a uniform distribution of [0,3]
    
    end
   
    R = Float64[]
    if :amp in datatypes
        R = [R; fill(σamp^2, npaths)]
    end
    if :phase in datatypes
        R = [R; fill(σphase^2, npaths)]
    end
    #R = [fill(σamp^2, npaths); fill(σphase^2, npaths)] # Vector is transformed into a diagonal matrix in MVIA

    # Run model
    ym = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, ntimes+1); 
        field=[:amp, :phase], path=pathname.(paths), ens=1:ens_size, t=0:ntimes)
    
    if precondition_rx
        @unpack precon_ens_size = parameters()
        start_time = Dates.now()
        @info "Preconditioning RX offsets: " start_time
        ym_precondition = deepcopy(ym(t=0))
        basic_foward_model = (z -> model(itp, z.xy_state, paths, dt; pathstep))
        prior_state = (; xy_state = xy_state(t=0))
        ensemble_model!(ym_precondition, basic_foward_model, prior_state)

        ϕ_statistics = KeyedArray(fill(NaN,4,npaths,ens_size,ntimes+1), ϕ_off=0:3, path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
        
        @showprogress Threads.@threads for e in ym_precondition.ens
            #outer set of ensembles, split to become new "forward models"

            #ens_size here is the inner set of ensembles used on each of the outer set.
            loc_rx_phi_offset = KeyedArray(fill(NaN,npaths,precon_ens_size,ntimes+1), path = pathname.(paths), ens=1:precon_ens_size, t=0:ntimes)
            loc_rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, precon_ens_size)) 
            #Pulls a new set of samples for each outer ensemble member, IE the inner ensembles should be different across outer ensemble members
            loc_ym_precondition = deepcopy(ym)

            for t in ym.t
                loc_ym_precondition(t=t) .= ym_precondition
            end

            for i in 1:maximum(loc_ym_precondition.t)
                
                ### G(b)
                for ee in loc_rx_phi_offset.ens
                    loc_ym_precondition(field=:phase, ens=ee, t=i) .+= loc_rx_phi_offset(ens=ee,t=i-1) .* (π/2)
                end

                ybar = mean(loc_ym_precondition(t=i), dims=:ens)

                Y = similar(loc_ym_precondition(t=i))
                Y(:amp) .= loc_ym_precondition(t=i,field=:amp) .- ybar(:amp)
                Y(:phase) .= phasediff.(loc_ym_precondition(t=i,field=:phase), ybar(:phase))
                #Currently, MVIA.rx_phi_offset() assumes amplitude and phase data. 

                #Y = MVIA.phasediff.(loc_ym_precondition(:phase), ybar(:phase)) #Exclude amplitude in constructing Y matrix
                xnew_phi = MVIA.rx_phi_update(loc_rx_phi_offset(t=i-1), data(t=i), ybar, Y, R; ρ=ρ)
                loc_rx_phi_offset(t=i) .= mod.(round.(xnew_phi), 4)
                for p in ϕ_statistics.path

                    vals = parent(loc_rx_phi_offset(t=i, path=p))
                    N = length(vals)

                    ratios = [count(==(off), vals) / N for off in ϕ_statistics.ϕ_off]

                    ϕ_statistics(path=p, ens=e, t=i) .= ratios
                end
            end

        end
        final_ϕ_statistics = dropdims(mean(ϕ_statistics(t=maximum(ϕ_statistics.t)), dims=:ens),dims=:ens)
        jldsave(joinpath(resdir(scenario), "$(scenario)_preconOffsets.jld2"); ϕ_statistics, final_ϕ_statistics, ym_precondition)

        offs = collect(final_ϕ_statistics.ϕ_off)  # e.g. 0:3

        rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size,ntimes+1), path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
        #^shouldn't be necessary now after renaming local variables. Need to confirm.
        for p in final_ϕ_statistics.path
            probs = Array(parent(parent(final_ϕ_statistics(path=p))))           
            dist  = Categorical(probs)

            rx_phi_offset(path=p, t=0) .= rand(dist, ens_size) .- 1
        end

        @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
    end 

    if :tx in statetypes && :xy in statetypes && !(:rx in statetypes)
        state = (; xy_state, tx_pwrs)
    elseif :rx in statetypes && :xy in statetypes && !(:tx in statetypes)
        state = (; xy_state, rx_phi_offset)
    elseif :rx in statetypes && :xy in statetypes && :tx in statetypes
        state = (; xy_state, tx_pwrs, rx_phi_offset)
    else 
        state = (; xy_state) 
    end

    forward_model =
        (:tx in statetypes && :xy in statetypes) ?
            (z -> model(itp, z.xy_state, rebuildpaths(paths, z.tx_pwrs), dt; pathstep)) :
            (z -> model(itp, z.xy_state, paths, dt; pathstep))

    H!(x, t) = ensemble_model!(ym(t=t), forward_model, x)


    for i in 1:ntimes
        start_time = Dates.now()
        @info "Iteration" i=i start_time

        if :tx in statetypes
            if shuffle_tx #We shuffle at input to forward model rather than output so saved outputs are correctly associated
                @info "Shuffling TX Elements of State Vector"
                shuffled = deepcopy(state.tx_pwrs(t=i-1))
                for tx in shuffled.pwrs 
                    shuffled(pwrs = tx) .= shuffle(shuffled(pwrs = tx))
                end
                tx_pwrs = shuffled
            else
                tx_pwrs = deepcopy(state.tx_pwrs(t=i-1))
            end
        end

        if shuffle_xy
            @info "Shuffling XY Grid Elements of State Vector"
            shuffled_xy = deepcopy(state.xy_state(t=i-1))
            for x in shuffled_xy.x
                for y in shuffled_xy.y
                    for f in shuffled_xy.field
                        var = std(shuffled_xy(x=x,y=y,field=f))
                        if (f == :h && var < 1.5) || (f == :b && var < 0.024) #TODO add counter for percentage of points shuffled
                            shuffled_xy(x=x,y=y,field=f) .= shuffle(shuffled_xy(x=x,y=y,field=f))
                        end
                    end
                end
            end
            xy_state = shuffled_xy
        else
            xy_state = state.xy_state(t=i-1)
        end

        if :rx in statetypes
            if shuffle_rx
                @info "Shuffling RX Elements of State Vector"
                shuffled = deepcopy(state.rx_phi_offset(t=i-1))
                for p in shuffled.path
                    shuffled(path = p) .= shuffle(shuffled(path = p))
                end
                rx_phi_offset = shuffled
            else
                rx_phi_offset = state.rx_phi_offset(t=i-1)
            end
        end
        if :tx in statetypes && :rx in statetypes
            xold = (; xy_state, tx_pwrs, rx_phi_offset)
        elseif :tx in statetypes
            xold = (; xy_state, tx_pwrs)
        elseif :rx in statetypes
            xold = (; xy_state, rx_phi_offset)
        elseif :xy in statetypes
            xold = (; xy_state)
        else
            error("Incomplete state types for LETKF")
        end

        xnew = LETKF_measupdate(x->H!(x,i-1), xold, data(t=i), R; ρ=ρ,
            localization=localization, datatypes=datatypes)

        # Constrain beta and h' to physically realistic values
        #TODO This performs many findfirst() computations. Rewrite for computational efficiency 
        for e in xnew.xy_state.ens
            l = findfirst(==(e),  axiskeys(xnew.xy_state, :ens))

            for x in xnew.xy_state.x
                k = findfirst(==(x),  axiskeys(xnew.xy_state, :x))

                for y in xnew.xy_state.y
                    j = findfirst(==(y),  axiskeys(xnew.xy_state, :y))

                    h = xnew.xy_state(x=x, y=y, ens=e, field=:h)
                    b = xnew.xy_state(x=x, y=y, ens=e, field=:b)

                    if h < MIN_H
                        idx = findfirst(==(:h), axiskeys(xnew.xy_state, :field))
                        if h < (MIN_H - (2 * median(hB)))
                            # If more than 2 sigma of the initial distribution beyond the physical bounds imposed, resample at random from the initial distribution
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = rand(state.xy_state(t=0, x=x, y=y, field=:h))
                        else #otherwise, snap to boundary value
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = MIN_H
                        end
                    elseif h > MAX_H
                        idx = findfirst(==(:h), axiskeys(xnew.xy_state, :field))
                        if h > (MAX_H + (2 * median(hB)))
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = rand(state.xy_state(t=0, x=x, y=y, field=:h))
                        else
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = MAX_H
                        end
                    end 

                    if b < MIN_BETA
                        idx = findfirst(==(:b), axiskeys(xnew.xy_state, :field))
                        if b < (MIN_BETA - (2 * median(bB)))
                            # If more than 2 sigma of the initial distribution beyond the physical bounds imposed, resample at random from the initial distribution
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = rand(state.xy_state(t=0, x=x, y=y, field=:b))
                        else
                            # Otherwise, snap to boundary value
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = MIN_BETA
                        end
                    elseif b > MAX_BETA
                        idx = findfirst(==(:b), axiskeys(xnew.xy_state, :field))
                        if b > (MAX_BETA + (2 * median(bB)))
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = rand(state.xy_state(t=0, x=x, y=y, field=:b))
                        else
                            @inbounds parent(xnew.xy_state)[idx, j, k, l] = MAX_BETA
                        end
                    end 

                end
            end
        end

        #=xnew.xy_state(:b)[xnew.xy_state(:b) .< MIN_BETA] .= MIN_BETA
        xnew.xy_state(:b)[xnew.xy_state(:b) .> MAX_BETA] .= MAX_BETA
        xnew.xy_state(:h)[xnew.xy_state(:h) .< MIN_H] .= MIN_H
        xnew.xy_state(:h)[xnew.xy_state(:h) .> MAX_H] .= MAX_H =#

        state.xy_state(t=i) .= xnew.xy_state

        if :tx in statetypes #TODO perhaps change this to pass TX_range to parameters rather than assign 4 consts.
            xnew.tx_pwrs(:NLK)[xnew.tx_pwrs(:NLK) .< NLK_LOWER] .= NLK_LOWER
            xnew.tx_pwrs(:NLK)[xnew.tx_pwrs(:NLK) .> NLK_UPPER] .= NLK_UPPER
            xnew.tx_pwrs(:NML)[xnew.tx_pwrs(:NML) .< NML_LOWER] .= NML_LOWER
            xnew.tx_pwrs(:NML)[xnew.tx_pwrs(:NML) .> NML_UPPER] .= NML_UPPER
            state.tx_pwrs(t=i) .= xnew.tx_pwrs
        end

        # Round to nearest int and modulus to [0, 3]
        if :rx in statetypes
            state.rx_phi_offset(t=i) .= mod.(round.(xnew.rx_phi_offset), 4)
        end

        @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
        jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)
    end

    # Compute ym with final estimate
    xy_state = state.xy_state(t=ntimes)
    if :tx in statetypes
        tx_pwrs = state.tx_pwrs(t=ntimes)
    end
    if :rx in statetypes
        rx_phi_offset = state.rx_phi_offset(t=ntimes)
    end
    if :tx in statetypes && :rx in statetypes
        xfinal = (; xy_state, tx_pwrs, rx_phi_offset)
    elseif :tx in statetypes
        xfinal = (; xy_state, tx_pwrs)
    elseif :rx in statetypes
        xfinal = (; xy_state, rx_phi_offset)
    else
        xfinal = (; xy_state)
    end

    H!(xfinal, ntimes)

    jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)

    return state, data, ym
end
