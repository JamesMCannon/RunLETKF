function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, ρ,
        localization, statetypes, xy_file = parameters()
    @unpack modelproj = common_simulation()

    shuffle_xy = get(parameters(),:shuffle_xy, false)

    lengthscale = only(modelsteps).lengthscale

    paths = buildtruepaths()
    npaths = length(paths)
    
    xy_state = KeyedArray(fill(NaN, 2, length(y_grid), length(x_grid), ens_size, ntimes+1),
        field=[:h, :b], y=y_grid, x=x_grid, ens=1:ens_size, t=0:ntimes)

    if xy_file == "false"
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

        xy_state(:h)(t=0) .= h_init
        xy_state(:b)(t=0) .= b_init
    else
        background_state = load(xy_file,"state")
        if background_state.xy_state.ens != 1:ens_size
            error("Ensemble size mismatch")
        end
        init_t = first_t_with_agreement(background_state.rx_phi_offset, 0.85)
        xy_state(t=0) .= background_state.xy_state(t=init_t+3)
        @info "Stacked state vector starting with $(init_t+3)th iteration from file"
    end

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
    if :rx in statetypes
        if precondition_rx
            @unpack precon_ens_size, precon_itrs = parameters()
            if precon_itrs > length(data.t)
                error("More preconditioning iterations specified than data time steps available")
            end
            start_time = Dates.now()
            @info "Preconditioning RX offsets: " start_time
            ym_precondition = deepcopy(ym(t=0))
            basic_forward_model = (z -> model(itp, z.xy_state, paths, dt; pathstep))
            prior_state = (; xy_state = xy_state(t=0))
            ensemble_model!(ym_precondition, basic_forward_model, prior_state)

            ϕ_statistics = KeyedArray(fill(NaN,4,npaths,ens_size,precon_itrs+1), ϕ_off=0:3, path = pathname.(paths), ens=1:ens_size, t=0:precon_itrs)
            
            @showprogress Threads.@threads for e in ym_precondition.ens
                #outer set of ensembles, split to become new "forward models"

                #ens_size here is the inner set of ensembles used on each of the outer set.
                loc_rx_phi_offset = KeyedArray(fill(NaN,npaths,precon_ens_size,precon_itrs+1), path = pathname.(paths), ens=1:precon_ens_size, t=0:precon_itrs)
                loc_rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, precon_ens_size)) 
                #Pulls a new set of samples for each outer ensemble member, IE the inner ensembles should be different across outer ensemble members
                loc_ym_precondition = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, precon_itrs+1); 
                    field=[:amp, :phase], path=pathname.(paths), ens=1:ens_size, t=0:precon_itrs)

                for t in loc_ym_precondition.t
                    for ee in loc_ym_precondition.ens
                        #Set each of the inner ensemble members to the current out ensemble member's ym. Each of these will be modified by the local RX offset, ie, G(b)
                        loc_ym_precondition(t=t, ens=ee) .= ym_precondition(ens=e)
                    end
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

            for p in final_ϕ_statistics.path
                probs = Array(parent(parent(final_ϕ_statistics(path=p))))           
                dist  = Categorical(probs)

                rx_phi_offset(path=p, t=0) .= rand(dist, ens_size) .- 1
            end

            @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
        end 
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


function rundualletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, ρ,
        localization, statetypes, xy_file = parameters()
    @unpack modelproj = common_simulation()

    shuffle_xy = get(parameters(),:shuffle_xy, false)

    lengthscale = only(modelsteps).lengthscale

    paths = buildtruepaths()
    npaths = length(paths)
    pathnames = pathname.(paths)

    xy_state = KeyedArray(fill(NaN, 2, length(y_grid), length(x_grid), ens_size, ntimes+1),
        field=[:h, :b], y=y_grid, x=x_grid, ens=1:ens_size, t=0:ntimes)

    if xy_file == "false"
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

        xy_state(:h)(t=0) .= h_init
        xy_state(:b)(t=0) .= b_init
    else
        background_state = load(xy_file,"state")
        if background_state.xy_state.ens != 1:ens_size
            error("Ensemble size mismatch")
        end
        init_t = first_t_with_agreement(background_state.rx_phi_offset, 0.85)
        xy_state(t=0) .= background_state.xy_state(t=init_t+3)
        @info "Starting with $(init_t+3)th iteration from file"
    end

    state = (; xy_state)

    if :tx in statetypes || :rx in statetypes
        @unpack precon_ens_size, precon_itrs = parameters()
            if precon_itrs > length(data.t)
                error("More dual filter iterations specified than data time steps available")
            end
    end

    if :tx in statetypes
        @unpack σNLK, σNML, shuffle_tx, NLKb, NMLb = parameters()

        NMLdistribution = Distributions.Normal(NMLb, σNML) #Nominal values from FWM dictionary
        NLKdistribution = Distributions.Normal(NLKb, σNLK)

        NML_init = rand(rng, NMLdistribution, ens_size)
        NLK_init = rand(rng, NLKdistribution, ens_size)

        tx_pwrs = KeyedArray(fill(NaN,2, ens_size, ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, t=0:ntimes)

        dual_tx_pwrs = KeyedArray(fill(NaN,2, ens_size, precon_ens_size, precon_itrs*ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, dual_ens=1:precon_ens_size, t=0:precon_itrs*ntimes)

        for e in dual_tx_pwrs.ens
            dual_tx_pwrs(pwrs=:NML, ens=e, t=0) .= rand(NMLdistribution, precon_ens_size)
            dual_tx_pwrs(pwrs=:NLK, ens=e, t=0) .= rand(NLKdistribution, precon_ens_size)
        end

        tx_pwrs(:NML)(t=0) .= NML_init
        tx_pwrs(:NLK)(t=0) .= NLK_init

        state = merge(state, (; tx_pwrs))
    end

    if :rx in statetypes
        @unpack precondition_rx, shuffle_rx = parameters()
        rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size,ntimes+1), path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
        rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, ens_size)) #with only 4 possible values, we initialize with a uniform distribution of [0,3]
        if precondition_rx #Assumed true if :rx in statetypes for dual LETKF. If statement kept in for future flexibility.
            @unpack precon_ens_size, precon_itrs = parameters()
            μ_ϕ_statistics = KeyedArray(fill(NaN,4,npaths,precon_itrs*ntimes+1), ϕ_off=0:3, path = pathname.(paths), t=0:precon_itrs*ntimes)
        end
        state = merge(state, (; rx_phi_offset))
    end
   
    R = Float64[]
    if :amp in datatypes
        R = [R; fill(σamp^2, npaths)]
    end
    if :phase in datatypes
        R = [R; fill(σphase^2, npaths)]
    end

    # Define model
    ym = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, ntimes+1); 
        field=[:amp, :phase], path=pathname.(paths), ens=1:ens_size, t=0:ntimes)

    forward_model =
        (:tx in statetypes && :xy in statetypes) ?
            (z -> model(itp, z.xy_state, rebuildpaths(paths, z.tx_pwrs), dt; pathstep)) :
            (z -> model(itp, z.xy_state, paths, dt; pathstep))

    basic_forward_model = (z -> model(itp, z.xy_state, paths, dt; pathstep))

    H!(x, t) = ensemble_model!(ym(t=t), forward_model, x)


    for i in 1:ntimes
        start_time = Dates.now()
        @info "Outer Iteration" i=i start_time

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

        ensemble_model!(ym(t=i-1), basic_forward_model, (; xy_state)) #ym(t-1) updated here
        if :tx in statetypes || :rx in statetypes
            data_start_idx = window_start(length(data.t), i, precon_itrs)
        end

        if :tx in statetypes
            start_tx_time = Dates.now()
            @info "Filtering TX powers: " start_tx_time
            ym_precondition = deepcopy(ym(t=i-1))

            @showprogress Threads.@threads for e in ym_precondition.ens
                #outer set of ensembles, split to become new "forward models"
                #mutates dual_tx_pwrs with final estimate being stored at t = precon_itrs*ntimes
                inner_tx_filter!(ym_precondition(ens=e), dual_tx_pwrs(ens=e), data(t=data_start_idx:data_start_idx+precon_itrs), 
                    precon_ens_size, precon_itrs, R, ρ, i, paths=paths)
            end

            jldsave(joinpath(resdir(scenario), "$(scenario)_dual_TXOffsets.jld2"); dual_tx_pwrs)

            for tx in dual_tx_pwrs.pwrs
                state.tx_pwrs(pwrs=tx, t=i) .= sample_values(dual_tx_pwrs(pwrs=tx, t=precon_itrs*i), ens_size)
            end

            if :tx in statetypes #TODO perhaps change this to pass TX_range to parameters rather than assign 4 consts.
                state.tx_pwrs(:NLK)(t=i)[state.tx_pwrs(pwrs = :NLK, t=i) .< NLK_LOWER] .= NLK_LOWER
                state.tx_pwrs(:NLK)(t=i)[state.tx_pwrs(pwrs = :NLK, t=i) .> NLK_UPPER] .= NLK_UPPER
                state.tx_pwrs(:NML)(t=i)[state.tx_pwrs(pwrs = :NML, t=i) .< NML_LOWER] .= NML_LOWER
                state.tx_pwrs(:NML)(t=i)[state.tx_pwrs(pwrs = :NML, t=i) .> NML_UPPER] .= NML_UPPER
            end
            
            ## G(b): apply TX power offsets
            for e in state.tx_pwrs.ens
                for tx in state.tx_pwrs.pwrs
                    txpaths = pathnames[startswith.(pathnames, String(tx) * "-")]

                    Δpwr_log = log10(state.tx_pwrs(pwrs=tx, ens=e, t=i) / txpower(paths, String(tx)))
                    ym(field=:amp, ens=e, t=i-1, path=txpaths) .+= Δpwr_log * 10  #10 dB per decade
                end
            end

            @info "Elapsed" Δt=canonicalize(Dates.now() - start_tx_time)

        end


        if :rx in statetypes
            #No additional checks necessary; if :rx is in statetypes and this function is called, precondition_rx must be true
            start_rx_time = Dates.now()
            @info "Filtering RX offsets: " start_rx_time
    
            ym_precondition = deepcopy(ym(t=i-1))            

            ϕ_statistics = KeyedArray(fill(NaN,4,npaths,ens_size,precon_itrs+1), ϕ_off=0:3, path = pathname.(paths), ens=1:ens_size, t=0:precon_itrs)

            initial_ϕ_statistics = KeyedArray(fill(NaN,4,npaths), ϕ_off=0:3, path = pathname.(paths))

            for p in xold.rx_phi_offset.path
                    vals = parent(xold.rx_phi_offset(path=p))
                    N = length(vals)

                    ratios = [count(==(off), vals) / N for off in initial_ϕ_statistics.ϕ_off]

                    initial_ϕ_statistics(path=p) .= ratios
            end

            if i == 1
                μ_ϕ_statistics(t=0) .= initial_ϕ_statistics
            end

            for ee in ϕ_statistics.ens
                ϕ_statistics(ens=ee, t=0) .= initial_ϕ_statistics
            end
            #TODO simplify by removing initial_ϕ_statistics variable

            @showprogress Threads.@threads for e in ym_precondition.ens
                #outer set of ensembles, split to become new "forward models"
                inner_rx_filter!(ym_precondition(ens=e), initial_ϕ_statistics, ϕ_statistics(ens=e), data(t=data_start_idx:data_start_idx+precon_itrs), 
                    precon_ens_size, precon_itrs, R, ρ)
            end

            for tt in 1:precon_itrs
                μ_ϕ_statistics(t = ((i-1) * precon_itrs) + tt) .= dropdims(mean(ϕ_statistics(t=tt), dims=:ens),dims=:ens)
            end

            jldsave(joinpath(resdir(scenario), "$(scenario)_dual_RXOffsets.jld2"); μ_ϕ_statistics) #TODO save ym_precondition as well?
    
            for p in μ_ϕ_statistics.path
                probs = Array(parent(parent(μ_ϕ_statistics(path=p, t=(i*precon_itrs)))))           
                dist  = Categorical(probs)
                #Sample from the current final estimated distribution of offsets
                state.rx_phi_offset(path=p, t=i) .= rand(dist, ens_size) .- 1
            end

            ### G(b)
            # Apply final estimated RX offsets to ym prior to measurement update for then calculating the XY update
            for e in state.rx_phi_offset.ens
                ym(field=:phase, ens=e, t=i-1) .+= state.rx_phi_offset(ens=e,t=i) .* (π/2)
            end

            @info "Elapsed" Δt=canonicalize(Dates.now() - start_rx_time)

        end

        ybar = mean(ym(t=i-1), dims=:ens)
        yb = ym(t=i-1)
        if :amp in datatypes && :phase in datatypes
            Y = similar(yb)
            Y(:amp) .= yb(:amp) .- ybar(:amp)
            Y(:phase) .= MVIA.phasediff.(yb(:phase), ybar(:phase))
        elseif :amp in datatypes
            Y = yb(:amp) .- ybar(:amp)
        elseif :phase in datatypes
            Y = MVIA.phasediff.(yb(:phase), ybar(:phase))
        else
            error("Unknown datatypes: $datatypes")
        end

        xnew_xy = MVIA.xy_state_update(xold.xy_state, data(t=i), ybar, Y, R;
            ρ=ρ, localization=localization, datatypes=datatypes)

        xnew = (; xy_state=xnew_xy) 

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



        @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
        jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)
    end

    #=
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
    =#
    
    jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)

    return state, data, ym
end



function inner_rx_filter!(ym_precondition, initial_ϕ_statistics, ϕ_statistics, data, 
    precon_ens_size, precon_itrs, R, ρ) 
    #TODO Refactor ϕ_statistics such that initial_ϕ_statistics is not needed as a separate variable
    #TODO Rename variables away from "precondition" since this is now part of dual filter rather than preconditioning step
    
    paths = buildpaths()
    npaths = length(paths)

    # Allocate local RX offsets (inner ensembles)
    loc_rx_phi_offset = KeyedArray(fill(NaN, npaths, precon_ens_size, precon_itrs + 1),
        path = pathname.(paths),
        ens  = 1:precon_ens_size,
        t    = 0:precon_itrs
    )

    # Sample initial offsets for inner ensembles
    for p in initial_ϕ_statistics.path
        probs = Array(parent(parent(initial_ϕ_statistics(path=p))))
        dist  = Categorical(probs)
        loc_rx_phi_offset(path=p, t=0) .= rand(dist, precon_ens_size) .- 1
    end

    # Local preconditioned forward model
    loc_ym_precondition = KeyedArray(
        Array{Float64,4}(undef, 2, npaths, precon_ens_size, precon_itrs + 1),
        field = [:amp, :phase],
        path  = pathname.(paths),
        ens   = 1:precon_ens_size,
        t     = 0:precon_itrs
    )

    # Initialize all inner ensembles from outer ensemble member `e`
    for t in loc_ym_precondition.t
        for ee in loc_ym_precondition.ens
            loc_ym_precondition(t=t, ens=ee) .= ym_precondition
        end
    end

    # Dual filter iterations
    for dual_i in 1:maximum(loc_ym_precondition.t)

        ## G(b): apply RX phase offsets
        for ee in loc_rx_phi_offset.ens
            loc_ym_precondition(field=:phase, ens=ee, t=dual_i-1) .+=
                loc_rx_phi_offset(ens=ee, t=dual_i-1) .* (π/2)
        end

        ybar = mean(loc_ym_precondition(t=dual_i-1), dims=:ens)

        Y = similar(loc_ym_precondition(t=dual_i-1))
        Y(:amp)   .= loc_ym_precondition(t=dual_i-1, field=:amp)   .- ybar(:amp)
        Y(:phase) .= phasediff.(loc_ym_precondition(t=dual_i-1, field=:phase),
            ybar(:phase))

        #Currently, MVIA.rx_phi_offset() assumes amplitude and phase data. 
        #Y = MVIA.phasediff.(loc_ym_precondition(:phase), ybar(:phase)) #Exclude amplitude in constructing Y matrix

        tkey = data.t[dual_i] #assumes windowed data of size at least precon_itrs
        xnew_phi = MVIA.rx_phi_update(loc_rx_phi_offset(t=dual_i-1),
            data(t=tkey), ybar, Y, R; ρ = ρ)

        loc_rx_phi_offset(t=dual_i) .= mod.(round.(xnew_phi), 4)

        # Update statistics
        for p in ϕ_statistics.path
            vals = parent(loc_rx_phi_offset(t=dual_i, path=p))
            N = length(vals)
            ratios = [count(==(off), vals) / N for off in ϕ_statistics.ϕ_off]
            ϕ_statistics(path=p, t=dual_i) .= ratios
        end
    end

    return nothing
end


function inner_tx_filter!(ym_precondition, dual_tx_pwrs, data, 
    precon_ens_size, precon_itrs, R, ρ, i; paths=buildtruepaths()) 
    #dual_tx_pwrs: KeyedArray of size (2, precon_ens_size, precon_itrs+1) with pwrs=:NML,:NLK
    
    npaths = length(paths)

    # Local preconditioned forward model
    loc_ym_precondition = KeyedArray(
        Array{Float64,4}(undef, 2, npaths, precon_ens_size, precon_itrs + 1),
        field = [:amp, :phase],
        path  = pathname.(paths),
        ens   = 1:precon_ens_size,
        t     = 0:precon_itrs
    )

    # Initialize all inner ensembles from outer ensemble member `e`
    for t in loc_ym_precondition.t
        for ee in loc_ym_precondition.ens
            loc_ym_precondition(t=t, ens=ee) .= ym_precondition
        end
    end

    loc_tx_pwrs = KeyedArray(fill(NaN,2, precon_ens_size, precon_itrs+1), pwrs = [:NML, :NLK], ens=1:precon_ens_size, t=0:precon_itrs)
    #needed because tx_pwrs_update requires structure with dimesion ens, not dual_ens. Currently keeping both so dual_tx_pwrs can be mutated for all ensembles.
    pathnames = pathname.(paths)

    # Dual filter iterations
    for dual_i in 1:maximum(loc_ym_precondition.t)
        dual_tx_t_idx = (i-1)*precon_itrs+dual_i
        ## G(b): apply RX phase offsets
        for ee in dual_tx_pwrs.dual_ens
            for tx in dual_tx_pwrs.pwrs
                txpaths = pathnames[startswith.(pathnames, String(tx) * "-")]

                Δpwr_log = log10(dual_tx_pwrs(pwrs=tx, dual_ens=ee, t=dual_tx_t_idx-1) / txpower(paths, String(tx)))
                loc_ym_precondition(field=:amp, ens=ee, t=dual_i-1, path=txpaths) .+= Δpwr_log * 10  #10 dB per decade

            end
        end

        ybar = mean(loc_ym_precondition(t=dual_i-1), dims=:ens)

        Y = similar(loc_ym_precondition(t=dual_i-1))
        Y(:amp)   .= loc_ym_precondition(t=dual_i-1, field=:amp)   .- ybar(:amp)
        Y(:phase) .= phasediff.(loc_ym_precondition(t=dual_i-1, field=:phase),
            ybar(:phase))

        #Currently, MVIA.rx_phi_offset() assumes amplitude and phase data. 

        for tx in dual_tx_pwrs.pwrs
            loc_tx_pwrs(pwrs=tx, t=dual_i-1) .= log10.(strip(dual_tx_pwrs(pwrs=tx, t=dual_tx_t_idx-1)))
        end

        tkey = data.t[dual_i] #assumes windowed data of size at least precon_itrs
        xnew_amp = MVIA.tx_pwrs_update(loc_tx_pwrs(t=dual_i-1), data(t=tkey), ybar, Y, R; ρ = ρ)

        for tx in dual_tx_pwrs.pwrs
            dual_tx_pwrs(pwrs=tx, t=dual_tx_t_idx) .= 10 .^(strip(xnew_amp(pwrs=tx)))
        end
    end

    return nothing
end