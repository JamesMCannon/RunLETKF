function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, σamp, σphase, σNML, σNLK, data, itp, 
        localization, statetypes = parameters()
    @unpack modelproj = common_simulation()

    ρ = get(parameters(), :ρ, 1.1) 

    shuffle_tx = get(parameters(),:shuffle_tx, false)
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
                                                     # possible corerction: b_init[b_init .< MIN_BETA] .= MIN_BETA
 
    locmask = anylocal(localization)
    h_init[CI[.!locmask],:] .= NaN
    b_init[CI[.!locmask],:] .= NaN

    xy_state = KeyedArray(fill(NaN, 2, length(y_grid), length(x_grid), ens_size, ntimes+1),
    field=[:h, :b], y=y_grid, x=x_grid, ens=1:ens_size, t=0:ntimes)

    xy_state(:h)(t=0) .= h_init
    xy_state(:b)(t=0) .= b_init

    if :tx in statetypes

        NLKb = get(parameters(), :NLKb, 250000)
        NMLb = get(parameters(), :NMLb, 233000)

        NMLdistribution = Distributions.Normal(NMLb, σNML) #Nominal values from FWM dictionary
        NLKdistribution = Distributions.Normal(NLKb, σNLK)

        NML_init = rand(rng, NMLdistribution, ens_size)
        NLK_init = rand(rng, NLKdistribution, ens_size)

        tx_pwrs = KeyedArray(fill(NaN,2, ens_size, ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, t=0:ntimes)

        tx_pwrs(:NML)(t=0) .= NML_init
        tx_pwrs(:NLK)(t=0) .= NLK_init
    end

    if :rx in statetypes

        rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size,ntimes+1), path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
        rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, ens_size)) #with only 4 possible values, we initialize with a uniform distribution of [0,3]

    end
   

    if :tx in statetypes && :xy in statetypes && !(:rx in statetypes)
        state = (; xy_state, tx_pwrs)
    elseif :rx in statetypes && :xy in statetypes && !(:tx in statetypes)
        state = (; xy_state, rx_phi_offset)
    elseif :rx in statetypes && :xy in statetypes && :tx in statetypes
        state = (; xy_state, tx_pwrs, rx_phi_offset)
    else 
        state = (; xy_state) #TODO refactor to (; xy_state) - requires functions in MVIA
    end

    R = [fill(σamp^2, npaths); fill(σphase^2, npaths)] # Vector is transformed into a diagonal matrix in MVIA

    # Run model
    ym = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, ntimes+1);
        field=[:amp, :phase], path=pathname.(paths), ens=state.xy_state.ens, t=0:ntimes)


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
            rx_phi_offset = state.rx_phi_offset(t=i-1)
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

        # Floor to β = 0.16
        xnew.xy_state(:b)[xnew.xy_state(:b) .< MIN_BETA] .= MIN_BETA
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