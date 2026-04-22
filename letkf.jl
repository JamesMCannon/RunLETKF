function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, data, itp, ρ,
        localization, statetypes, xy_file, paths, R = parameters()
    @unpack modelproj = common_simulation()

    shuffle_xy = get(parameters(),:shuffle_xy, false)

    lengthscale = only(modelsteps).lengthscale

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
        if :rx in statetypes
            init_t = first_t_with_agreement(background_state.rx_phi_offset, 0.85) + 3
            @info "Stacked state vector starting with $(init_t+3)th iteration from file"
        else
            init_t = 10
            @info "Starting with final t of xy_file $(xy_file)"
        end
        xy_state(t=0) .= background_state.xy_state(t=init_t)
        
    end

    if :tx in statetypes || :rx in statetypes
        # split/dual filter parameters are required by the LETKF rx/tx paths but
        # not by the :categorical rx path. Default sensibly so users can run
        # `RX_METHOD=categorical DO_RX=true` without also setting `DO_DUAL=true`.
        rx_method_check = get(parameters(), :rx_method, :letkf)
        split_ens_size = get(parameters(), :split_ens_size, ens_size)
        split_itrs     = get(parameters(), :split_itrs, 1)
        filtertype     = get(parameters(), :filtertype, :stacked)
        if rx_method_check != :categorical && split_itrs > length(data.t)
            error("More split filter iterations specified than data time steps available")
        end
    end

    if :tx in statetypes
        @unpack σNLK, σNML, shuffle_tx, NLKb, NMLb = parameters()

        NMLdistribution = Distributions.Normal(NMLb, σNML) #Nominal values from FWM dictionary
        NLKdistribution = Distributions.Normal(NLKb, σNLK)

        if filtertype != :split
            tx_pwrs = KeyedArray(fill(NaN,2, ens_size, ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, t=0:ntimes)

            tx_pwrs(:NML)(t=0) .= rand(rng, NMLdistribution, ens_size)
            tx_pwrs(:NLK)(t=0) .= rand(rng, NLKdistribution, ens_size)
        else
            tx_pwrs = KeyedArray(fill(NaN,2, ens_size, split_ens_size, ntimes+1), pwrs = [:NML, :NLK], ens=1:ens_size, split_ens=1:split_ens_size, t=0:ntimes)
            
            for e in tx_pwrs.ens
                tx_pwrs(pwrs=:NML, ens=e, t=0) .= rand(rng, NMLdistribution, split_ens_size)
                tx_pwrs(pwrs=:NLK, ens=e, t=0) .= rand(rng, NLKdistribution, split_ens_size)
            end

        end
        
    end

    if :rx in statetypes
        @unpack precondition_rx, shuffle_rx, rx_offsets = parameters()
        rx_method = get(parameters(), :rx_method, :letkf)

        if rx_method == :categorical
            # Categorical path: per-path log-posterior over k ∈ {0,1,2,3} accumulates
            # evidence across iterations. `rx_phi_offset` is kept (sampled per-member
            # from the posterior each iteration) so heatmaps and saves stay compatible.
            rx_phi_logpost = KeyedArray(
                fill(0.0, npaths, 4, ntimes+1);
                path = pathname.(paths), k = 0:3, t = 0:ntimes)
            # Uniform prior at t=0 (zeros = log of unnormalized uniform; normalize on read)
            rx_phi_offset = KeyedArray(
                fill(NaN, npaths, ens_size, ntimes+1);
                path = pathname.(paths), ens = 1:ens_size, t = 0:ntimes)
            rx_phi_offset(t=0) .= rand(rng, 0:3, npaths, ens_size)
            @info "rx_phi_offset using :categorical update path"
        elseif filtertype != :split
            rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size,ntimes+1), path = pathname.(paths), ens=1:ens_size, t=0:ntimes)
            rx_phi_offset(t=0) .= rand(rng, 0:3, npaths, ens_size) #with only 4 possible values, we initialize with a uniform distribution of [0,3]
        else
            rx_phi_offset = KeyedArray(fill(NaN,npaths,ens_size, split_ens_size, ntimes+1), path = pathname.(paths), ens=1:ens_size, split_ens=1:split_ens_size, t=0:ntimes)
            rx_phi_offset(t=0) .= rand(rng, 0:3, npaths, ens_size, split_ens_size)
            #rx_phi_offset(t=0) .= round.(rand(Distributions.Uniform(0,3), npaths, ens_size, split_ens_size))
            #=
            # Sample a circular guassian aka a Von Mises distribution at random giving ~10% to the furthest bin
            κ_rx = 0.25 #log(sqrt(10) - 1)  # ~10% probability at opposite value
            period = 4
            offsets = 0:3
            for e in 1:ens_size
                for (pi, p) in enumerate(pathname.(paths))
                    center = rand(rng, offsets)
                    θ_center = 2π * center / period
                    w = [exp(κ_rx * cos(2π * v / period - θ_center)) for v in offsets]
                    rx_phi_offset(path=p, ens=e, t=0) .= sample(rng, collect(offsets), Weights(w), split_ens_size)
                end
            end
            =#
        end 
        jldsave(joinpath(resdir(scenario), "rx_offsets_$scenario.jld2"); rx_offsets)

        #TODO add option to 'precondition' rx offsets by running a loop split_itrs times without updating the xy_state
    end
   
    ym = KeyedArray(Array{Float64,4}(undef, 2, npaths, ens_size, ntimes+1); 
        field=[:amp, :phase], path=pathname.(paths), ens=1:ens_size, t=0:ntimes)

    state = (; xy_state)
    if :tx in statetypes
        state = merge(state, (; tx_pwrs))
    end
    if :rx in statetypes
        state = merge(state, (; rx_phi_offset))
        if get(parameters(), :rx_method, :letkf) == :categorical
            state = merge(state, (; rx_phi_logpost))
        end
    end

    # Intermediate bias history — only for windowed split filter (split_itrs > 1).
    #   tx_pwrs_history       (pwrs, ens, split_ens, t=1:ntimes, wstep=1:split_itrs)
    #   rx_phi_offset_history (path, ens, split_ens, t=1:ntimes, wstep=1:split_itrs)
    #   window_obs            Matrix{Int}(ntimes, split_itrs)  — data time index j
    #                         consumed at each (main iteration i, window step k)
    if (:tx in statetypes || :rx in statetypes) && filtertype == :split && split_itrs > 1
        if :tx in statetypes
            tx_pwrs_history = KeyedArray(
                fill(NaN, length(state.tx_pwrs.pwrs), ens_size, split_ens_size,
                     ntimes, split_itrs);
                pwrs      = state.tx_pwrs.pwrs,
                ens       = 1:ens_size,
                split_ens = 1:split_ens_size,
                t         = 1:ntimes,
                wstep     = 1:split_itrs)
            state = merge(state, (; tx_pwrs_history))
        end
        if :rx in statetypes
            rx_phi_offset_history = KeyedArray(
                fill(NaN, npaths, ens_size, split_ens_size, ntimes, split_itrs);
                path      = state.rx_phi_offset.path,
                ens       = 1:ens_size,
                split_ens = 1:split_ens_size,
                t         = 1:ntimes,
                wstep     = 1:split_itrs)
            state = merge(state, (; rx_phi_offset_history))
        end
        window_obs = zeros(Int, ntimes, split_itrs)
        state = merge(state, (; window_obs))
    end

    # TODO if tx_pwrs contains dims :split_ens, rebuildpaths should take 
    # dropdims(mean(z.tx_pwrs, dims=:split_ens), dims=:split_ens) as input rather 
    # than z.tx_pwrs itself.
    if :tx in statetypes && :xy in statetypes
         @info "Using forward model with TX power offsets"
         if filtertype == :split
            forward_model = (z -> model(itp, z.xy_state, rebuildpaths(paths, dropdims(mean(z.tx_pwrs, dims=:split_ens), dims=:split_ens)), dt; pathstep))
         else
            forward_model = (z -> model(itp, z.xy_state, rebuildpaths(paths, z.tx_pwrs), dt; pathstep))
         end
    else
        @info "Using forward model without TX power offsets"
        forward_model = (z -> model(itp, z.xy_state, paths, dt; pathstep))
    end

    H!(x, t) = ensemble_model!(ym(t=t), forward_model, x)

    # Categorical RX path uses a different H! that samples per-member k from the
    # current log-posterior. Closure captures `rng` and the live `state.rx_phi_logpost`.
    rx_method = (:rx in statetypes) ? get(parameters(), :rx_method, :letkf) : :letkf
    rx_commit_threshold = (rx_method == :categorical) ?
        get(parameters(), :rx_commit_threshold, 1.0) : 1.0
    H_cat!(x, t, log_post_t) = ensemble_model!(ym(t=t), forward_model, x, log_post_t, rng;
                                               commit_threshold=rx_commit_threshold)

    for i in 1:ntimes
        start_time = Dates.now()
        @info "Iteration" i=i start_time

        # ── Categorical RX branch ─────────────────────────────────────────────
        # Self-contained: forward model with per-member sampled k, accumulate
        # log-posterior, then do an xy_state-only LETKF update on the resulting yb.
        # Skips all LETKF rx_phi_offset machinery (tie-breaking, split paths, etc.).
        if rx_method == :categorical
            # Carry posterior forward (will be incremented in-place by accumulation)
            state.rx_phi_logpost(t=i) .= state.rx_phi_logpost(t=i-1)

            xy_state_in = shuffle_xy ? begin
                @info "Shuffling XY Grid Elements of State Vector"
                shuffled_xy = deepcopy(state.xy_state(t=i-1))
                for x in shuffled_xy.x, y in shuffled_xy.y, f in shuffled_xy.field
                    var = std(shuffled_xy(x=x,y=y,field=f))
                    if (f == :h && var < 1.5) || (f == :b && var < 0.024)
                        shuffled_xy(x=x,y=y,field=f) .= shuffle(shuffled_xy(x=x,y=y,field=f))
                    end
                end
                shuffled_xy
            end : state.xy_state(t=i-1)

            # Local buffer for this iteration's per-member sampled offsets.
            # Avoids relying on AxisKeys view-vs-copy semantics in `state.rx_phi_offset(t=i)`.
            rx_buf = KeyedArray(
                fill(NaN, npaths, ens_size);
                path = state.rx_phi_offset.path, ens = state.rx_phi_offset.ens)

            xold = (; xy_state = xy_state_in, rx_phi_offset = rx_buf)
            if :tx in statetypes
                tx_pwrs_in = deepcopy(state.tx_pwrs(t=i-1))
                if shuffle_tx
                    for tx in tx_pwrs_in.pwrs
                        tx_pwrs_in(pwrs=tx) .= shuffle(tx_pwrs_in(pwrs=tx))
                    end
                end
                xold = merge(xold, (; tx_pwrs = tx_pwrs_in))
            end

            # Run forward model with per-member sampled offsets. H_cat! writes the
            # drawn offsets into xold.rx_phi_offset (== rx_buf) and returns yb.
            yb = H_cat!(xold, i-1, state.rx_phi_logpost(t=i))

            # Persist what was actually used at this iteration
            state.rx_phi_offset(t=i) .= rx_buf

            # Accumulate evidence into log-posterior at t=i
            MVIA.categorical_rx_update!(state.rx_phi_logpost(t=i), yb,
                                        rx_buf, data(t=i), R)

            # xy_state LETKF update on the offset-corrected yb (no extra forward call)
            xnew_xy = MVIA.xy_only_update(yb, xold.xy_state, data(t=i), R;
                ρ=ρ, localization=localization, datatypes=datatypes)

            xnew_xy(:b)[xnew_xy(:b) .< MIN_BETA] .= MIN_BETA
            xnew_xy(:b)[xnew_xy(:b) .> MAX_BETA] .= MAX_BETA
            xnew_xy(:h)[xnew_xy(:h) .< MIN_H]    .= MIN_H
            xnew_xy(:h)[xnew_xy(:h) .> MAX_H]    .= MAX_H
            state.xy_state(t=i) .= xnew_xy

            if :tx in statetypes
                # tx_pwrs not updated in categorical-rx-only branch; carry forward.
                state.tx_pwrs(t=i) .= xold.tx_pwrs
            end

            @info "Elapsed" Δt=canonicalize(Dates.now() - start_time)
            jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)
            continue
        end
        # ── End categorical branch; standard LETKF path follows ───────────────

        #TODO adapt shuffle to handle split_ens and move to dedicated function for readability
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

            
            # Guard: break mode ties in rx_phi_offset before update
            n_ties_broken = 0
            for p in rx_phi_offset.path
                for e in rx_phi_offset.ens
                    if filtertype == :split
                        vals = strip(rx_phi_offset(path=p, ens=e))  # 1D over split_ens
                    else
                        vals = strip(rx_phi_offset(path=p))          # 1D over ens
                    end

                    counts = [count(==(b), vals) for b in 0:3]
                    max_count = maximum(counts)
                    tied_values = findall(==(max_count), counts) .- 1  # back to 0-indexed {0,1,2,3}

                    if length(tied_values) > 1
                        winner, loser = tied_values[randperm(rng, length(tied_values))[1:2]]
                        loser_idx = findfirst(==(loser), vals)
                        if !isnothing(loser_idx)
                            if filtertype == :split
                                rx_phi_offset(path=p, ens=e)[loser_idx] = winner
                            else
                                rx_phi_offset(path=p)[loser_idx] = winner
                            end
                            n_ties_broken += 1
                        end
                    end

                    # For !split, the outer ens loop is redundant — break after first iteration
                    filtertype != :split && break
                end
            end
            if n_ties_broken > 0
                @info "Broke mode ties in rx_phi_offset" n_ties_broken filtertype
            end
        
        end

        xold = (; xy_state)
        if :tx in statetypes
            xold = merge(xold, (; tx_pwrs))
        end

        if :rx in statetypes
            xold = merge(xold, (; rx_phi_offset))
        end

        # ── Kalman update ─────────────────────────────────────────────────────────
        if :rx in statetypes || :tx in statetypes
            if filtertype == :split && split_itrs > 1 
                # Windowed pre-update: bias refined over split_itrs observations before
                # the single xy_state update.  After each bias_only_update! call, `yb`
                # reflects H(xb) under the current bias means; accumulated corrections
                # telescope to a single correction from prior to final estimate.

                # yb aliases ym(t=i-1); bias corrections accumulate in-place so that
                # ym(t=i-1) holds the bias-corrected predictions after the loop.
                yb = H!(xold, i-1)

                win_tx_pwrs       = haskey(xold, :tx_pwrs)       ? deepcopy(xold.tx_pwrs)       : nothing
                win_rx_phi_offset = haskey(xold, :rx_phi_offset) ? deepcopy(xold.rx_phi_offset) : nothing

                # TX clamped after each step so out-of-bounds values never feed forward
                # as a prior.  RX is already constrained to {0,1,2,3} in split_rx_update!.
                for (k, j) in enumerate(centered_window(i, split_itrs, ntimes))
                    @info "  Windowed bias pre-update" main_t=i wstep=k window_obs=j
                    new_tx, new_rx = MVIA.bias_only_update!(yb, win_tx_pwrs, win_rx_phi_offset,
                                                    data(t=j), R; ρ=ρ)
                    if !isnothing(new_tx)
                        new_tx(:NLK)[new_tx(:NLK) .< NLK_LOWER] .= NLK_LOWER
                        new_tx(:NLK)[new_tx(:NLK) .> NLK_UPPER] .= NLK_UPPER
                        new_tx(:NML)[new_tx(:NML) .< NML_LOWER] .= NML_LOWER
                        new_tx(:NML)[new_tx(:NML) .> NML_UPPER] .= NML_UPPER
                        win_tx_pwrs = new_tx
                        state.tx_pwrs_history(t=i, wstep=k) .= win_tx_pwrs
                    end
                    if !isnothing(new_rx)

                        # Guard: break mode ties in rx_phi_offset before update
                        n_ties_broken = 0
                        for p in new_rx.path
                            for e in new_rx.ens
                                if filtertype == :split
                                    vals = strip(new_rx(path=p, ens=e))  # 1D over split_ens
                                else
                                    vals = strip(new_rx(path=p))          # 1D over ens
                                end

                                counts = [count(==(b), vals) for b in 0:3]
                                max_count = maximum(counts)
                                tied_values = findall(==(max_count), counts) .- 1  # back to 0-indexed {0,1,2,3}

                                if length(tied_values) > 1
                                    winner, loser = tied_values[randperm(rng, length(tied_values))[1:2]]
                                    loser_idx = findfirst(==(loser), vals)
                                    if !isnothing(loser_idx)
                                        if filtertype == :split
                                            new_rx(path=p, ens=e)[loser_idx] = winner
                                        else
                                            new_rx(path=p)[loser_idx] = winner
                                        end
                                        n_ties_broken += 1
                                    end
                                end

                                # For !split, the outer ens loop is redundant — break after first iteration
                                filtertype != :split && break
                            end
                        end
                        if n_ties_broken > 0
                            @info "Broke mode ties in rx_phi_offset" n_ties_broken filtertype
                        end

                        win_rx_phi_offset = new_rx
                        state.rx_phi_offset_history(t=i, wstep=k) .= win_rx_phi_offset
                    end
                    state.window_obs[i, k] = j
                end

                xnew_xy = MVIA.xy_only_update(yb, xold.xy_state, data(t=i), R;
                    ρ=ρ, localization=localization, datatypes=datatypes)

                xnew = (; xy_state = xnew_xy)
                if !isnothing(win_tx_pwrs)
                    xnew = merge(xnew, (; tx_pwrs = win_tx_pwrs))
                end
                if !isnothing(win_rx_phi_offset)
                    xnew = merge(xnew, (; rx_phi_offset = win_rx_phi_offset))
                end

            else
                xnew = LETKF_measupdate(x->H!(x,i-1), xold, data(t=i), R; ρ=ρ,
                    localization=localization, datatypes=datatypes, filtertype=filtertype)
            end
        else
            xnew = LETKF_measupdate(x->H!(x,i-1), xold, data(t=i), R; ρ=ρ,
                localization=localization, datatypes=datatypes)
        end

        # Constrain beta and h' to physically realistic values
        xnew.xy_state(:b)[xnew.xy_state(:b) .< MIN_BETA] .= MIN_BETA
        xnew.xy_state(:b)[xnew.xy_state(:b) .> MAX_BETA] .= MAX_BETA
        xnew.xy_state(:h)[xnew.xy_state(:h) .< MIN_H] .= MIN_H
        xnew.xy_state(:h)[xnew.xy_state(:h) .> MAX_H] .= MAX_H 

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

    jldsave(joinpath(resdir(scenario), "$scenario.jld2"); state, data, ym)
    =#
    return state, data, ym
end