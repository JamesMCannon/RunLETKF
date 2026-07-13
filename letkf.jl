function runletkf(parameters)
    @unpack scenario, ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
        datatypes, h0, b0, hB, bB, rng, data, itp, ρ,
        localization, statetypes, xy_file, paths, R = parameters()
    @unpack modelproj = common_simulation()

    shuffle_xy = get(parameters(),:shuffle_xy, false)
    filtertype = get(parameters(), :filtertype, :stacked)

    lengthscale = only(modelsteps).lengthscale

    npaths = length(paths)
    nfields = length(datatypes)
    
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
        file_t = parse(Int, get(ENV, "FILE_T","-1"))
        background_state = load(xy_file,"state")
        @info "Loaded background state from file =  $xy_file"
        if background_state.xy_state.ens != 1:ens_size
            error("Ensemble size mismatch")
        end
        if :rx in statetypes && file_t == -1
            init_t = first_t_with_agreement(background_state.rx_phi_offset, 0.85) + 3
            @info "RX agreement found, starting with $(init_t)th iteration from file"
        elseif file_t != -1
            init_t = file_t
            @info "Starting with t of xy_file t = $(init_t)"
        else
            init_t = DATALENGTH  
        end

        src = background_state.xy_state(t=init_t)
    
        old_x_grid = axiskeys(src, :x)
        old_y_grid = axiskeys(src, :y)
    
        if old_x_grid == x_grid && old_y_grid == y_grid
            xy_state(t=0) .= src
            @info "Grids match; directly assigned xy_state from file at iteration $(init_t)"
        else
            @info "Grids don't match; interpolating new grid values from old grid"

            old_locmask = .!isnan.(vec(parent(src(:h)(ens=1))))
            old_itppts = MVIA.build_xygrid(old_locmask, old_x_grid, old_y_grid)
            old_itp = MVIA.ScatteredInterpolant(SI.GeneralizedPolyharmonic(1,1), modelproj, old_itppts)
    
            locmask = anylocal(localization)
            CI = CartesianIndices((length(y_grid), length(x_grid)))
    
            for e in 1:ens_size
                for field in (:h, :b)
                    vals = vec(parent(src(field)(ens=e)))
                    vgrid = MVIA.dense_grid(old_itp, vals, x_grid, y_grid)
                    xy_state(field)(t=0, ens=e) .= vgrid
                end
            end
    
            xy_state(t=0)[:, CI[.!locmask], :] .= NaN
    
            # Clamp to physical bounds after interpolation
            xy_state(:h)(t=0)[xy_state(:h)(t=0) .< MIN_H] .= MIN_H
            xy_state(:h)(t=0)[xy_state(:h)(t=0) .> MAX_H] .= MAX_H
            xy_state(:b)(t=0)[xy_state(:b)(t=0) .< MIN_BETA] .= MIN_BETA
            xy_state(:b)(t=0)[xy_state(:b)(t=0) .> MAX_BETA] .= MAX_BETA
    
            @info "Interpolated xy_state from old grid onto new grid; starting from iteration $(init_t)"
        end        
    end

    if :tx in statetypes || :rx in statetypes
        split_ens_size = get(parameters(), :split_ens_size, ens_size)
        split_itrs     = get(parameters(), :split_itrs, 1)
        if filtertype == :split && :tx in statetypes && split_itrs > length(data.t)
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
        @unpack rx_offsets, η = parameters()

        # Persistent per-path log-posterior over Bϕ ∈ {0,1,2,3}. Saved at every t —
        # this is the lossless statistic for downstream analysis.
        rx_phi_logpost = KeyedArray(
            fill(0.0, npaths, 4, ntimes+1);
            path = pathname.(paths), Bϕ = 0:3, t = 0:ntimes)

        # Per-iteration record of the offsets actually baked into yb for that
        # iteration's xy_state update. For :stacked this is the prior draw; for
        # :dual / :split it is the post-evidence draw.
        rx_phi_offset = KeyedArray(
            fill(NaN, npaths, ens_size, ntimes+1);
            path = pathname.(paths), ens = 1:ens_size, t = 0:ntimes)
        rx_phi_offset(t=0) .= rand(rng, 0:3, npaths, ens_size)

        jldsave(joinpath(resdir(scenario), "rx_offsets_$scenario.jld2"); rx_offsets)
    end
   
    # Ensemble prediction record: one row-block per observable field, in the
    # canonical datatypes order shared with data and R.
    ym = KeyedArray(Array{Float64,4}(undef, nfields, npaths, ens_size, ntimes+1); 
        field=collect(datatypes), path=pathname.(paths), ens=1:ens_size, t=0:ntimes)

    state = (; xy_state)
    if :tx in statetypes
        state = merge(state, (; tx_pwrs))
    end
    if :rx in statetypes
        state = merge(state, (; rx_phi_offset, rx_phi_logpost))
    end

    # Windowed split history — TX only. RX is categorical, accumulates evidence
    # across iterations natively in rx_phi_logpost, so it has no per-window record.
    if :tx in statetypes && filtertype == :split && split_itrs > 1
        tx_pwrs_history = KeyedArray(
            fill(NaN, length(state.tx_pwrs.pwrs), ens_size, split_ens_size,
                ntimes, split_itrs);
            pwrs = state.tx_pwrs.pwrs, ens = 1:ens_size, split_ens = 1:split_ens_size,
            t = 1:ntimes, wstep = 1:split_itrs)
        window_obs = zeros(Int, ntimes, split_itrs)
        state = merge(state, (; tx_pwrs_history, window_obs))
    end

    # ── Resume from existing checkpoint ───────────────────────────────────────
    savefile  = joinpath(resdir(scenario), "$scenario.jld2")
    start_iter = 1
    if isfile(savefile)
        saved       = load(savefile)
        saved_state = saved["state"]

        # bail out if the saved file doesn't match current parameters
        if keys(saved_state) != keys(state) ||
        size(saved_state.xy_state) != size(state.xy_state)
            error("Saved file for $scenario is incompatible with current parameters")
        end

        # highest iteration that actually has data written
        last_done = something(
            findlast(t -> !all(isnan, saved_state.xy_state(t=t)), 1:ntimes), 0)

        if last_done == ntimes
            @info "Scenario $scenario already complete; returning saved results"
            return saved_state, saved["data"], saved["ym"]
        elseif last_done > 0
            @info "Resuming $scenario from iteration $(last_done+1)/$ntimes"
            state = saved_state
            ym    = saved["ym"]
            start_iter = last_done + 1
        end
    end

    # TODO if tx_pwrs contains dims :split_ens, rebuildpaths should take 
    # dropdims(mean(z.tx_pwrs, dims=:split_ens), dims=:split_ens) as input rather 
    # than z.tx_pwrs itself.
    # The forward model returns a (field × path) KeyedArray of the requested
    # observables from a single LMP run per member (see model_observables).
    if :tx in statetypes && :xy in statetypes
         @info "Using forward model with TX power offsets"
         if filtertype == :split
            forward_model = (z -> model_observables(itp, z.xy_state, rebuildpaths(paths, dropdims(mean(z.tx_pwrs, dims=:split_ens), dims=:split_ens)), dt; pathstep, datatypes))
         else
            forward_model = (z -> model_observables(itp, z.xy_state, rebuildpaths(paths, z.tx_pwrs), dt; pathstep, datatypes))
         end
    else
        @info "Using forward model without TX power offsets"
        forward_model = (z -> model_observables(itp, z.xy_state, paths, dt; pathstep, datatypes))
    end

    H!(x, t) = ensemble_model!(ym(t=t), forward_model, x)

    # Defaults for the unified call to LETKF_measupdate; these are ignored if :rx or :tx aren't in statetypes
    rx_commit_threshold = (:rx in statetypes) ?
        get(parameters(), :rx_commit_threshold, 1.0) : 1.0
    η = (:rx in statetypes) ? get(parameters(), :η, 1.0) : 1.0
    
    for i in start_iter:ntimes
        start_time = Dates.now()
        @info "Iteration" i=i start_time

        # Epoch slice of the per-path, per-epoch observation-error variance,
        # flattened to the stacked field-major layout consumed by the updates.
        R_i = MVIA.stack_R(R(t=i), datatypes)

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

        # ── Assemble xold ────────────────────────────────────────────────────────
        xold = (; xy_state)
        if :tx in statetypes
            xold = merge(xold, (; tx_pwrs))
        end
        if :rx in statetypes
            # Carry posterior forward into the t=i slot. categorical_rx_update!
            # mutates this in place, so post-iteration state.rx_phi_logpost(t=i)
            # already reflects the evidence folded in.
            state.rx_phi_logpost(t=i) .= state.rx_phi_logpost(t=i-1)

            # Fresh per-iteration offset buffer; LETKF_*_update fills it with the
            # prior draw (stacked) or overwrites with the post-evidence draw (dual/split).
            rx_phi_offset_buf = KeyedArray(
                fill(NaN, npaths, ens_size);
                path = state.rx_phi_offset.path, ens = state.rx_phi_offset.ens)
            xold = merge(xold, (;
                rx_phi_offset  = rx_phi_offset_buf,
                rx_phi_logpost = state.rx_phi_logpost(t=i)))
        end

        # ── Measurement update ──────────────────────────────────────────────────
        if :tx in statetypes && filtertype == :split && split_itrs > 1
            # Windowed TX bias pre-update. RX (if present) gets a single categorical
            # update against data(t=i) — its posterior already accumulates across
            # iterations natively, so no windowing is needed.

            if :rx in statetypes
                MVIA._draw_prior_rx_offsets!(xold.rx_phi_offset, xold.rx_phi_logpost,
                    rng; commit_threshold=rx_commit_threshold)
            end

            yb = H!(xold, i-1)

            win_tx_pwrs = deepcopy(xold.tx_pwrs)
            for (k, j) in enumerate(centered_window(i, split_itrs, ntimes))
                @info "  Windowed TX bias pre-update" main_t=i wstep=k window_obs=j
                # Each window observation is weighted by its own epoch's R.
                new_tx = MVIA.bias_only_update!(yb, win_tx_pwrs, data(t=j),
                    MVIA.stack_R(R(t=j), datatypes); ρ=ρ, datatypes=datatypes)
                new_tx(:NLK)[new_tx(:NLK) .< NLK_LOWER] .= NLK_LOWER
                new_tx(:NLK)[new_tx(:NLK) .> NLK_UPPER] .= NLK_UPPER
                new_tx(:NML)[new_tx(:NML) .< NML_LOWER] .= NML_LOWER
                new_tx(:NML)[new_tx(:NML) .> NML_UPPER] .= NML_UPPER
                win_tx_pwrs = new_tx
                state.tx_pwrs_history(t=i, wstep=k) .= win_tx_pwrs
                state.window_obs[i, k] = j
            end

            if :rx in statetypes
                MVIA.categorical_rx_measupdate!(xold.rx_phi_logpost, yb,
                    xold.rx_phi_offset, data(t=i), R_i, rng;
                    η=η, commit_threshold=rx_commit_threshold, correct_yb=true,
                    datatypes=datatypes)
            end

            xnew_xy = MVIA.xy_only_update(yb, xold.xy_state, data(t=i), R_i;
                ρ=ρ, localization=localization, datatypes=datatypes)

            xnew = (; xy_state = xnew_xy, tx_pwrs = win_tx_pwrs)
            if :rx in statetypes
                xnew = merge(xnew, (; rx_phi_offset  = xold.rx_phi_offset, 
                    rx_phi_logpost = xold.rx_phi_logpost))
            end

        else 
            xnew = LETKF_measupdate(x->H!(x,i-1), xold, data(t=i), R_i; ρ=ρ,
                localization=localization, datatypes=datatypes, filtertype=filtertype,
                rng=rng, η=η, commit_threshold=rx_commit_threshold)
        end

        # ── Constrain and persist ───────────────────────────────────────────────
        xnew.xy_state(:b)[xnew.xy_state(:b) .< MIN_BETA] .= MIN_BETA
        xnew.xy_state(:b)[xnew.xy_state(:b) .> MAX_BETA] .= MAX_BETA
        xnew.xy_state(:h)[xnew.xy_state(:h) .< MIN_H]    .= MIN_H
        xnew.xy_state(:h)[xnew.xy_state(:h) .> MAX_H]    .= MAX_H
        state.xy_state(t=i) .= xnew.xy_state

        if :tx in statetypes
            xnew.tx_pwrs(:NLK)[xnew.tx_pwrs(:NLK) .< NLK_LOWER] .= NLK_LOWER
            xnew.tx_pwrs(:NLK)[xnew.tx_pwrs(:NLK) .> NLK_UPPER] .= NLK_UPPER
            xnew.tx_pwrs(:NML)[xnew.tx_pwrs(:NML) .< NML_LOWER] .= NML_LOWER
            xnew.tx_pwrs(:NML)[xnew.tx_pwrs(:NML) .> NML_UPPER] .= NML_UPPER
            state.tx_pwrs(t=i) .= xnew.tx_pwrs
        end

        if :rx in statetypes
            # rx_phi_logpost was mutated in place via state.rx_phi_logpost(t=i).
            state.rx_phi_offset(t=i) .= xnew.rx_phi_offset
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