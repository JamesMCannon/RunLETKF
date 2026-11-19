include(joinpath(@__DIR__, "setup.jl"))

harness_config() = (;
    mode    = get(ENV, "H_MODE", "bench"),         # bench | fidelity
    state   = get(ENV, "H_STATE", "g1/state_t1"),
    t       = parse(Int, get(ENV, "H_T", "1")),
    members = parse(Int, get(ENV, "H_MEMBERS", "40")),
    pathset = get(ENV, "H_PATHSET", "bench4"),
    calls   = parse(Int, get(ENV, "H_CALLS", "3")),
    out     = get(ENV, "H_OUT", "g1/run"),
    gclog   = parse(Bool, get(ENV, "H_GCLOG", "false")))

"""
    timed_forward(p, paths, rows, lk, t0) → f

`f(z)` is the LETKF forward model for one member; each call appends
`(thread, start, stop)` (s since `t0[]`) to `rows`.
"""
function timed_forward(p, paths, rows, lk, t0)
    return function (z)
        tid = threadid()
        a = time_ns()
        r = model_observables(p.itp, z.xy_state, paths, p.dt; pathstep=p.pathstep, datatypes=p.datatypes)
        b = time_ns()
        lock(lk) do
            push!(rows, (; thread=tid, start=(a - t0[])/1e9, stop=(b - t0[])/1e9))
        end
        return r
    end
end

function summarize(k, m, rows)
    d = [r.stop - r.start for r in rows]
    busy = sum(d) / (nthreads() * m.wall)
    @printf("call %d: wall %.1f s, GC %.1f s (%.1f%%), %.1f GiB, pauses %d, full %d, maxrss %.1f GiB; member s min/med/max %.1f/%.1f/%.1f, busy %.1f%%\n",
        k, m.wall, m.gctime, 100m.gctime/m.wall, m.bytes/2^30, m.pauses, m.full, m.maxrss/2^30,
        minimum(d), median(d), maximum(d), 100busy)
end

"""
    fidelity_check(ym, ref)

Compare harness `ym` with the LETKF's saved `ym` for the same members and paths: linear
fields must be identical; phase fields may differ by whole turns from ensemble wrapping.
"""
function fidelity_check(ym, ref)
    ok = true
    for fld in axiskeys(ym, :field)
        a = raw(ym(field=fld))
        b = raw(ref(field=fld, path=axiskeys(ym, :path), ens=axiskeys(ym, :ens)))
        if MVIA.is_phase_field(fld)
            δ = maximum(abs, rem2pi.(a .- b, RoundNearest))
            good = δ ≤ 1e-12
            @printf("%-6s max |Δ mod 2π| = %.3g  %s\n", fld, δ, good ? "ok" : "MISMATCH")
        else
            good = isequal(a, b)
            @printf("%-6s %s\n", fld, good ? "identical" : "MISMATCH, max |Δ| = $(maximum(abs, a .- b))")
        end
        ok &= good
    end
    println(ok ? "fidelity: PASS" : "fidelity: FAIL")
    return ok
end

function main()
    c = harness_config()
    outdir = joinpath(ROOT, c.out)
    isdir(outdir) && error("$outdir exists; refusing to overwrite")
    c.gclog && GC.enable_logging(true)
    prov = provenance()
    @info "harness" c prov

    p = init_params()
    saved, file = load_state(c.state)
    xy = saved["state"].xy_state
    (axiskeys(xy, :x) == p.x_grid && axiskeys(xy, :y) == p.y_grid) ||
        error("state grid in $file differs from init_params grid")
    c.members ≤ length(axiskeys(xy, :ens)) || error("state has fewer than $(c.members) members")

    x = (; xy_state = xy(t=c.t)[ens=1:c.members])
    paths = selectpaths(p.paths, c.pathset)
    rows = NamedTuple{(:thread, :start, :stop), Tuple{Int,Float64,Float64}}[]
    lk = ReentrantLock()
    t0 = Ref(time_ns())

    wpaths = selectpaths(p.paths, "warm")
    fw = timed_forward(p, wpaths, rows, lk, t0)
    tw = @elapsed ensemble_model!(newym(p.datatypes, wpaths, 1), fw, (; xy_state = x.xy_state[ens=1:1]))
    @printf("warm-up (1 member, NML-PIN): %.1f s\n", tw)

    f = timed_forward(p, paths, rows, lk, t0)
    results, timelines, yms = NamedTuple[], Vector{eltype(rows)}[], KeyedArray[]
    for k in 1:c.calls
        ym = newym(p.datatypes, paths, c.members)
        empty!(rows)
        t0[] = time_ns()
        m = measure(() -> ensemble_model!(ym, f, x); mark=c.gclog ? "call $k" : "")
        summarize(k, m, rows)
        push!(results, (; k, m.wall, m.gctime, m.bytes, m.pauses, m.full, m.maxrss))
        push!(timelines, copy(rows))
        push!(yms, ym)
    end

    same = [isequal(raw(yms[1]), raw(yms[k])) for k in 2:length(yms)]
    isempty(same) || println("ym identical across calls in this process: ", all(same), " ", same)
    fid = c.mode == "fidelity" ? fidelity_check(yms[1], saved["ym"](t=c.t)) : missing

    mkpath(outdir)
    jldsave(joinpath(outdir, "harness.jld2"); config=c, provenance=prov, statefile=file,
        paths=MVIA.pathname.(paths), warmup=tw, results, timelines, ym=yms[1], same, fidelity=fid)
    @info "saved" joinpath(outdir, "harness.jld2")
end

main()
