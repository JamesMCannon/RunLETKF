include(joinpath(@__DIR__, "setup.jl"))
using Distributed

harness_config() = (;
    mode    = get(ENV, "H_MODE", "bench"),         # bench | fidelity
    state   = get(ENV, "H_STATE", "g1/state_t1"),
    t       = parse(Int, get(ENV, "H_T", "1")),
    first   = parse(Int, get(ENV, "H_FIRST", "1")),
    members = parse(Int, get(ENV, "H_MEMBERS", "40")),
    pathset = get(ENV, "H_PATHSET", "bench4"),
    calls   = parse(Int, get(ENV, "H_CALLS", "3")),
    out     = get(ENV, "H_OUT", "g1/run"),
    gclog   = parse(Bool, get(ENV, "H_GCLOG", "false")),
    numtasks = parse(Int, get(ENV, "H_NUMTASKS", "0")),     # 0: lmpparams() default
    start_at = parse(Float64, get(ENV, "H_START_AT", "0")),  # epoch s; 0: no barrier
    procs    = parse(Int, get(ENV, "H_PROCS", "0")),         # 0: this process's threads only
    proc_threads = parse(Int, get(ENV, "H_PROC_THREADS", "10")),
    proc_flags = get(ENV, "H_PROC_FLAGS", ""))

"""
    timed_forward(p, paths, rows, lk, t0, params) → f

`f(z)` is the LETKF forward model for one member; each call appends
`(thread, start, stop)` (s since `t0[]`) to `rows`.
"""
function timed_forward(p, paths, rows, lk, t0, params)
    return function (z)
        tid = threadid()
        a = time_ns()
        r = model_observables(p.itp, z.xy_state, paths, p.dt; pathstep=p.pathstep, datatypes=p.datatypes,
            params)
        b = time_ns()
        lock(lk) do
            push!(rows, (; thread=tid, start=(a - t0[])/1e9, stop=(b - t0[])/1e9))
        end
        return r
    end
end

"""
    start_workers(c) → pids

Add `c.procs` worker processes with `c.proc_threads` threads and `c.proc_flags`, and load the
packages the forward model needs on them.
"""
function start_workers(c)
    pids = addprocs(c.procs; exeflags=`--project=$ROOT --threads=$(c.proc_threads) $(split(c.proc_flags))`)
    Distributed.remotecall_eval(Main, pids, :(using ModifiedVLFInversionAlgorithms, AxisKeys, StatsBase))
    return pids
end

"""
    measure_procs(f, pids)

As `measure`, with GC and allocation totals taken from the workers `pids`. `gctime` is the
mean GC time per worker, comparable with one process's GC time; `maxrss` sums all processes.
"""
function measure_procs(f, pids)
    GC.gc()
    @sync for w in pids
        @async remotecall_wait(GC.gc, w)
    end
    g0 = [remotecall_fetch(Base.gc_num, w) for w in pids]
    t = @timed f()
    ds = [Base.GC_Diff(remotecall_fetch(Base.gc_num, w), g) for (w, g) in zip(pids, g0)]
    gcs = [d.total_time/1e9 for d in ds]
    return (; wall=t.time, gctime=mean(gcs), gcmax=maximum(gcs), bytes=sum(d -> d.allocd, ds),
        pauses=sum(d -> d.pause, ds), full=sum(d -> d.full_sweep, ds),
        maxrss=Sys.maxrss() + sum(w -> remotecall_fetch(Sys.maxrss, w), pids), value=t.value)
end

function summarize(k, m, rows)
    if isempty(rows)
        @printf("call %d: wall %.1f s, GC per worker mean %.1f s (%.1f%%) max %.1f s, %.1f GiB, pauses %d, full %d, maxrss %.1f GiB\n",
            k, m.wall, m.gctime, 100m.gctime/m.wall, m.gcmax, m.bytes/2^30, m.pauses, m.full, m.maxrss/2^30)
        return
    end
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
    c.procs > 0 && (c.numtasks > 0 || c.gclog) && error("H_NUMTASKS and H_GCLOG apply only without H_PROCS")
    c.gclog && GC.enable_logging(true)
    prm = finderparams(c.numtasks)
    pids = c.procs > 0 ? start_workers(c) : nothing
    wk = pids === nothing ? missing : [(; pid=w, nthreads=remotecall_fetch(Threads.nthreads, w),
        ngcthreads=remotecall_fetch(Threads.ngcthreads, w)) for w in pids]
    prov = (; provenance()..., numtasks=pids === nothing ? prm.finderparams.numtasks : c.proc_threads,
        workers=wk)
    @info "harness" c prov

    p = init_params()
    saved, file = load_state(c.state)
    xy = saved["state"].xy_state
    (axiskeys(xy, :x) == p.x_grid && axiskeys(xy, :y) == p.y_grid) ||
        error("state grid in $file differs from init_params grid")
    members = c.first:(c.first + c.members - 1)
    last(members) ≤ length(axiskeys(xy, :ens)) || error("state has fewer than $(last(members)) members")

    x = (; xy_state = xy(t=c.t)[ens=members])
    paths = selectpaths(p.paths, c.pathset)
    rows = NamedTuple{(:thread, :start, :stop), Tuple{Int,Float64,Float64}}[]
    lk = ReentrantLock()
    t0 = Ref(time_ns())

    # Workers run the forward model with lmpparams() evaluated on the worker (numtasks = its threads).
    function forward(pths)
        pids === nothing && return timed_forward(p, pths, rows, lk, t0, prm)
        itp, dt, pathstep, datatypes = p.itp, p.dt, p.pathstep, p.datatypes
        return z -> model_observables(itp, z.xy_state, pths, dt; pathstep, datatypes)
    end

    wpaths = selectpaths(p.paths, "warm")
    nw = pids === nothing ? 1 : length(pids)
    xw = (; xy_state = x.xy_state[ens=1:nw])
    tw = @elapsed ensemble_model!(newym(p.datatypes, wpaths, axiskeys(xw.xy_state, :ens)), forward(wpaths), xw; pids)
    @printf("warm-up (1 member per process, NML-PIN): %.1f s\n", tw)

    if c.start_at > 0
        wait_s = c.start_at - time()
        if wait_s > 0
            sleep(wait_s)
        else
            @warn "start barrier missed" late_s = -wait_s
        end
    end

    f = forward(paths)
    results, timelines, yms = NamedTuple[], Vector{eltype(rows)}[], KeyedArray[]
    for k in 1:c.calls
        ym = newym(p.datatypes, paths, members)
        empty!(rows)
        t0[] = time_ns()
        ta = time()
        m = pids === nothing ? measure(() -> ensemble_model!(ym, f, x); mark=c.gclog ? "call $k" : "") :
            measure_procs(() -> ensemble_model!(ym, f, x; pids), pids)
        tb = time()
        summarize(k, m, rows)
        push!(results, (; k, m.wall, m.gctime, gcmax=get(m, :gcmax, m.gctime), m.bytes, m.pauses, m.full,
            m.maxrss, t_start=ta, t_end=tb))
        push!(timelines, copy(rows))
        push!(yms, ym)
    end

    same = [isequal(raw(yms[1]), raw(yms[k])) for k in 2:length(yms)]
    isempty(same) || println("ym identical across calls in this process: ", all(same), " ", same)
    fid = c.mode == "fidelity" ? fidelity_check(yms[1], saved["ym"](t=c.t)) : missing

    mkpath(outdir)
    jldsave(joinpath(outdir, "harness.jld2"); config=c, provenance=prov, statefile=file,
        paths=MVIA.pathname.(paths), members, warmup=tw, results, timelines, ym=yms[1], same, fidelity=fid)
    @info "saved" joinpath(outdir, "harness.jld2")
end

main()
