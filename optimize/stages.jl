include(joinpath(@__DIR__, "setup.jl"))

"""
    TimedMesh(inner)

Mesh wrapper whose `findmodes` calls `findmodes(me, inner)` and adds its time, bytes and GC
time to `SEARCH`. Single-threaded use only.
"""
struct TimedMesh{M}
    inner::M
end

const SEARCH = (n=Ref(0), time=Ref(0.0), bytes=Ref(0), gctime=Ref(0.0))

function LMP.findmodes(me::LMP.ModeEquation, mesh::TimedMesh; params=LMPParams())
    t = @timed LMP.findmodes(me, mesh.inner; params)
    SEARCH.n[] += 1
    SEARCH.time[] += t.time
    SEARCH.bytes[] += t.bytes
    SEARCH.gctime[] += t.gctime
    return t.value
end

stagestats(t) = (; time=t.time, bytes=t.bytes, gctime=t.gctime)

"""
    staged(p, z, paths) → (stats, observables)

`model_observables(p.itp, z, paths, p.dt; ...)` in its four stages, each timed.
"""
function staged(p, z, paths; mesh=SAFEMESH, params=lmpparams())
    xv = [vec(z(:h)); vec(z(:b))]
    npts = length(xv) ÷ 2
    hprimes, betas = xv[1:npts], xv[npts+1:end]

    s1 = @timed (SI.interpolate(p.itp.method, p.itp.coords, filter(!isnan, hprimes)),
                 SI.interpolate(p.itp.method, p.itp.coords, filter(!isnan, betas)))
    hitp, bitp = s1.value

    s2 = @timed begin
        batch = BatchInput{ExponentialInput}()
        batch.name = "estimate"
        batch.description = ""
        batch.datetime = Dates.now()
        batch.inputs = Vector{ExponentialInput}(undef, length(paths))
        for i in eachindex(paths)
            tx, rx = paths[i]
            batch.inputs[i] = MVIA.model_observation(p.itp, hitp, bitp, tx, rx, p.dt;
                pathstep=p.pathstep, fieldcomponent=Fields.H)
        end
        batch
    end

    SEARCH.n[], SEARCH.time[], SEARCH.bytes[], SEARCH.gctime[] = 0, 0.0, 0, 0.0
    s3 = @timed LMP.buildrun(s2.value; mesh=TimedMesh(mesh), params)
    search = (; n=SEARCH.n[], time=SEARCH.time[], bytes=SEARCH.bytes[], gctime=SEARCH.gctime[])

    s4 = @timed MVIA._collect_observables(s3.value, paths, p.datatypes; Hunits_pT=true)

    nseg = sum(length(input.segment_ranges) for input in s2.value.inputs)
    stats = (; nseg, interpolants=stagestats(s1), inputs=stagestats(s2), buildrun=stagestats(s3),
        search, propagation=(; time=s3.time - search.time, bytes=s3.bytes - search.bytes,
            gctime=s3.gctime - search.gctime), collect=stagestats(s4))
    return stats, s4.value
end

function report(stats)
    total = stats.interpolants.time + stats.inputs.time + stats.buildrun.time + stats.collect.time
    tb = stats.interpolants.bytes + stats.inputs.bytes + stats.buildrun.bytes + stats.collect.bytes
    @printf("%-28s %10s %7s %12s %10s %12s\n", "stage", "time (s)", "share", "alloc (MiB)", "GC (s)", "MiB/segment")
    for (name, s) in (("interpolants", stats.interpolants), ("path inputs", stats.inputs),
                      ("buildrun", stats.buildrun), ("  root search", stats.search),
                      ("  everything else", stats.propagation), ("collect", stats.collect))
        @printf("%-28s %10.2f %6.1f%% %12.1f %10.2f %12.2f\n", name, s.time, 100s.time/total,
            s.bytes/2^20, s.gctime, s.bytes/2^20/stats.nseg)
    end
    @printf("%-28s %10.2f %7s %12.1f\n", "total", total, "", tb/2^20)
    @printf("segments %d, searches %d, %.2f s and %.1f MiB per segment\n",
        stats.nseg, stats.search.n, total/stats.nseg, tb/2^20/stats.nseg)
end

function main()
    nthreads() == 1 || @warn "stages.jl is meant for one thread; got $(nthreads())"
    state = get(ENV, "H_STATE", "g1/state_t1")
    t = parse(Int, get(ENV, "H_T", "1"))
    member = parse(Int, get(ENV, "S_MEMBER", "1"))
    pathset = get(ENV, "H_PATHSET", "bench4")
    outdir = joinpath(ROOT, get(ENV, "S_OUT", "g2/M3"))
    isdir(outdir) && error("$outdir exists; refusing to overwrite")
    prov = provenance()
    @info "stages" state t member pathset prov

    p = init_params()
    saved, file = load_state(state)
    z = saved["state"].xy_state(t=t, ens=member)
    paths = selectpaths(p.paths, pathset)

    tw = @elapsed staged(p, z, selectpaths(p.paths, "warm"))
    @printf("warm-up (NML-PIN): %.1f s\n", tw)

    GC.gc()
    stats, obs = staged(p, z, paths)
    report(stats)

    GC.gc()
    ref = model_observables(p.itp, z, paths, p.dt; pathstep=p.pathstep, datatypes=p.datatypes)
    same = isequal(raw(ref), raw(obs))
    println("staged result identical to model_observables: ", same)

    mkpath(outdir)
    jldsave(joinpath(outdir, "stages.jld2"); provenance=prov, statefile=file, member, t,
        paths=MVIA.pathname.(paths), stats, same)
    @info "saved" joinpath(outdir, "stages.jld2")
end

main()
