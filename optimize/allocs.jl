include(joinpath(@__DIR__, "setup.jl"))
using Profile

const OWNERS = ("LongwaveModePropagator", "RootsAndPoles", "ModifiedVLFInversionAlgorithms")

pkgof(file) = (m = match(r"/packages/([^/]+)/", file); m === nothing ? nothing : String(m.captures[1]))

frame(fr, pkg) = "$(fr.func) $(pkg)/$(basename(string(fr.file))):$(fr.line)"

"""
    attribute(st) → (site, owner, package)

`site`: innermost frame in any package; `owner`: innermost frame in LMP, RootsAndPoles or MVIA;
`package`: the package of `site`.
"""
function attribute(st)
    site, owner, pkg = "(no package frame)", "(no LMP/RaP/MVIA frame)", "(none)"
    for fr in st
        p = pkgof(string(fr.file))
        p === nothing && continue
        if site == "(no package frame)"
            site, pkg = frame(fr, p), p
        end
        if p in OWNERS
            owner = frame(fr, p)
            break
        end
    end
    return site, owner, pkg
end

function tally(allocs, key, rate)
    d = Dict{String,Tuple{Float64,Float64}}()
    for a in allocs
        k = key(a)
        b, n = get(d, k, (0.0, 0.0))
        d[k] = (b + a.size/rate, n + 1/rate)
    end
    return sort!(collect(d); by=x -> -x[2][1])
end

function printtop(title, rows, total; n=20)
    println("\n== ", title)
    @printf("%10s %7s %12s  %s\n", "MiB", "share", "allocs", "location")
    for (k, (b, c)) in first(rows, n)
        @printf("%10.1f %6.1f%% %12.3g  %s\n", b/2^20, 100b/total, c, k)
    end
end

function main()
    nthreads() == 1 || @warn "allocs.jl is meant for one thread; got $(nthreads())"
    state = get(ENV, "H_STATE", "g1/state_t1")
    t = parse(Int, get(ENV, "H_T", "1"))
    member = parse(Int, get(ENV, "S_MEMBER", "1"))
    pathset = get(ENV, "A_PATHSET", "warm")
    rate = parse(Float64, get(ENV, "A_RATE", "0.001"))
    outdir = joinpath(ROOT, get(ENV, "A_OUT", "g2/M4"))
    isdir(outdir) && error("$outdir exists; refusing to overwrite")
    prov = provenance()
    @info "allocs" state t member pathset rate prov

    p = init_params()
    saved, file = load_state(state)
    z = saved["state"].xy_state(t=t, ens=member)
    paths = selectpaths(p.paths, pathset)
    fwd() = model_observables(p.itp, z, paths, p.dt; pathstep=p.pathstep, datatypes=p.datatypes)

    fwd()
    GC.gc()
    m = @timed fwd()
    @printf("unprofiled call: %.1f s, %.1f MiB allocated\n", m.time, m.bytes/2^20)

    GC.gc()
    Profile.Allocs.clear()
    tp = @elapsed Profile.Allocs.@profile sample_rate=rate fwd()
    allocs = Profile.Allocs.fetch().allocs
    est = sum(a.size for a in allocs; init=0)/rate
    @printf("profiled call: %.1f s, %d samples at rate %g, estimated %.1f MiB (%.0f%% of unprofiled)\n",
        tp, length(allocs), rate, est/2^20, 100est/m.bytes)

    bysite = tally(allocs, a -> attribute(a.stacktrace)[1], rate)
    byowner = tally(allocs, a -> attribute(a.stacktrace)[2], rate)
    bypkg = tally(allocs, a -> attribute(a.stacktrace)[3], rate)
    bytype = tally(allocs, a -> string(a.type), rate)
    printtop("by innermost package frame", bysite, est)
    printtop("by innermost LMP / RootsAndPoles / MVIA frame", byowner, est)
    printtop("by package of innermost frame", bypkg, est; n=15)
    printtop("by type", bytype, est; n=15)

    mkpath(outdir)
    jldsave(joinpath(outdir, "allocs.jld2"); provenance=prov, statefile=file, member, t,
        paths=MVIA.pathname.(paths), rate, unprofiled=(; time=m.time, bytes=m.bytes),
        bysite, byowner, bypkg, bytype)
    @info "saved" joinpath(outdir, "allocs.jld2")
end

main()
