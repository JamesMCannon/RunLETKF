using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys, Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
raw(a) = collect(AxisKeys.keyless_unname(a))
const MESHLIMIT = Dict(:amp => 0.15, :phase => 1.5, :s2 => 0.019, :s3 => 0.019)    # dB, deg, -
const SIGMA = Dict(:amp => 0.435, :phase => 4.35, :s2 => 0.055, :s3 => 0.055)
const DROPS = (:filter_drop, :guard_drop, :dup_drop)

q(v, p) = isempty(v) ? NaN : quantile(v, p)

function fielddiff(ym, ref, fld)
    a, b = raw(ym(field=fld)), raw(ref(field=fld, path=axiskeys(ym, :path), ens=axiskeys(ym, :ens)))
    return fld == :phase ? rad2deg.(abs.(rem2pi.(a .- b, RoundNearest))) : abs.(a .- b)
end

"""
    p1_report(p1dir, refdir, costdir)

Effect of polishing on `ym` (polished `p1dir` against unpolished `refdir`, same members, paths
and seeded searches), root outcomes, aborted searches, raw roots beyond 89.9°, and cost
(`costdir`: one member on one thread).
"""
function p1_report(p1dir, refdir, costdir)
    P = load(joinpath(ROOT, p1dir, "p1.jld2"))
    ref = load(joinpath(ROOT, refdir, "harness.jld2"), "ym")
    ym, ev, se, segs = P["ym"], P["events"], P["searches"], P["segments"]
    paths = P["paths"]
    bounds = cumsum(segs)
    pathof(seq) = searchsortedfirst(bounds, seq)
    perm = Dict(m => count(s -> s.member == m, se) for m in unique(s.member for s in se))
    counts = collect(values(perm))
    if !all(==(last(bounds)), counts)
        @warn "searches per member differ from the segment count" expected = last(bounds) observed = extrema(counts)
    end

    dropped = Set((e.member, pathof(e.seq)) for e in ev if e.outcome in DROPS)
    println("== effect on ym, polished vs unpolished ($(length(axiskeys(ym, :ens))) members × $(length(paths)) paths)")
    @printf("%-6s %-10s %6s %10s %10s %10s %10s %10s\n", "field", "subset", "n", "median", "p99", "max", ">mesh lim", ">σ")
    for fld in axiskeys(ym, :field)
        d = fielddiff(ym, ref, fld)
        mask = [((m, k) in dropped) for k in eachindex(paths), m in axiskeys(ym, :ens)]
        for (name, sel) in (("all", trues(size(d))), ("no drop", .!mask), ("with drop", mask))
            v = d[sel]
            @printf("%-6s %-10s %6d %10.3g %10.3g %10.3g %10d %10d\n", fld, name, length(v), q(v, 0.5), q(v, 0.99),
                isempty(v) ? NaN : maximum(v), count(>(MESHLIMIT[fld]), v), count(>(SIGMA[fld]), v))
        end
    end
    println("(units: amp dB, phase degrees modulo 360, s2/s3 dimensionless)")

    println("\n== root outcomes ($(length(ev)) roots in $(length(se)) searches)")
    for o in (:converged, :retry_converged, :unconverged_kept, :filter_drop, :guard_drop, :dup_drop)
        @printf("%-17s %8d\n", o, count(e -> e.outcome == o, ev))
    end
    moved = [abs(e.polished - e.raw) for e in ev if e.outcome in (:converged, :retry_converged)]
    @printf("|θ polished − θ raw| (rad): median %.2g, p99 %.2g, max %.2g; > 3e-5: %d\n",
        q(moved, 0.5), q(moved, 0.99), isempty(moved) ? NaN : maximum(moved), count(>(3e-5), moved))

    println("\n== aborted searches")
    for a in (:maxiters, :maxnodes)
        s = [x for x in se if x.abort == a]
        isempty(s) && (println("$a: none"); continue)
        es = [e for e in ev if e.abort == a]
        @printf("%s: %d searches, %d roots: ", a, length(s), length(es))
        println(join(["$o $(count(e -> e.outcome == o, es))" for o in unique(e.outcome for e in es)], ", "))
        for x in s
            @printf("    member %d, %s segment %d\n", x.member, paths[pathof(x.seq)], x.seq - (pathof(x.seq) == 1 ? 0 : bounds[pathof(x.seq)-1]))
        end
    end

    println("\n== events other than plain convergence, and all raw roots with Re θ > 89.9°")
    @printf("%7s %-8s %4s %12s %12s %12s %12s %-17s %8s %10s\n", "member", "path", "seg", "Re raw °", "Im raw °",
        "Re pol °", "Im pol °", "outcome", "winding", "dB/Mm")
    for e in ev
        (e.outcome in (:converged, :retry_converged) && real(e.raw) ≤ deg2rad(89.9)) && continue
        k = pathof(e.seq)
        @printf("%7d %-8s %4d %12.6f %12.6f %12.6f %12.6f %-17s %8.2f %10.3f\n", e.member, paths[k],
            e.seq - (k == 1 ? 0 : bounds[k-1]), rad2deg(real(e.raw)), rad2deg(imag(e.raw)),
            rad2deg(real(e.polished)), rad2deg(imag(e.polished)), e.outcome, e.winding, e.loss)
    end

    r = P["result"]
    @printf("\n== cost\n100-member call: wall %.1f s, GC %.1f s, %.1f GiB, maxrss %.1f GiB\n",
        r.wall, r.gctime, r.bytes/2^30, r.maxrss/2^30)
    C = load(joinpath(ROOT, costdir, "p1.jld2"))
    cs = C["searches"]
    ts, tp, bp = sum(x.t_search for x in cs), sum(x.t_polish for x in cs), sum(x.bytes_polish for x in cs)
    @printf("one member, one thread: wall %.1f s, %.1f GiB; search %.1f s, polishing %.1f s (%.1f%% of search), polishing %.1f MiB (%.2f MiB/segment)\n",
        C["result"].wall, C["result"].bytes/2^30, ts, tp, 100tp/ts, bp/2^20, bp/2^20/length(cs))
end

p1_report(ARGS...)
