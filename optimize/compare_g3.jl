using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PHASE_FIELDS = (:phase,)  # as MVIA.PHASE_FIELDS
raw(a) = collect(AxisKeys.keyless_unname(a))

"""
    procdirs(dir)

The harness output directories of one configuration: `dir` itself for a single process, or its
subdirectories (`p0`, `p1`, …) for concurrent processes.
"""
procdirs(dir) = isfile(joinpath(ROOT, dir, "harness.jld2")) ? [joinpath(ROOT, dir)] :
    sort(filter(isdir, readdir(joinpath(ROOT, dir); join=true)))

function loadconfig(dir)
    hs = [load(joinpath(d, "harness.jld2")) for d in procdirs(dir)]
    r = [h["results"][1] for h in hs]
    p = hs[1]["provenance"]
    ts = [x.t_start for x in r]
    return (; dir, procs=length(hs), threads=p.nthreads, gcthreads=p.ngcthreads,
        heap=p.heap_size_hint, numtasks=p.numtasks, members=sum(length(h["members"]) for h in hs),
        span=maximum(x.t_end for x in r) - minimum(ts), startspread=maximum(ts) - minimum(ts),
        gc=sum(x.gctime for x in r), rss=sum(x.maxrss for x in r), yms=[h["ym"] for h in hs])
end

"""
    exactness(cfg, ref) → String

Amp, s2 and s3 must equal `ref` for the same members and paths. Phase must equal `ref` when one
process held all of `ref`'s members; otherwise it may differ by whole turns from wrapping.
"""
function exactness(cfg, ref)
    whole = cfg.procs == 1 && axiskeys(cfg.yms[1], :ens) == axiskeys(ref, :ens)
    δphase = 0.0
    for y in cfg.yms, fld in axiskeys(y, :field)
        a = raw(y(field=fld))
        b = raw(ref(field=fld, path=axiskeys(y, :path), ens=axiskeys(y, :ens)))
        if fld in PHASE_FIELDS && !whole
            δphase = max(δphase, maximum(abs, rem2pi.(a .- b, RoundNearest)))
        elseif !isequal(a, b)
            return "MISMATCH in $fld, max |Δ| = $(maximum(abs, filter(isfinite, a .- b)))"
        end
    end
    whole && return "identical"
    return δphase ≤ 1e-12 ? @sprintf("identical (phase mod 2π, %.1g)", δphase) :
        @sprintf("MISMATCH in phase mod 2π, %.3g", δphase)
end

function main(refdir, basedir, dirs...)
    ref = load(joinpath(ROOT, refdir, "harness.jld2"), "ym")
    base = loadconfig(basedir)
    rate(c) = c.span / c.members
    @printf("%-22s %8s %5s %5s %6s %6s %8s %9s %7s %8s %9s %7s  %s\n", "config", "procs×thr", "gc", "nt",
        "heap", "memb", "span s", "s/member", "speedup", "GC s", "RSS GiB", "skew s", "ym vs $refdir")
    for d in (basedir, dirs...)
        c = loadconfig(d)
        @printf("%-22s %8s %5d %5d %6s %6d %8.1f %9.2f %7.3f %8.1f %9.1f %7.1f  %s\n", d,
            "$(c.procs)×$(c.threads)", c.gcthreads, c.numtasks, c.heap == 0 ? "none" : "$(round(Int, c.heap/2^30))G",
            c.members, c.span, rate(c), rate(base)/rate(c), c.gc, c.rss/2^30, c.startspread, exactness(c, ref))
    end
end

main(ARGS...)
