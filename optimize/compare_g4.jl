using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys, Printf, Statistics, Distributions

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BASE = "L1x40"
raw(a) = collect(AxisKeys.keyless_unname(a))

"""`ym` saved in `dir` by the harness (`harness.jld2`) or by P1 (`p1.jld2`)."""
function loadym(dir)
    f = joinpath(ROOT, dir, "harness.jld2")
    return load(isfile(f) ? f : joinpath(ROOT, dir, "p1.jld2"), "ym")
end

"""Every field of `ym` must equal `ref` for the same members and paths, phase included."""
function exact(ym, ref)
    for fld in axiskeys(ym, :field)
        a = raw(ym(field=fld))
        b = raw(ref(field=fld, path=axiskeys(ym, :path), ens=axiskeys(ym, :ens)))
        isequal(a, b) || return "MISMATCH in $fld, max |Δ| = $(maximum(abs, filter(isfinite, a .- b)))"
    end
    return "identical"
end

layout(h) = h["config"].procs == 0 ? "1×$(h["provenance"].nthreads)" :
    "$(h["config"].procs)×$(h["config"].proc_threads)"

"""
    compare_g4(refdir, repdirs)

Per run: layout, wall, GC (per-worker mean), RSS and `ym` against `refdir`. Across replicates:
speedup of each configuration over `L1x40` on the same node, from log wall times, with a t-based
95% CI on the replicate-paired differences and the two-way residual SD.
"""
function compare_g4(refdir, repdirs)
    ref = loadym(refdir)
    cfgs = sort(filter(d -> isdir(joinpath(ROOT, repdirs[1], d)), readdir(joinpath(ROOT, repdirs[1]))))
    BASE in cfgs || error("$BASE missing in $(repdirs[1])")
    cfgs = [BASE; filter(!=(BASE), cfgs)]
    R, C = length(repdirs), length(cfgs)
    y = fill(NaN, R, C)
    @printf("%-14s %-6s %-26s %8s %8s %7s %8s  %s\n", "run", "layout", "host", "wall s", "GC s", "GC %",
        "RSS GiB", "ym vs $refdir")
    for (i, d) in enumerate(repdirs), (j, c) in enumerate(cfgs)
        f = joinpath(ROOT, d, c, "harness.jld2")
        isfile(f) || (@warn "missing" f; continue)
        h = load(f)
        r = h["results"][1]
        y[i, j] = log(r.wall)
        @printf("%-14s %-6s %-26s %8.1f %8.1f %7.1f %8.1f  %s\n", "$(basename(d))/$c", layout(h),
            h["provenance"].host, r.wall, r.gctime, 100r.gctime/r.wall, r.maxrss/2^30, exact(h["ym"], ref))
    end
    R < 2 && return
    any(isnan, y) && (println("\nincomplete design; no replicate analysis"); return)
    resid = y .- mean(y; dims=2) .- mean(y; dims=1) .+ mean(y)
    σ = sqrt(sum(abs2, resid) / ((R - 1) * (C - 1)))
    tq = quantile(TDist(R - 1), 0.975)
    @printf("\nreplicates %d; run-to-run SD %.2f%%\n", R, 100(exp(σ) - 1))
    @printf("%-8s %9s %9s %18s\n", "config", "wall s", "speedup", "95% CI")
    for (j, c) in enumerate(cfgs)
        d = y[:, 1] .- y[:, j]
        m, se = mean(d), std(d) / sqrt(R)
        @printf("%-8s %9.1f %9.3f %8.3f – %-7.3f\n", c, exp(mean(y[:, j])), exp(m),
            j == 1 ? 1.0 : exp(m - tq*se), j == 1 ? 1.0 : exp(m + tq*se))
    end
end

compare_g4(ARGS[1], ARGS[2:end])
