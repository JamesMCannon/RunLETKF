using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys, Printf, Statistics, Distributions

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIGS = ["base", "gc40", "gc20s1", "gc40s1", "h150", "gc40_h150", "gc20s1_h150", "gc40s1_h150"]
raw(a) = collect(AxisKeys.keyless_unname(a))

"""
    r1r2_analysis(refdir, repdirs)

Two-way (replicate × configuration) analysis of log wall time. Each replicate ran on one node, so
replicate effects absorb node differences. The residual SD is the run-to-run SD of a single run;
a configuration is kept if its mean improvement is at least max(5%, 2 × that SD).
"""
function r1r2_analysis(refdir, repdirs)
    ref = raw(load(joinpath(ROOT, refdir, "harness.jld2"), "ym"))
    R, C = length(repdirs), length(CONFIGS)
    y = fill(NaN, R, C); gc = similar(y); pauses = similar(y); full = similar(y); rss = similar(y)
    hosts = String[]
    for (i, d) in enumerate(repdirs)
        for (j, c) in enumerate(CONFIGS)
            f = joinpath(ROOT, d, c, "harness.jld2")
            isfile(f) || (@warn "missing" f; continue)
            h = load(f)
            r = h["results"][1]
            y[i, j], gc[i, j], pauses[i, j], full[i, j], rss[i, j] = log(r.wall), r.gctime, r.pauses, r.full, r.maxrss
            j == 1 && push!(hosts, h["provenance"].host)
            isequal(raw(h["ym"]), ref) || @warn "ym differs from $refdir" d c
        end
    end
    any(isnan, y) && error("incomplete design; cannot fit the two-way model")

    resid = y .- mean(y; dims=2) .- mean(y; dims=1) .+ mean(y)
    σ = sqrt(sum(abs2, resid) / ((R - 1) * (C - 1)))
    threshold = max(log(1.05), 2σ)
    tq = quantile(TDist(R - 1), 0.975)
    @printf("replicates %d (%s); run-to-run SD %.2f%%; keep threshold %.2f%%\n\n",
        R, join(unique(hosts), ", "), 100(exp(σ) - 1), 100(exp(threshold) - 1))
    @printf("%-13s %9s %9s %18s %9s %9s %8s %9s  %s\n", "config", "wall s", "speedup", "95% CI", "GC s",
        "pauses", "full", "RSS GiB", "verdict")
    for (j, c) in enumerate(CONFIGS)
        d = y[:, 1] .- y[:, j]
        m, se = mean(d), std(d) / sqrt(R)
        verdict = j == 1 ? "" : m ≥ threshold ? "keep" : "within noise"
        @printf("%-13s %9.1f %9.3f %8.3f – %-7.3f %9.1f %9.0f %8.0f %9.1f  %s\n", c, exp(mean(y[:, j])),
            exp(m), exp(m - tq*se), exp(m + tq*se), mean(gc[:, j]), mean(pauses[:, j]), mean(full[:, j]),
            mean(rss[:, j])/2^30, verdict)
    end
end

r1r2_analysis(ARGS[1], ARGS[2:end])
