# Usage: compare_g5.jl CANDIDATE:BASELINE ...  (run directories relative to RunLETKF)
# Each pair must give exactly equal state, ym (filled time steps) and data; exits 1 otherwise.
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys, Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
raw(a) = collect(AxisKeys.keyless_unname(a))
same(x::KeyedArray, y::KeyedArray) = axiskeys(x) == axiskeys(y) && isequal(raw(x), raw(y))
same(x, y) = isequal(x, y)

function loadrun(dir)
    res = joinpath(ROOT, dir, "results")
    sub = only(filter(d -> isdir(joinpath(res, d)), readdir(res)))
    d = load(joinpath(res, sub, sub * ".jld2"))
    g = load(joinpath(ROOT, dir, "g5run.jld2"))
    pm = joinpath(ROOT, dir, "peakmem.txt")
    mem = isfile(pm) ? parse.(Int, split(read(pm, String))) : Int[]
    return (; dir, state=d["state"], ym=d["ym"], data=d["data"], run=g["result"], config=g["config"],
        prov=g["provenance"], mem)
end

layout(r) = r.config.DISTRIBUTED == "true" ? "$(r.config.NPROCS)×$(r.config.PROC_THREADS)" : "1×$(r.prov.nthreads)"

"""Names of the parts of `a` and `b` that differ; the last `ym` time step is unfilled by design."""
function differences(a, b)
    out = String[]
    keys(a.state) == keys(b.state) || return ["state fields $(keys(a.state)) vs $(keys(b.state))"]
    for k in keys(a.state)
        same(a.state[k], b.state[k]) && continue
        x, y = a.state[k], b.state[k]
        δ = x isa KeyedArray && size(x) == size(y) ? maximum(abs, filter(isfinite, raw(x) .- raw(y)); init=0.0) : NaN
        push!(out, "state.$k (max |Δ| $δ)")
    end
    axiskeys(a.ym) == axiskeys(b.ym) && isequal(raw(a.ym)[:, :, :, 1:end-1], raw(b.ym)[:, :, :, 1:end-1]) ||
        push!(out, "ym")
    same(a.data, b.data) || push!(out, "data")
    return out
end

function main(args)
    pairs = [Tuple(split(a, ':')) for a in args]
    all(p -> length(p) == 2, pairs) || error("arguments must be CANDIDATE:BASELINE")
    runs = Dict{String,Any}()
    for p in pairs, d in p
        haskey(runs, d) || (runs[d] = loadrun(d))
    end
    @printf("%-22s %-6s %4s %-8s %4s %-30s %9s %22s %7s %7s %8s %8s\n", "run", "layout", "ens", "paths", "itrs",
        "host", "total s", "iterations s", "GC main", "GC wkr", "RSS GiB", "node GiB")
    for d in sort(collect(keys(runs)))
        r = runs[d]
        w = r.run.workers
        gw = isempty(w) ? "—" : @sprintf("%.1f%%", 100sum(x -> x.gctime, w) / length(w) / r.run.wall)
        rss = (r.run.maxrss + sum(x -> x.maxrss, w; init=0)) / 2^30
        node = length(r.mem) == 2 ? @sprintf("%.1f", (r.mem[2] - r.mem[1]) / 1024) : "—"
        @printf("%-22s %-6s %4s %-8s %4s %-30s %9.1f %22s %6.1f%% %7s %8.1f %8s\n", d, layout(r), r.config.ENS_SIZE,
            r.config.PATHS, r.config.ITRS, r.prov.host, r.run.wall, join(round.(Int, r.run.iterations), ", "),
            100r.run.gctime / r.run.wall, gw, rss, node)
    end
    println()
    ok = true
    for (c, b) in pairs
        rc, rb = runs[c], runs[b]
        diff = differences(rc, rb)
        ok &= isempty(diff)
        @printf("%s vs %s: %s; speedup per iteration %.3f, total %.3f\n", c, b,
            isempty(diff) ? "identical" : "MISMATCH in " * join(diff, ", "),
            sum(rb.run.iterations) / sum(rc.run.iterations), rb.run.wall / rc.run.wall)
    end
    ok || exit(1)
end

main(ARGS)
