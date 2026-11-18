using Pkg
Pkg.activate(@__DIR__)
using JLD2, AxisKeys

raw(a) = collect(AxisKeys.keyless_unname(a))

function load_run(dir)
    res = joinpath(@__DIR__, dir, "results")
    sub = only(readdir(res; join=true))
    file = only(filter(f -> endswith(f, ".jld2") && !startswith(basename(f), "rx_offsets"),
        readdir(sub; join=true)))
    return load(file)
end

function report(a, b, label)
    A, B = raw(a), raw(b)
    axiskeys(a) == axiskeys(b) || (println(rpad(label, 18), "DIFFERENT axis keys"); return false)
    isequal(A, B) && (println(rpad(label, 18), "identical"); return true)
    d = filter(isfinite, A .- B)
    nanmis = count(isnan.(A) .!= isnan.(B))
    println(rpad(label, 18), "DIFFERENT: max |Δ| = ", isempty(d) ? NaN : maximum(abs, d),
        ", NaN-pattern mismatches = ", nanmis)
    return false
end

runs = ARGS
res = [load_run(r) for r in runs]
println("NaN in ym of ", runs[1], ": ", count(isnan, raw(res[1]["ym"])), " of ", length(res[1]["ym"]))

ok = true
for k in 2:length(runs)
    println("== ", runs[1], " vs ", runs[k])
    a, b = res[1], res[k]
    global ok &= report(a["data"], b["data"], "data")
    global ok &= report(a["ym"], b["ym"], "ym")
    for f in keys(a["state"])
        global ok &= report(a["state"][f], b["state"][f], "state.$f")
    end
end
println(ok ? "G0 reproducibility: PASS" : "G0 reproducibility: FAIL")