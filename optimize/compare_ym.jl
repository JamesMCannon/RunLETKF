using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using JLD2, AxisKeys

const ROOT = normpath(joinpath(@__DIR__, ".."))
raw(a) = collect(AxisKeys.keyless_unname(a))
loadym(dir) = load(joinpath(ROOT, dir, "harness.jld2"), "ym")

ref = loadym(ARGS[1])
ok = true
for d in ARGS[2:end]
    y = loadym(d)
    println("== ", ARGS[1], " vs ", d)
    if axiskeys(ref) != axiskeys(y)
        println("DIFFERENT axis keys"); global ok = false; continue
    end
    for fld in axiskeys(ref, :field)
        a, b = raw(ref(field=fld)), raw(y(field=fld))
        same = isequal(a, b)
        println(rpad(string(fld), 8), same ? "identical" : "DIFFERENT, max |Δ| = $(maximum(abs, filter(isfinite, a .- b)))")
        global ok &= same
    end
end
println(ok ? "ym: PASS" : "ym: FAIL")
