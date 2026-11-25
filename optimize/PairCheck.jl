"""
    PairCheck

Checks whether close raw roots of the modal equation are one mode or several, and whether a
production mesh's root list is missing any verified mode. Depends only on LongwaveModePropagator,
so the same file serves any environment whose regions have the fields `re`, `im`, `owner`,
`maxadaptivenodes`, `maxnodes` and `seed` (degrees; finite `im` bounds).

For one segment, [`checksegment`](@ref):
1. searches every region with RootsAndPoles at the production settings but no duplicate removal,
   keeping each raw root's region and the poles;
2. groups raw roots closer than `radius` (single linkage) into clusters;
3. Newton-polishes each cluster member twice, once with its step capped at 0.4 of the distance to
   its nearest neighbour and once uncapped (with a Vern7 retry), and keeps the distinct results;
4. counts the zeros in a circle around the cluster with the argument principle, adding the poles
   inside, and runs Newton on F deflated by the roots already found while fewer are located;
5. checks each located root with a circle of radius min(gap/10, pole distance/2, 1e-4) whose
   winding number must be +1;
6. lists every converged root that the production root list lacks (roots within `guard90` of 90°
   excepted).

[`AugmentedMesh`](@ref) adds given roots to a mesh's own, segment by segment, for field comparisons.
"""
module PairCheck

using LongwaveModePropagator, Random, Logging
const LMP = LongwaveModePropagator

export checksegment, AugmentedMesh

const MERGE = 1e-7    # rad; polished results closer than this are one root (production dedupe)
const MATCH = 1e-6    # rad; a root within this of a production root is present there

function rawroots(me, regions; tol=3e-5, remax=89.9, params=LMPParams())
    roots, boxes, poles = ComplexF64[], Int[], ComplexF64[]
    for (k, r) in pairs(regions)
        isinf(r.im[1]) && error("PairCheck needs finite region bounds")
        re, im = deg2rad.(r.re), deg2rad.(r.im)
        pts = ComplexF64[complex(re[1], im[1]), complex(re[1], im[2]), complex(re[2], im[1]), complex(re[2], im[2])]
        fp = LMP.FinderParams(tol=tol, maxadaptivenodes=r.maxadaptivenodes, maxnodes=r.maxnodes, numtasks=1)
        p = LMPParams(params; finderparams=fp)
        f = θ -> LMP.solvemodalequation(θ, me; params=p)
        run() = (Random.seed!(r.seed);
                 LMP.RootsAndPoles.rootsandpoles(f, LMP.RootsAndPoles.ComplexMesh(pts), fp; dedupe=false))
        saved = copy(Random.default_rng())
        zr, zp = try
            with_logger(run, NullLogger())
        finally
            copy!(Random.default_rng(), saved)
        end
        lo, hi = deg2rad(r.owner[1]), deg2rad(r.owner[2])
        for x in zr
            imag(x) < 0 && real(x) ≤ deg2rad(remax) && lo ≤ real(x) < hi || continue
            push!(roots, x); push!(boxes, k)
        end
        append!(poles, zp)
    end
    return roots, boxes, poles
end

"""
    newton(me, θ0; params, maxit=8, xtol=1e-8, maxstep=5e-4, deflate=()) → (θ, ok)

Newton on F/∏(θ - d) for d in `deflate`; fails when a step is not finite or θ moves more than
`maxstep` from θ0. Without `deflate` the steps are those of production polishing.
"""
function newton(me, θ0; params=LMPParams(), maxit=8, xtol=1e-8, maxstep=5e-4, deflate=())
    θ = θ0
    try
        for _ in 1:maxit
            dF, R, Rg = LMP.solvedmodalequation(θ, me; params)
            F = LMP.modalequation(R, Rg)
            δ = isempty(deflate) ? F/dF : 1/(dF/F - sum(1/(θ - d) for d in deflate))
            isfinite(δ) || return θ0, false
            θ -= δ
            abs(θ - θ0) > maxstep && return θ0, false
            abs(δ) < xtol && return θ, true
        end
    catch err
        err isa InterruptException && rethrow()
    end
    return θ0, false
end

tightparams(params) = LMPParams(params; integrationparams=LMP.IntegrationParams(solver=LMP.Vern7(), tolerance=1e-9))

"""Production polishing: Newton capped at `maxstep`, then a Vern7 retry at 1e-9."""
function polish(me, θ0; params=LMPParams(), maxstep=5e-4)
    θ, ok = newton(me, θ0; params, maxstep)
    ok && return θ, true
    return newton(me, θ0; params=tightparams(params), maxit=20, maxstep)
end

"""
    winding(me, c, r; params, n0=64, nmax=2048) → (w, n, reliable)

Winding number of F on the circle |θ - c| = r, doubling the points until no step of arg F
exceeds π/3.
"""
function winding(me, c, r; params=LMPParams(), n0=64, nmax=2048)
    n = n0
    while true
        z = c .+ r .* cis.(range(0, 2π; length=n+1))
        f = try
            [LMP.solvemodalequation(θ, me; params) for θ in z]
        catch err
            err isa InterruptException && rethrow()
            return NaN, n, false
        end
        d = angle.(f[2:end] ./ f[1:end-1])
        all(isfinite, d) || return NaN, n, false
        maximum(abs, d) < π/3 && return sum(d)/2π, n, true
        n ≥ nmax && return sum(d)/2π, n, false
        n *= 2
    end
end

function addunique!(D, θ)
    any(d -> abs(d - θ) < MERGE, D) || push!(D, θ)
    return D
end

function clusters(roots, radius)
    n = length(roots)
    label = collect(1:n)
    find(i) = (while label[i] != i; i = label[i]; end; i)
    for i in 1:n, j in i+1:n
        abs(roots[i] - roots[j]) < radius && (label[find(j)] = find(i))
    end
    groups = Dict{Int,Vector{Int}}()
    for i in 1:n
        push!(get!(groups, find(i), Int[]), i)
    end
    return sort!(collect(values(groups)); by=first)
end

function checkcluster(me, roots, boxes, poles, idx; params)
    raw = roots[idx]
    gaps = [minimum(abs(raw[i] - raw[j]) for j in eachindex(raw) if j != i) for i in eachindex(raw)]
    capped = [newton(me, raw[i]; params, maxstep=0.4gaps[i]) for i in eachindex(raw)]
    uncapped = [polish(me, x; params) for x in raw]
    D = ComplexF64[]
    for (θ, ok) in Iterators.flatten((capped, uncapped))
        ok && addunique!(D, θ)
    end

    c = sum(raw)/length(raw)
    margin = max(0.5minimum(gaps), 1e-6)
    radius = maximum(abs(x - c) for x in [raw; D]) + margin
    pin = count(p -> abs(p - c) < radius, poles)
    pnear = count(p -> abs(abs(p - c) - radius) < 0.1radius, poles)
    w, npts, reliable = winding(me, c, radius; params)
    Z = isfinite(w) ? round(Int, w) + pin : -1

    for start in [raw; c]
        length(D) ≥ Z && break
        θ, ok = newton(me, start; params, maxit=30, maxstep=2radius, deflate=D)
        ok && abs(θ - c) < radius && addunique!(D, θ)
    end

    verified = map(D) do d
        others = [abs(d - e) for e in D if e != d]
        pd = isempty(poles) ? Inf : minimum(abs(p - d) for p in poles)
        r = min(isempty(others) ? Inf : minimum(others)/10, pd/2, 1e-4)
        wd, _, rel = winding(me, d, r; params)
        rel && isfinite(wd) && round(Int, wd) == 1
    end

    class = !reliable || pnear > 0 || Z < 1 ? "inconclusive" :
            Z == 1 && length(D) == 1 && all(verified) ? "single" :
            Z ≥ 2 && length(D) == Z && all(verified) ? "distinct" : "inconclusive"
    return (; members=idx, raw, boxes=boxes[idx], gaps, capped, uncapped, roots=D, verified,
            center=c, radius, winding=w, npoints=npts, reliable, poles_inside=pin,
            poles_near_circle=pnear, zeros=Z, class)
end

"""
    checksegment(me, regions; params, production, radius=1e-3, tol=3e-5, remax=89.9, guard90=1e-3)

Pair check of one segment (see [`PairCheck`](@ref)). `production` is the production mesh's root
list for the same `me` and `params` (radians). Returns a NamedTuple.
"""
function checksegment(me, regions; params=LMPParams(), production, radius=1e-3, tol=3e-5,
        remax=89.9, guard90=1e-3)
    roots, boxes, poles = rawroots(me, regions; tol, remax, params)
    groups = clusters(roots, radius)
    multi = [checkcluster(me, roots, boxes, poles, g; params) for g in groups if length(g) > 1]

    found = ComplexF64[]
    for g in groups
        length(g) == 1 || continue
        θ, ok = polish(me, roots[only(g)]; params)
        ok && addunique!(found, θ)
    end
    for cl in multi, d in cl.roots
        addunique!(found, d)
    end
    prod = ComplexF64.(production)
    absent = [θ for θ in found if abs(θ - π/2) ≥ guard90 && all(p -> abs(p - θ) ≥ MATCH, prod)]
    unmatched = [p for p in prod if all(θ -> abs(θ - p) ≥ MATCH, found)]
    return (; raw=roots, boxes, poles, clusters=multi, found, production=prod, absent, unmatched)
end

"""
    AugmentedMesh(inner, extra)

`inner`'s roots plus each root in `extra` (radians) that Newton confirms within 1e-6 rad on the
segment being solved and that `inner` lacks.
"""
struct AugmentedMesh{M}
    inner::M
    extra::Vector{ComplexF64}
end

function LMP.findmodes(me::LMP.ModeEquation, mesh::AugmentedMesh; params=LMPParams())
    modes = LMP.findmodes(me, mesh.inner; params)
    θs = ComplexF64[m.θ for m in modes]
    for e in mesh.extra
        θ, ok = newton(me, e; params, maxstep=1e-6)
        ok && all(x -> abs(x - θ) ≥ MATCH, θs) && push!(θs, θ)
    end
    return LMP.EigenAngle.(θs)
end

end
