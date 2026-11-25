# H2 inputs for Mesh_Testing's lmp_truth_states.jl: one record per (member, path) of a saved LETKF
# state, with the production F1r roots that the pair check stored for each segment.
# Env: H2_STATE (state dir, as PC_STATE), H2_LABEL (pair-check label; output in h2/$H2_LABEL/),
# H2_T (default 1), H2_MEMBERS (optional, e.g. "1:10" or "3,17,42"), H2_NAME (file name, default inputs).
include(joinpath(@__DIR__, "setup.jl"))
using Serialization

const INPUT_FIELDS = (:segment_ranges, :hprimes, :betas, :b_mags, :b_dips, :b_azs, :ground_sigmas,
    :ground_epsrs, :frequency)

const HTOL, BTOL = 1e-6, 1e-8   # km, km⁻¹: far below physical significance, far above round-off

parsemembers(s) = occursin(':', s) ? (r = parse.(Int, split(s, ':')); collect(r[1]:r[2])) :
    parse.(Int, split(s, ','))

"""Plain-typed copy of `input`'s fields, so the record deserializes under another LMP commit."""
plain(input) = NamedTuple{INPUT_FIELDS}(Tuple(copy(getfield(input, f)) for f in INPUT_FIELDS))

function record(p, xm, path, m, pcdir)
    tx, rx = path
    h, b = xm(:h), xm(:b)
    isnan.(h) == isnan.(b) || error("h′ and β have different active cells")
    coords = MVIA.build_xygrid(h)
    hitp = SI.interpolate(p.itp.method, coords, filter(!isnan, vec(h)))
    bitp = SI.interpolate(p.itp.method, coords, filter(!isnan, vec(b)))
    input = MVIA.model_observation(p.itp, hitp, bitp, tx, rx, p.dt; pathstep=p.pathstep,
        fieldcomponent=LMP.Fields.H)
    H = MVIA.waittopheight(input.hprimes, input.betas)
    name = MVIA.pathname(path)
    pc = load(joinpath(pcdir, @sprintf("m%03d_%s.jld2", m, name)), "result")
    n = length(input.segment_ranges)
    (length(pc.hprimes) == n && length(pc.segments) == n) ||
        error("member $m $name: $n segments here, $(length(pc.hprimes)) in the pair check")
    dh, db = maximum(abs, pc.hprimes .- input.hprimes), maximum(abs, pc.betas .- input.betas)
    (dh ≤ HTOL && db ≤ BTOL) ||
        error(@sprintf("member %d %s: pair check's h′ differs by %.3g km, β by %.3g km⁻¹", m, name, dh, db))
    # h′, β and top height exactly as the pair check found the roots for
    input.hprimes .= pc.hprimes
    input.betas .= pc.betas
    roots = [ComplexF64.(s.production) for s in pc.segments]
    return (; member=m, path=name, rxdist=Float64(range(tx, rx)), input=plain(input), roots,
        H=pc.topheight, dh, db, dH=abs(pc.topheight - H))
end

function main()
    state, label = ENV["H2_STATE"], ENV["H2_LABEL"]
    T = parse(Int, get(ENV, "H2_T", "1"))
    pcdir = joinpath(ROOT, "pc", label, "units")
    isdir(pcdir) || error("no pair-check units in $pcdir")
    p = init_params()
    saved, file = load_state(state)
    xy = saved["state"].xy_state
    (axiskeys(xy, :x) == p.x_grid && axiskeys(xy, :y) == p.y_grid) ||
        error("state grid in $file differs from init_params grid")
    allm = collect(axiskeys(xy, :ens))
    members = haskey(ENV, "H2_MEMBERS") ? parsemembers(ENV["H2_MEMBERS"]) : allm
    issubset(members, allm) || error("members not in the state: $(setdiff(members, allm))")

    records = [record(p, xy(t=T)(ens=m), path, m, pcdir) for m in members for path in p.paths]
    outdir = joinpath(ROOT, "h2", label)
    mkpath(outdir)
    out = joinpath(outdir, get(ENV, "H2_NAME", "inputs") * ".jls")
    serialize(out, records)

    nseg = sum(r -> length(r.roots), records)
    Hs = sort([r.H for r in records])
    seglen = [diff([r.input.segment_ranges; r.rxdist]) for r in records]
    worst = maximum(s -> maximum(abs.(s .- s[1])), seglen)
    @printf("%s: %d records (%d members × %d paths), %d segments, from %s\n", out, length(records),
        length(members), length(p.paths), nseg, file)
    @printf("top height %.1f–%.1f km, median %.1f; h′ %.1f–%.1f km; β %.3f–%.3f\n", Hs[1]/1e3, Hs[end]/1e3,
        Hs[(end + 1) ÷ 2]/1e3, extrema(reduce(vcat, [r.input.hprimes for r in records]))...,
        extrema(reduce(vcat, [r.input.betas for r in records]))...)
    @printf("largest segment-length spread within a path: %.3g m\n", worst)
    @printf("largest rebuild difference from the pair check: h′ %.3g km, β %.3g km⁻¹, top height %.3g m\n",
        maximum(r -> r.dh, records), maximum(r -> r.db, records), maximum(r -> r.dH, records))
    @printf("provenance: %s\n", provenance())
end

main()
