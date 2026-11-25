# Pair check of one LETKF state (optimize/PairCheck.jl): this array task's (member, path) units,
# one result file each in pc/$PC_LABEL/units/. Finished units are skipped, so a requeued task resumes.
# Env: PC_STATE (state dir, as H_STATE), PC_LABEL, PC_T (default 1), PC_NTASKS, PC_RADIUS (rad,
# default 1e-3), SLURM_ARRAY_TASK_ID.
include(joinpath(@__DIR__, "setup.jl"))
include(joinpath(@__DIR__, "PairCheck.jl"))

const FIELDS = (:amp, :phase, :s2, :s3)

batchof(input) = (b = BatchInput{ExponentialInput}(); b.name = "paircheck"; b.description = "";
    b.datetime = now(); b.inputs = [input]; b)

observables(input, mesh, path) =
    vec(raw(MVIA._collect_observables(MVIA.avidrun(batchof(input); mesh), [path], FIELDS)))

"""
    rununit(p, xm, path; radius) → NamedTuple

Builds the path's input from member state `xm` as `model_observables` does, checks every
segment against `MVIA.AVIDMESH`, and, when a verified root is absent from it, compares the
receiver observables with and without the absent roots.
"""
function rununit(p, xm, path; radius)
    tx, rx = path
    v = [vec(xm(:h)); vec(xm(:b))]
    n = length(v) ÷ 2
    hitp = SI.interpolate(p.itp.method, p.itp.coords, filter(!isnan, v[1:n]))
    bitp = SI.interpolate(p.itp.method, p.itp.coords, filter(!isnan, v[n+1:end]))
    input = MVIA.model_observation(p.itp, hitp, bitp, tx, rx, p.dt; pathstep=p.pathstep,
        fieldcomponent=LMP.Fields.H)
    H = MVIA.waittopheight(input.hprimes, input.betas)
    params = MVIA.pathparams(H)
    freq = LMP.Frequency(input.frequency)
    segments = map(eachindex(input.segment_ranges)) do k
        me = LMP.PhysicalModeEquation(freq, LMP.buildwaveguide(input, k))
        production = [m.θ for m in LMP.findmodes(me, MVIA.AVIDMESH; params)]
        PairCheck.checksegment(me, MVIA.avidregions(); params, production, radius)
    end
    extra = unique(reduce(vcat, [s.absent for s in segments]; init=ComplexF64[]))
    fields = missing
    if !isempty(extra)
        base = observables(input, MVIA.AVIDMESH, path)
        aug = observables(input, PairCheck.AugmentedMesh(MVIA.AVIDMESH, extra), path)
        Δ = aug .- base
        Δ[2] = rad2deg(rem2pi(Δ[2], RoundNearest))
        fields = (; base, augmented=aug, delta=NamedTuple{FIELDS}(Tuple(Δ)))
    end
    return (; topheight=H, hprimes=input.hprimes, betas=input.betas, segments, extra, fields)
end

function main()
    state, label = ENV["PC_STATE"], ENV["PC_LABEL"]
    T = parse(Int, get(ENV, "PC_T", "1"))
    ntasks = parse(Int, get(ENV, "PC_NTASKS", "1"))
    task = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))
    radius = parse(Float64, get(ENV, "PC_RADIUS", "1e-3"))
    outdir = joinpath(ROOT, "pc", label, "units")
    mkpath(outdir)
    prov = (; provenance()..., cpu=Sys.cpu_info()[1].model, task, ntasks)
    @info "paircheck" state label T radius prov

    p = init_params()
    saved, file = load_state(state)
    xy = saved["state"].xy_state
    (axiskeys(xy, :x) == p.x_grid && axiskeys(xy, :y) == p.y_grid) ||
        error("state grid in $file differs from init_params grid")
    members = collect(axiskeys(xy, :ens))
    units = [(m, i) for m in members for i in eachindex(p.paths)]
    mine = [u for (k, u) in enumerate(units) if (k - 1) % ntasks == task - 1]
    unitfile((m, i)) = joinpath(outdir, @sprintf("m%03d_%s.jld2", m, MVIA.pathname(p.paths[i])))
    todo = filter(u -> !isfile(unitfile(u)), mine)
    @printf("task %d/%d: %d units, %d already done, %d to run on %d threads\n",
        task, ntasks, length(mine), length(mine) - length(todo), length(todo), nthreads())

    queue = Channel{Tuple{Int,Int}}(length(todo))
    foreach(u -> put!(queue, u), todo)
    close(queue)
    failed = Threads.Atomic{Int}(0)
    @sync for _ in 1:nthreads()
        Threads.@spawn for (m, i) in queue
            path = p.paths[i]
            try
                t = @elapsed r = rununit(p, xy(t=T)(ens=m), path; radius)
                f = unitfile((m, i))
                tmp = f * ".tmp" * string(rand(UInt32); base=16)
                jldsave(tmp; result=r, member=m, path=MVIA.pathname(path), state=file, t=T, radius,
                    seconds=t, host=gethostname(), cpu=prov.cpu, commit=prov.commit, finished=now())
                mv(tmp, f; force=true)
                @printf("member %d %s: %.0f s, %d segments, %d clusters, %d absent\n", m,
                    MVIA.pathname(path), t, length(r.segments), sum(s -> length(s.clusters), r.segments),
                    length(r.extra))
            catch err
                err isa InterruptException && rethrow()
                Threads.atomic_add!(failed, 1)
                @error "unit failed" member=m path=MVIA.pathname(path) exception=(err, catch_backtrace())
            end
        end
    end
    failed[] > 0 && error("$(failed[]) units failed; rerun the task to retry them")
end

main()
