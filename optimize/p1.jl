include(joinpath(@__DIR__, "setup.jl"))
include(joinpath(@__DIR__, "polishmesh.jl"))

"""
    p1_forward(p, paths, mesh, memberof) → f

LETKF forward model for one member with `mesh`; each call tags its searches with the member
number and a per-member search counter.
"""
function p1_forward(p, paths, mesh, memberof)
    return z -> task_local_storage(:p1, (member=memberof[hash(raw(z.xy_state))], seq=Ref(0))) do
        model_observables(p.itp, z.xy_state, paths, p.dt; pathstep=p.pathstep, datatypes=p.datatypes, mesh)
    end
end

function main()
    n = parse(Int, get(ENV, "P_MEMBERS", "100"))
    pathset = get(ENV, "H_PATHSET", "bench4")
    outdir = joinpath(ROOT, get(ENV, "P_OUT", "p1/run"))
    isdir(outdir) && error("$outdir exists; refusing to overwrite")
    prov = provenance()
    @info "p1" n pathset prov

    p = init_params()
    saved, file = load_state(get(ENV, "H_STATE", "g1/state_t1"))
    x = (; xy_state = saved["state"].xy_state(t=1)[ens=1:n])
    paths = selectpaths(p.paths, pathset)
    segments = [length(MVIA.pathpts(tx, rx; dist=p.pathstep)[2]) for (tx, rx) in paths]
    memberof = Dict(hash(raw(x.xy_state(ens=e))) => e for e in axiskeys(x.xy_state, :ens))
    length(memberof) == n || error("member states are not distinct")
    mesh = PolishMesh(SAFEMESH)

    wpaths = selectpaths(p.paths, "warm")
    xw = (; xy_state = x.xy_state[ens=1:1])
    tw = @elapsed ensemble_model!(newym(p.datatypes, wpaths, axiskeys(xw.xy_state, :ens)),
        p1_forward(p, wpaths, mesh, memberof), xw)
    @printf("warm-up (1 member, NML-PIN): %.1f s\n", tw)
    empty!(EVENTS); empty!(SEARCHES)

    ym = newym(p.datatypes, paths, axiskeys(x.xy_state, :ens))
    m = measure(() -> ensemble_model!(ym, p1_forward(p, paths, mesh, memberof), x))
    @printf("call: wall %.1f s, GC %.1f s (%.1f%%), %.1f GiB, pauses %d, full %d, maxrss %.1f GiB; searches %d, roots %d\n",
        m.wall, m.gctime, 100m.gctime/m.wall, m.bytes/2^30, m.pauses, m.full, m.maxrss/2^30,
        length(SEARCHES), length(EVENTS))

    mkpath(outdir)
    jldsave(joinpath(outdir, "p1.jld2"); provenance=prov, statefile=file, members=1:n,
        paths=MVIA.pathname.(paths), segments, warmup=tw,
        result=(; m.wall, m.gctime, m.bytes, m.pauses, m.full, m.maxrss),
        ym, events=copy(EVENTS), searches=copy(SEARCHES))
    @info "saved" joinpath(outdir, "p1.jld2")
end

main()
