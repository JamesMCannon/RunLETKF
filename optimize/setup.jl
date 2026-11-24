using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JLD2, Dates, Printf, Statistics
using LongwaveModePropagator
using AxisKeys, StableRNGs, Parameters, Random
using ScatteredInterpolation, Proj
using LinearAlgebra, Distributions, StatsBase
using Base.Threads, ProgressMeter
using LMPTools, WaitProfileEstimators
using ModifiedVLFInversionAlgorithms

const MVIA = ModifiedVLFInversionAlgorithms
const LMP = LongwaveModePropagator
const SI = ScatteredInterpolation
const ROOT = normpath(joinpath(@__DIR__, ".."))
const DATALENGTH = 10  # as in Main.jl; construct_R sizes R with it

include(joinpath(ROOT, "common.jl"))

const PATHSETS = Dict(
    "warm"   => ("NML-PIN",),
    "bench4" => ("NML-PIN", "NLK-FSM", "NLK-GIL", "NML-WHI"),
    "bench2" => ("NLK-WHI", "NML-FSM"))

raw(a) = collect(AxisKeys.keyless_unname(a))

selectpaths(paths, set) = set == "all" ? paths : filter(p -> MVIA.pathname(p) in PATHSETS[set], paths)

newym(datatypes, paths, ens) = KeyedArray(fill(NaN, length(datatypes), length(paths), length(ens));
    field=collect(datatypes), path=MVIA.pathname.(paths), ens=collect(ens))

function load_state(dir)
    res = joinpath(ROOT, dir, "results")
    sub = only(readdir(res; join=true))
    file = only(filter(f -> endswith(f, ".jld2") && !startswith(basename(f), "rx_offsets"),
        readdir(sub; join=true)))
    return load(file), file
end

"""
    measure(f; mark="")

Run `f()` after a full collection and return wall, GC and allocation totals. A non-empty `mark`
writes `== GCLOG start/end <mark>` to stderr around the timed region, for `gcpauses.jl`.
"""
function measure(f; mark="")
    GC.gc()
    isempty(mark) || (println(stderr, "\n== GCLOG start ", mark); flush(stderr))
    g0 = Base.gc_num()
    t = @timed f()
    d = Base.GC_Diff(Base.gc_num(), g0)
    isempty(mark) || (println(stderr, "\n== GCLOG end ", mark); flush(stderr))
    return (; wall=t.time, gctime=t.gctime, bytes=t.bytes, pauses=d.pause, full=d.full_sweep,
        maxrss=Sys.maxrss(), value=t.value)
end

"""
    finderparams(numtasks) → LMPParams

`lmpparams()` with `numtasks` replaced; `numtasks ≤ 0` returns `lmpparams()` unchanged.
"""
function finderparams(numtasks)
    p = lmpparams()
    numtasks ≤ 0 && return p
    f = p.finderparams
    return LMPParams(p; finderparams=FinderParams(; tol=f.tol, maxiters=f.maxiters,
        maxadaptivenodes=f.maxadaptivenodes, maxnodes=f.maxnodes, numtasks))
end

function provenance()
    jo = Base.JLOptions()
    return (; commit=readchomp(`git -C $ROOT rev-parse HEAD`),
        dirty=!isempty(readchomp(`git -C $ROOT status --porcelain`)),
        host=gethostname(), nthreads=nthreads(), ngcthreads=Threads.ngcthreads(),
        heap_size_hint=hasfield(typeof(jo), :heap_size_hint) ? Int(jo.heap_size_hint) : -1,
        julia=string(VERSION), cpu_target=get(ENV, "JULIA_CPU_TARGET", ""), started=now())
end
