# Loads what Main.jl loads, here and on one worker, and runs the AxisKeys symbol lookup that fails
# when packages are evaluated without their precompile caches.
using JLD2, LongwaveModePropagator, AxisKeys, StableRNGs, Parameters, Random, ScatteredInterpolation, Proj,
    LinearAlgebra, Distributions, StatsBase, Dates, ProgressMeter, LMPTools, WaitProfileEstimators,
    ModifiedVLFInversionAlgorithms, Distributed
k = KeyedArray(zeros(2, 3); field=[:amp, :phase], x=1:3)
size(k(:amp)) == (3,) || error("AxisKeys lookup failed")
w = only(addprocs(1; exeflags="--project=$(Base.active_project())"))
Distributed.remotecall_eval(Main, [w], :(using ModifiedVLFInversionAlgorithms, AxisKeys, StatsBase))
remotecall_fetch(() -> KeyedArray([1.0, 2.0]; field=[:a, :b])(:b), w) == 2.0 || error("worker AxisKeys lookup failed")
rmprocs(w)
println("load test ok: ", ModifiedVLFInversionAlgorithms.AVIDMESH isa PolishedMesh{TiledMesh})
