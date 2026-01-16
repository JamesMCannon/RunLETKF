using Pkg
Pkg.activate(@__DIR__)

using JLD2
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, waitsparameter
using AxisKeys, StableRNGs, Parameters, Random
using ScatteredInterpolation, Proj 
using LinearAlgebra, Distributions
using Dates

using Base.Threads
using ProgressMeter

using LMPTools

using ModifiedVLFInversionAlgorithms

const MVIA = ModifiedVLFInversionAlgorithms
const SI = ScatteredInterpolation

include("common.jl")
include("letkf.jl")

# Path at which to write output.
const RESDIR = Ref(abspath(joinpath(@__DIR__, "results")))
resdir() = RESDIR[]
resdir(d) = joinpath(resdir(), d)
resdir!(d) = isdir(d) ? RESDIR[] = d : throw(ArgumentError("$(abspath(d)) must be a directory"))
resdir!() = RESDIR[] = normpath(joinpath(@__DIR__))

"Constrain values to physically realistic (but wide) bounds"
const MIN_BETA = 0.22
const MAX_BETA = 0.425
const MIN_H = 65
const MAX_H = 80

"Sets size of simulated data - effectively sets the maximum number of LETKF iterations."
const DATALENGTH = 10
dt = DateTime(2020, 3, 1, 20, 00) #using Dates

σamp, σphase = 0.1, deg2rad(1.0)

σTX = parse(Float64, get(ENV, "STDDEV_TX_KW", "50")) * 1000 #convert to watts

shuffle_tx = parse(Bool, get(ENV, "SHUFFLE_TX", "false"))
shuffle_xy = parse(Bool, get(ENV, "SHUFFLE_XY", "false"))
shuffle_rx = parse(Bool, get(ENV, "SHUFFLE_RX", "false"))

precondition_rx = parse(Bool, get(ENV, "PRECONDITION_RX", "false"))

σNLK, σNML = σTX, σTX

datatypes=(:amp, :phase)

statetypes=(:xy, :rx)

ens_size = parse(Int, get(ENV, "ENS_SIZE", "50"))
ntimes = parse(Int, get(ENV, "ITRS", "1"))

ρ = parse(Float64, get(ENV, "RHO","1.1"))

rng = reset_rng()

data = observations("Inputs/day1",σamp, σphase)

σTXkw = Int(σTX/1000)

if shuffle_rx || shuffle_xy
    scenario = "rx_offset_$(ntimes)itr_$(ens_size)ens_$(ρ)_daytime_constrAndResamp_shuffled"
else
    scenario = "rx_offset_$(ntimes)itr_$(ens_size)ens_$(ρ)_daytime_constrAndResamp"
end

if precondition_rx
    scenario = scenario * "_preconditioned_v2"
end

parameters() = merge(init_params(), (;data, σamp, σphase, dt, rng, scenario, σNLK, σNML, datatypes, ens_size, ntimes, ρ, statetypes, shuffle_rx, shuffle_tx, shuffle_xy, precondition_rx))

isdir(resdir(scenario)) || mkdir(resdir(scenario))

### Run LETKF here
state, data, ym = runletkf(parameters)