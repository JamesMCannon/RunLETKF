#Testing
using Pkg
Pkg.activate(@__DIR__)

using JLD2
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, waitsparameter
using AxisKeys, StableRNGs, Parameters
using ScatteredInterpolation, Proj 
using LinearAlgebra, Distributions
using Dates

using Base.Threads

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

"LWPC sometimes fails with low β. Clip the minimum β value to forward models at MIN_BETA."
const MIN_BETA = 0.22

"Sets size of simulated data - effectively sets the maximum number of LETKF iterations."
const DATALENGTH = 10
dt = DateTime(2020, 3, 1, 20, 00) #using Dates

σamp, σphase = 0.1, deg2rad(1.0)

σTX = parse(Float64, get(ENV, "STDDEV_TX_KW", "50")) * 1000 #convert to watts

σNLK, σNML = σTX, σTX

datatypes=(:amp, :phase)

ens_size = parse(Int, get(ENV, "ENS_SIZE", "50"))
ntimes = parse(Int, get(ENV, "ITRS", "1"))

ρ = parse(Float64, get(ENV, "RHO","1.1"))

NLKb = parse(Float64, get(ENV, "NLKB","250")) * 1000
NMLb = parse(Float64, get(ENV, "NLKB","233")) * 1000

rng = reset_rng()

data = observations("Inputs/bigdif_tx_pwrs_params1")

σTXkw = Int(σTX/1000)

scenario = "large_delta_tx_pwrs_$(ntimes)itr_$(ens_size)ens_$(σTXkw)txkW_$(ρ)"

parameters() = merge(init_params(), (;data, σamp, σphase, dt, rng, scenario, σNLK, σNML, datatypes, ens_size, ntimes, ρ, NLKb, NMLb))

#@unpack ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
#datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, localization, scenario = parameters()


isdir(resdir(scenario)) || mkdir(resdir(scenario))

### Run LETKF here
state, data, ym = runletkf(parameters)