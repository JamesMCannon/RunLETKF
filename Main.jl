#Testing
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

NLKb = parse(Float64, get(ENV, "NLKB","250")) * 1000
NMLb = parse(Float64, get(ENV, "NLKB","233")) * 1000

σTX = parse(Float64, get(ENV, "STDDEV_TX_KW", "50")) * 1000 #convert to watts

σNLK, σNML = σTX, σTX

do_amp = parse(Bool, get(ENV, "DO_AMP", "true"))
do_phase = parse(Bool, get(ENV, "DO_PHASE", "true"))

datatypes = ()

if do_amp 
    datatypes = (datatypes..., :amp)
end

if do_phase
    datatypes = (datatypes..., :phase)
end

datatypes=(:amp, :phase)

do_tx = parse(Bool, get(ENV, "DO_TX", "true"))
do_xy = parse(Bool, get(ENV, "DO_XY", "true")) #These sets are for the tx power testing.
do_rx = parse(Bool, get(ENV, "DO_RX", "false"))

statetypes = ()

if do_xy
    statetypes = (statetypes..., :xy)
end

if do_tx
    statetypes = (statetypes..., :tx)
end

if do_rx
    statetypes = (statetypes..., :rx)
end

@info "statetypes: " statetypes

TX_range = parse(Float64, get(ENV, "TX_RANGE_KW", "500")) * 1000 #convert to watts

const NML_LOWER = max(1, NMLb - TX_range/2)
const NML_UPPER = NMLb + TX_range/2
const NLK_LOWER = max(1, NLKb - TX_range/2)
const NLK_UPPER = NLKb + TX_range/2

ens_size = parse(Int, get(ENV, "ENS_SIZE", "50"))
ntimes = parse(Int, get(ENV, "ITRS", "1"))

shuffle_tx = parse(Bool, get(ENV, "SHUFFLE_TX", "false"))
shuffle_xy = parse(Bool, get(ENV, "SHUFFLE_XY", "false"))

ρ = parse(Float64, get(ENV, "RHO","1.1"))

rng = reset_rng()

data = observations("Inputs/day1",σamp, σphase)

σTXkw = Int(σTX/1000)

if shuffle_tx && shuffle_xy
    scenario = "tx_pwrs_$(ntimes)itr_$(ens_size)ens_$(σTXkw)txkW_$(ρ)_shuffle_all_limited"
elseif shuffle_tx
    scenario = "tx_pwrs_$(ntimes)itr_$(ens_size)ens_$(σTXkw)txkW_$(ρ)_shuffle_$(shuffle_tx)"
else
    scenario = "baseline_$(ntimes)itr_$(ens_size)_shuffle_$(shuffle_xy)_limited"
end

if !(do_amp && do_phase)
    if do_amp
        scenario = scenario* "_amp_only"
    elseif do_phase
        scenario = scenario * "_phase_only"
    end 
end

if do_tx
    if TX_range != 500*1000.0
        scenario = scenario * "_tx_constrained_$(TX_range/1000)kW"
    end
end

@info "Starting " scenario

parameters() = merge(init_params(), (;data, σamp, σphase, dt, rng, scenario, σNLK, σNML, NLKb, NMLb,
datatypes, ens_size, ntimes, ρ, statetypes, shuffle_tx, shuffle_xy))

#@unpack ens_size, ntimes, dt, pathstep, x_grid, y_grid, modelsteps,
#datatypes, h0, b0, hB, bB, rng, σamp, σphase, data, itp, localization, scenario = parameters()


isdir(resdir(scenario)) || mkdir(resdir(scenario))

### Run LETKF here
state, data, ym = runletkf(parameters)