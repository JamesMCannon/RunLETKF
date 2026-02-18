### Package Management
using Pkg
Pkg.activate(@__DIR__)

using JLD2
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, waitsparameter
using AxisKeys, StableRNGs, Parameters, Random
using ScatteredInterpolation, Proj 
using LinearAlgebra, Distributions, StatsBase
using Dates

using Base.Threads
using ProgressMeter

using LMPTools

using ModifiedVLFInversionAlgorithms

const MVIA = ModifiedVLFInversionAlgorithms
const SI = ScatteredInterpolation

include("common.jl")
include("letkf.jl")

### Set Global Parameters
# Path at which to write output.
const RESDIR = Ref(abspath(joinpath(@__DIR__, "results")))
resdir() = RESDIR[]
resdir(d) = joinpath(resdir(), d)
resdir!(d) = isdir(d) ? RESDIR[] = d : throw(ArgumentError("$(abspath(d)) must be a directory"))
resdir!() = RESDIR[] = normpath(joinpath(@__DIR__))

# Constrain values to physically realistic (but wide) bounds
const MIN_BETA = 0.22
const MAX_BETA = 0.55
const MIN_H = 55
const MAX_H = 90

# Sets size of simulated data - effectively sets the maximum number of LETKF iterations.
const DATALENGTH = 10

### Set Common parameters

σamp, σphase = 0.1, deg2rad(1.0)

do_amp = parse(Bool, get(ENV, "DO_AMP", "true"))
do_phase = parse(Bool, get(ENV, "DO_PHASE", "true"))
xy_file = get(ENV, "XY_FILE","false")

datatypes = ()

if do_amp 
    datatypes = (datatypes..., :amp)
end

if do_phase
    datatypes = (datatypes..., :phase)
end
do_xy = parse(Bool, get(ENV, "DO_XY", "true")) # Currently, only "true" is supported

ens_size = parse(Int, get(ENV, "ENS_SIZE", "4"))
ntimes = parse(Int, get(ENV, "ITRS", "1"))
shuffle_xy = parse(Bool, get(ENV, "SHUFFLE_XY", "false"))

ρ = parse(Float64, get(ENV, "RHO","1.1"))

rng = reset_rng()

data = observations("Inputs/day1", σamp, σphase)

if !(:amp in datatypes && :phase in datatypes)
    if :phase in datatypes
        data = data[Key([:phase, :phase_noiseless]), :, :]
    elseif :amp in datatypes
        data = data[Key([:amp, :amp_noiseless]), :, :]
    else
        error("Whoops! Datatypes don't make sense")
    end
end

statetypes = ()

if do_xy
    statetypes = (statetypes..., :xy)
end

params = init_params()
parameters() = params

### Set TX parameters
do_tx = parse(Bool, get(ENV, "DO_TX", "false"))
precondition_tx = parse(Bool, get(ENV, "PRECONDITION_TX", "false"))
if do_tx
    
    NLKb = parse(Float64, get(ENV, "NLKB","250")) * 1000
    NMLb = parse(Float64, get(ENV, "NLKB","233")) * 1000

    σTX = parse(Float64, get(ENV, "STDDEV_TX_KW", "50")) * 1000 #convert to watts
    σNLK, σNML = σTX, σTX

    TX_range = parse(Float64, get(ENV, "TX_RANGE_KW", "500")) * 1000 #convert to watts

    const NML_LOWER = max(1, NMLb - TX_range/2)
    const NML_UPPER = NMLb + TX_range/2
    const NLK_LOWER = max(1, NLKb - TX_range/2)
    const NLK_UPPER = NLKb + TX_range/2

    shuffle_tx = parse(Bool, get(ENV, "SHUFFLE_TX", "false"))

    σTXkw = Int(σTX/1000)

    if precondition_tx #theres got to be a better way to do this...
        precon_ens_size = parse(Int, get(ENV, "PRECON_ENS_SIZE", "-1")) # Default value means use ens_size
        if precon_ens_size == -1
            precon_ens_size = ens_size
        end
        precon_itrs = parse(Int, get(ENV, "PRECON_ITRS", "-1"))
        if precon_itrs == -1
            precon_itrs = ntimes
        end
        do_dual = parse(Bool, get(ENV, "DO_DUAL", "false"))
        params = merge(params, (; precon_ens_size, precon_itrs))
    end

    statetypes = (statetypes..., :tx)

    params = merge(params, (; σNLK, σNML, shuffle_tx, NLKb, NMLb))
end

### Set RX Parameters
do_rx = parse(Bool, get(ENV, "DO_RX", "false"))
precondition_rx = parse(Bool, get(ENV, "PRECONDITION_RX", "false"))
if do_rx
    shuffle_rx = parse(Bool, get(ENV, "SHUFFLE_RX", "false"))
    statetypes = (statetypes..., :rx)
    if precondition_rx && !precondition_tx #theres got to be a better way to do this...
        precon_ens_size = parse(Int, get(ENV, "PRECON_ENS_SIZE", "-1")) # Default value means use ens_size
        if precon_ens_size == -1
            precon_ens_size = ens_size
        end
        precon_itrs = parse(Int, get(ENV, "PRECON_ITRS", "-1"))
        if precon_itrs == -1
            precon_itrs = ntimes
        end
        do_dual = parse(Bool, get(ENV, "DO_DUAL", "false"))
        params = merge(params, (; precon_ens_size, precon_itrs))
    end
    params = merge(params, (; precondition_rx, shuffle_rx))
end

### Bring it all together

@info "statetypes: " statetypes

if !do_tx && !do_rx
    scenario = "baseline_"
elseif do_tx && !do_rx
    scenario = "tx_pwrs_"
elseif !do_tx && do_rx
    scenario = "rx_offset_"   
else
    scenario = "tx_rx_"
end

scenario = scenario * "$(ntimes)itr_$(ens_size)ens_$(ρ)_"

check = false
if do_tx
    if shuffle_tx
        check=true
    end
end
if do_rx
    if shuffle_rx
        check=true
    end
end
if shuffle_xy
    check=true
end
if check
    scenario = scenario * "shuffle_"
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
    else
        scenario = scenario * "_tx"
    end
    scenario = scenario * "_log10"
end

if do_rx
    scenario = scenario * "_rx"
end

if precondition_rx || precondition_tx
    if do_dual
        scenario = scenario*"_dual_$(precon_itrs)dualitrs_$(precon_ens_size)dualens"
    else
        scenario = scenario * "_preconditioned_v2"
    end
end

if xy_file != "false"
    scenario = scenario * "_xy_file"
    @info "Background Ionosphere from file"
end

params = merge(params, (; data, σamp, σphase, dt, rng, scenario, datatypes, ens_size, ntimes, ρ, statetypes, shuffle_xy, xy_file))

isdir(resdir(scenario)) || mkdir(resdir(scenario))

### Run LETKF here


if do_dual
    @info "Running Dual LETKF"
    state, data, ym = rundualletkf(parameters)
else
    @info "Running LETKF"
    state, data, ym = runletkf(parameters)
end


