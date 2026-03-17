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
# Often functions are called by MVIA.foo() for transparency, even when those same functions are 
# exported explicitely by MVIA (as an example). This is only to aid in debugging and identifying 
# which package is at the root of errors should they arise.

include("common.jl")
include("letkf.jl")

### Set Global Parameters
# Constrain values to physically realistic (but wide) bounds
const MIN_BETA = 0.18 # Gasdia used 0.22 when running with LWPC
const MAX_BETA = 0.55
const MIN_H = 55
const MAX_H = 90

# Sets size of simulated data - effectively sets the maximum number of LETKF iterations.
const DATALENGTH = 10

### Set Common parameters
params = init_params()
if parse(Bool, get(ENV, "DO_RX", "false"))
    params = init_rx_params(params)
end
if parse(Bool, get(ENV, "DO_TX", "false"))
    params = init_tx_params(params)
end
if parse(Bool, get(ENV, "DO_DUAL", "false"))
    params = init_dual_params(params)
end

data = observations(params.datafile, params.σamp, params.σphase, paths=params.paths)

if !(:amp in params.datatypes && :phase in params.datatypes)
    if :phase in params.datatypes
        data = data[Key([:phase, :phase_noiseless]), :, :]
    elseif :amp in params.datatypes
        data = data[Key([:amp, :amp_noiseless]), :, :]
    else
        error("Whoops! Datatypes don't make sense")
    end
end

if :rx in params.statetypes
    data(:phase) .+= params.rx_offsets .* (π/2)
    data(:phase_noiseless) .+= params.rx_offsets .* (π/2)
end

params = merge(params, (; data))

parameters() = params

# Path at which to write output.
if params.new_folder == "false"
    const RESDIR = Ref(abspath(joinpath(@__DIR__, "results")))
else
    const RESDIR = Ref(abspath(joinpath(@__DIR__, params.new_folder*"/results")))
end
resdir() = RESDIR[]
resdir(d) = joinpath(resdir(), d)
resdir!(d) = isdir(d) ? RESDIR[] = d : throw(ArgumentError("$(abspath(d)) must be a directory"))
resdir!() = RESDIR[] = normpath(joinpath(@__DIR__))

if :tx in params.statetypes
    # TODO revise this whole logic using constants. 
    # Perhaps just set TX_range as the constant and caclulate the bounds in the LETKF code?
    const NML_LOWER = max(1, params.NMLb - params.TX_range/2)
    const NML_UPPER = params.NMLb + params.TX_range/2
    const NLK_LOWER = max(1, params.NLKb - params.TX_range/2)
    const NLK_UPPER = params.NLKb + params.TX_range/2
end

### Bring it all together
@info "statetypes: " params.statetypes

do_tx = :tx in params.statetypes
do_rx = :rx in params.statetypes

if !do_tx && !do_rx
    scenario = "baseline_"
elseif do_tx && !do_rx
    scenario = "tx_pwrs_"
elseif !do_tx && do_rx
    scenario = "rx_offset_"   
else
    scenario = "tx_rx_"
end

name_scenario!(scenario, parameters)

params = merge(params, (; scenario))

isdir(resdir(scenario)) || mkdir(resdir(scenario))

state, data, ym = runletkf(parameters)



