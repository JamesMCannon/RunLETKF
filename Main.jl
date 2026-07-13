### Package Management
using Pkg
Pkg.activate(@__DIR__)

using JLD2
using LongwaveModePropagator
using LongwaveModePropagator: QE, ME, C0, waitsparameter
using AxisKeys, StableRNGs, Parameters, Random
using ScatteredInterpolation, Proj 
using LinearAlgebra, Distributions, StatsBase
using Dates

using Base.Threads
using ProgressMeter

using LMPTools, WaitProfileEstimators

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
const MAX_BETA = 0.90
const MIN_H = 55
const MAX_H = 90

# Sets size of simulated data - effectively sets the maximum number of filter iterations.
const DATALENGTH = 10

### Set Common parameters
params = init_params()
do_rx = parse(Bool, get(ENV, "DO_RX", "false"))
do_tx = parse(Bool, get(ENV, "DO_TX", "false"))
if do_rx
    params = init_rx_params(params)
end
if do_tx
    params = init_tx_params(params)
end
if do_rx || do_tx
    params = init_filter_params(params)
end

### Cross-validate state and observable selections before any expensive work.
# TX power estimation updates against :amp observations from per-TX paths;
# categorical RX (Bϕ) estimation accumulates evidence from :phase observations.
if do_tx && !(:amp in params.datatypes)
    error("TX power estimation (DO_TX=true) requires :amp observations (DO_AMP=true)")
end
if do_rx && !(:phase in params.datatypes)
    error("RX offset estimation (DO_RX=true) requires :phase observations (DO_PHASE=true)")
end

# Observations are built for exactly the fields in datatypes (plus their
# _noiseless copies), so no post-hoc field slicing is performed.
data = observations(params.datafile, params.σobs, params.datatypes; paths=params.paths)

if :rx in params.statetypes
    # The per-path Bϕ offset is a property of the (single, synchronized)
    # demodulation trellis: it shifts the absolute Hy phase only. The ratio
    # observables (:s2, :s3) are invariant to any trellis offset common to the
    # two loop channels, so they receive no offset.
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
@info "datatypes: " params.datatypes

do_tx = :tx in params.statetypes
do_rx = :rx in params.statetypes
epp = params.epp

if !do_tx && !do_rx
    if epp == "none"
        scenario = "baseline_"
    else
        scenario = epp*"_"
    end
elseif do_tx && !do_rx
    scenario = "tx_pwrs_"
elseif !do_tx && do_rx
    scenario = "rx_sig_"   
else
    scenario = "tx_rx_"
end

scenario = name_scenario(scenario, parameters)

params = merge(params, (; scenario))

isdir(resdir(scenario)) || mkpath(resdir(scenario))

state, data, ym = runletkf(parameters)