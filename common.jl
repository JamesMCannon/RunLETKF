reset_rng() = StableRNG(1234)
reset_rng(seed) = StableRNG(seed)

"""
    buildpaths()

Return a vector of `(Transmitter, Receiver)` propagation paths used in the scenarios.
"""
function buildpaths()
    transmitters = [TRANSMITTER[:NLK], TRANSMITTER[:NML]]
    receivers = [
        Receiver("Whitehorse", 60.724, -135.043, 0.0, VerticalDipole()),
        Receiver("Churchill", 58.74, -94.085, 0.0, VerticalDipole()),
        Receiver("Stony Rapids", 59.253, -105.834, 0.0, VerticalDipole()),
        Receiver("Fort Smith", 60.006, -111.92, 0.0, VerticalDipole()),
        Receiver("Bella Bella", 52.1675508, -128.1545219, 0.0, VerticalDipole()),
        Receiver("Nahanni Butte", 61.0304412, -123.3926734, 0.0, VerticalDipole()),
        Receiver("Juneau", 58.32, -134.41, 0.0, VerticalDipole()),
        Receiver("Ketchikan", 55.35, -131.673, 0.0, VerticalDipole()),
        Receiver("Winnipeg", 49.8822, -97.1308, 0.0, VerticalDipole()),
        Receiver("IslandLake", 53.8626, -94.6658, 0.0, VerticalDipole()),
        Receiver("Gillam", 56.3477, -94.7093, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end

function buildtruepaths()
    transmitters = [
        Transmitter{VerticalDipole}(TRANSMITTER[:NLK].name, TRANSMITTER[:NLK].latitude, TRANSMITTER[:NLK].longitude, TRANSMITTER[:NLK].antenna, TRANSMITTER[:NLK].frequency,240e3),    
        Transmitter{VerticalDipole}(TRANSMITTER[:NML].name, TRANSMITTER[:NML].latitude, TRANSMITTER[:NML].longitude, TRANSMITTER[:NML].antenna, TRANSMITTER[:NML].frequency, 235e3)
    ]
    receivers = [
        Receiver("Whitehorse", 60.724, -135.043, 0.0, VerticalDipole()),
        Receiver("Churchill", 58.74, -94.085, 0.0, VerticalDipole()),
        Receiver("Stony Rapids", 59.253, -105.834, 0.0, VerticalDipole()),
        Receiver("Fort Smith", 60.006, -111.92, 0.0, VerticalDipole()),
        Receiver("Bella Bella", 52.1675508, -128.1545219, 0.0, VerticalDipole()),
        Receiver("Nahanni Butte", 61.0304412, -123.3926734, 0.0, VerticalDipole()),
        Receiver("Juneau", 58.32, -134.41, 0.0, VerticalDipole()),
        Receiver("Ketchikan", 55.35, -131.673, 0.0, VerticalDipole()),
        Receiver("Winnipeg", 49.8822, -97.1308, 0.0, VerticalDipole()),
        Receiver("IslandLake", 53.8626, -94.6658, 0.0, VerticalDipole()),
        Receiver("Gillam", 56.3477, -94.7093, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end


"""
    buildAVIDpaths()

Return a vector of `(Transmitter, Receiver)` propagation paths from the Array for VLF Imaging of the D-region (AVID).
"""
function buildAVIDpaths()
    transmitters = [
        Transmitter{VerticalDipole}(TRANSMITTER[:NLK].name, TRANSMITTER[:NLK].latitude, TRANSMITTER[:NLK].longitude, TRANSMITTER[:NLK].antenna, TRANSMITTER[:NLK].frequency,240e3),    
        Transmitter{VerticalDipole}(TRANSMITTER[:NML].name, TRANSMITTER[:NML].latitude, TRANSMITTER[:NML].longitude, TRANSMITTER[:NML].antenna, TRANSMITTER[:NML].frequency, 235e3)
    ]
    receivers = [
        Receiver("WHI", 60.751, -135.101, 0.0, VerticalDipole()),
        Receiver("CH", 58.763, -94.080, 0.0, VerticalDipole()),
        Receiver("RAB", 58.222, -103.680, 0.0, VerticalDipole()),
        Receiver("FSM", 60.026, -111.931, 0.0, VerticalDipole()),
        Receiver("CAL", 51.654, -128.133, 0.0, VerticalDipole()),
        Receiver("FSI", 61.757, -121.229, 0.0, VerticalDipole()),
        Receiver("JU", 58.590, -134.904, 0.0, VerticalDipole()),
        Receiver("PRU", 54.310, -130.326, 0.0, VerticalDipole()),
        Receiver("PIN", 50.258, -95.865, 0.0, VerticalDipole()),
        Receiver("ISL", 53.855, -94.660, 0.0, VerticalDipole()),
        Receiver("GIL", 56.377, -94.644, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end

"""
    buildAVIDpluspaths()

Return a vector of `(Transmitter, Receiver)` propagation paths from AVID including an AARDVARK receiver near Edmonton.
"""
function buildAVIDpluspaths()
    transmitters = [
        Transmitter{VerticalDipole}(TRANSMITTER[:NLK].name, TRANSMITTER[:NLK].latitude, TRANSMITTER[:NLK].longitude, TRANSMITTER[:NLK].antenna, TRANSMITTER[:NLK].frequency,240e3),    
        Transmitter{VerticalDipole}(TRANSMITTER[:NML].name, TRANSMITTER[:NML].latitude, TRANSMITTER[:NML].longitude, TRANSMITTER[:NML].antenna, TRANSMITTER[:NML].frequency, 235e3)
    ]
    receivers = [
        Receiver("WHI", 60.751, -135.101, 0.0, VerticalDipole()),
        Receiver("CH", 58.763, -94.080, 0.0, VerticalDipole()),
        Receiver("RAB", 58.222, -103.680, 0.0, VerticalDipole()),
        Receiver("FSM", 60.026, -111.931, 0.0, VerticalDipole()),
        Receiver("CAL", 51.654, -128.133, 0.0, VerticalDipole()),
        Receiver("FSI", 61.757, -121.229, 0.0, VerticalDipole()),
        Receiver("JU", 58.590, -134.904, 0.0, VerticalDipole()),
        Receiver("PRU", 54.310, -130.326, 0.0, VerticalDipole()),
        Receiver("PIN", 50.258, -95.865, 0.0, VerticalDipole()),
        Receiver("ISL", 53.855, -94.660, 0.0, VerticalDipole()),
        Receiver("GIL", 56.377, -94.644, 0.0, VerticalDipole()),
        Receiver("ED", 53.147, -113.343, 0.0, VerticalDipole())
    ]
    paths = [(tx, rx) for tx in transmitters for rx in receivers]

    return paths
end

function txpower(paths, txname::AbstractString)
    for (tx, _) in paths
        tx.name == txname && return tx.power
    end
    error("No transmitter found with name = $txname")
end

#TODO: Why is rng not defaulted to reset_rng()?
function sample_values(A::KeyedArray, n::Int; rng=Random.default_rng())
    data = parent(A)              # raw Array
    idx = sample(rng, eachindex(data), n, replace=false)
    return data[idx]
end

"""
observations(name, σamp, σphase)

From filename "name" read in a JLD2 file of amp and phase measurments and randomly add noise from provided σamp, σphase.
"""
function observations(name, σamp, σphase; paths=buildpaths(), rng=reset_rng())

    obs_fname = name*".jld2"
    f = jldopen(obs_fname, "r")
    obsamp, obsphase = f["obsamp"], f["obsphase"]

    npaths = length(paths)
    data = KeyedArray(Array{Float64,3}(undef, 4, npaths, DATALENGTH);
        field=[:amp, :phase, :amp_noiseless, :phase_noiseless], path=MVIA.pathname.(paths), t=1:DATALENGTH)
    data(:amp_noiseless) .= obsamp
    data(:phase_noiseless) .= obsphase
    data(:amp) .= obsamp .+ σamp.*randn(rng, npaths, DATALENGTH)
    data(:phase) .= obsphase .+ σphase.*randn(rng, npaths, DATALENGTH)

    return data
end

"""
observations(name, σamp, σphase)

From filename "name" read in a JLD2 file of amp and phase measurments and randomly add noise from provided σamp, σphase.
"""
function observations(name)
    
    obs_fname = name*".jld2"
    f = jldopen(obs_fname, "r")
    data = f["params"].data.contents

    return data
end

function init_params()
    ### Read parameters from environment variables, with defaults if not set.
    new_folder = get(ENV, "NEW_FOLDER", "false")

    ens_size = parse(Int, get(ENV, "ENS_SIZE", "4"))
    ntimes = parse(Int, get(ENV, "ITRS", "1"))
    shuffle_xy = parse(Bool, get(ENV, "SHUFFLE_XY", "false"))
    ρ = parse(Float64, get(ENV, "RHO","1.1"))
    xy_file = get(ENV, "XY_FILE","false")
    rng = reset_rng()

    statetypes = ()
    if parse(Bool, get(ENV, "DO_XY", "true")) # Currently, only "true" is supported
        statetypes = (statetypes..., :xy)
    end

    datatypes = ()
    if parse(Bool, get(ENV, "DO_AMP", "true")) 
        datatypes = (datatypes..., :amp)
    end
    if parse(Bool, get(ENV, "DO_PHASE", "true"))
        datatypes = (datatypes..., :phase)
    end

    ### Determine which data file to use based on time of day and exact array configuration.
    timeofday = get(ENV,"TOD","day")
    pathset = get(ENV, "PATHS", "AVID")

    if timeofday == "day" && pathset == "Standard"
        dt = DateTime(2020, 3, 1, 20, 00) #using Dates
    elseif timeofday == "day" && (pathset == "AVID" || pathset == "AVIDPLUS")
        dt = DateTime(2025, 6, 1, 19, 00)
    elseif timeofday == "morning" && (pathset == "AVID" || pathset == "AVIDPLUS")
        dt = DateTime(2025, 6, 1, 13, 00)
    elseif timeofday == "night" && (pathset == "AVID" || pathset == "AVIDPLUS")
        dt = DateTime(2025, 6, 2, 7)
    else
        error("Unknown time and pathset combination: ", timeofday, " ", pathset)
    end

    epp = get(ENV, "EPP", "none")

    if pathset == "Standard"
        paths = buildpaths()
        datafile = "Inputs/"*timeofday*"1"
    elseif pathset == "AVID"
        paths = buildAVIDpaths()
        if epp == "none"
            datafile = "Inputs/"*timeofday*"1_buildAVIDpaths"
        else
            datafile = "Inputs/"*timeofday*"1_"*epp
        end
    elseif pathset == "AVIDPLUS"
        paths = buildAVIDpluspaths()
        datafile = "Inputs/"*timeofday*"1_buildAVIDpluspaths"
    else
        error("Unknown Pathset: ", pathset)
    end

    ### Construct observations data structure, and related parameters.
    σamp, σphase = 0.1, deg2rad(1.0)

    npaths = length(paths)

    R = Float64[]
    if :amp in datatypes
        R = [R; fill(σamp^2, npaths)]
    end
    if :phase in datatypes
        R = [R; fill(σphase^2, npaths)]
    end

    ### Create geospatial grid and related parameters
    pathstep = 50e3 # defined in WGS84

    @unpack west, east, south, north, modelproj = common_simulation()

    # Setup Grid
    dr = parse(Float64, get(ENV, "DR_KM", "200")) * 1e3
    lengthscale = 600e3
    modelsteps = ((;dr, lengthscale),)

    x_grid, y_grid = MVIA.build_xygrid(west, east, south, north, MVIA.wgs84(), modelproj; dr) #build_xygrid() & wgs84() defined in MVIA

    trans = Proj.Transformation(modelproj, MVIA.wgs84())
    lola = trans.(parent(parent(MVIA.densify(x_grid, y_grid))))

    localization = MVIA.obs2grid_distance(lola, paths; r=lengthscale)
    localization_mask = get(ENV, "LOCALIZATION_MASK", "KRIGING")
    if localization_mask == "RECTANGLE"
        filterbounds!(localization, lola, west, east, south, north)
    elseif localization_mask =="KRIGING"
        varmap = MVIA.krigingmask(paths, modelproj, x_grid, y_grid; pathstep=pathstep)
        filterbounds!(localization, varmap, 0.2^2)
    else
        error("Unknown localization mask: "*localization_mask)
    end

    # Build interpolant (needed for plots)
    itppts = MVIA.build_xygrid(MVIA.anylocal(localization), x_grid, y_grid)
    itp = MVIA.ScatteredInterpolant(SI.GeneralizedPolyharmonic(1,1), modelproj, itppts) #MVIA extends SI.ScatteredInterpolant

    ### Build the initial background ensemble parameters. NOTE: for best results, ensure LMPTools.ferguson() has been updated to add 0.15 to Beta. 
    ncells = length(lola)
    hB = fill(2,ncells) # σ_h′
    bB = fill(0.04, ncells) # σ_β

    h_off = parse(Float64, get(ENV, "H_OFF", "0"))
    b_off = parse(Float64, get(ENV, "B_OFF", "0"))
    #These two commands allow brute force testing of initial background ensemble means. Generally should be left at

    estimator_name = get(ENV, "WAIT_ESTIMATOR", "THOMSON")
    if     estimator_name == "FIRI"
        back_estimator = FIRIFit()
    elseif estimator_name == "MCRAETHOMSON" || estimator_name == "THOMSON"
        back_estimator = McRaeThomson()
    elseif estimator_name == "FERGUSON"
        back_estimator = Ferguson()
    else
        throw(ArgumentError(
            "Unknown WAIT_ESTIMATOR='$estimator_name'. " *
            "Options: FERGUSON (default), FIRI, MCRAETHOMSON (alias THOMSON)."))
    end

    hb0 = [hprime_beta(back_estimator, ll, dt) for ll in lola]
    h0 = getindex.(hb0, 1) .+ h_off
    b0 = getindex.(hb0, 2) .+ b_off

    @assert length(h0) == length(hB) == ncells

    return(;new_folder, ens_size, ntimes, shuffle_xy, ρ, xy_file, rng, statetypes, datatypes, 
    timeofday, pathset, dt, epp, paths, datafile, σamp, σphase, R, pathstep, modelsteps, 
    x_grid, y_grid, localization, localization_mask, itp, hB, bB, h_off, b_off, estimator_name, h0, b0)
end

function init_rx_params(p)
    statetypes = (p.statetypes..., :rx)
    rx_scenario = parse(Int, get(ENV, "RX_SCENARIO", "0"))

    # Paths whose normalized posterior max exceeds this threshold use the MAP Bϕ
    # deterministically (all members same offset) instead of per-member sampling.
    # Collapses phase-ensemble spread once a path is confident, stabilizing the
    # xy_state update. Set to 1.0 to always sample.
    rx_commit_threshold = parse(Float64, get(ENV, "RX_COMMIT_THRESHOLD", "1.0"))

    # Optional learning rate that tempers per-iteration evidence to prevent
    # overconfidence in early iterations.
    η = parse(Float64, get(ENV, "RX_TEMPER", "1.0"))
    (0 < η ≤ 1) || error("RX_TEMPER (η) must be in (0, 1]; got $η")

    if rx_scenario == 0
        rx_offsets = rx_offsets0(p.paths)
    elseif rx_scenario == 1
        rx_offsets = rx_offsets1(p.paths)
    else
        error("Unknown RX_SCENARIO: ", rx_scenario)
    end
    return merge(p, (; statetypes, rx_scenario, rx_commit_threshold, η, rx_offsets))
end

function init_tx_params(p)
    statetypes = (p.statetypes..., :tx)
    shuffle_tx = parse(Bool, get(ENV, "SHUFFLE_TX", "false"))

    NLKb = parse(Float64, get(ENV, "NLKB","250")) * 1000
    NMLb = parse(Float64, get(ENV, "NLKB","233")) * 1000

    σTX = parse(Float64, get(ENV, "STDDEV_TX_KW", "50")) * 1000 #convert to watts
    σNLK, σNML = σTX, σTX

    TX_range = parse(Float64, get(ENV, "TX_RANGE_KW", "500")) * 1000 #convert to watts

    σTXkw = Int(σTX/1000)
    return merge(p, (; statetypes, shuffle_tx, NLKb, NMLb, σNLK, σNML, TX_range, σTXkw))
end

function init_filter_params(p)
    split_ens_size = parse(Int, get(ENV, "SPLIT_ENS_SIZE", "-1")) # Default value means use ens_size
    if split_ens_size == -1
        split_ens_size = p.ens_size
    end
    split_itrs = parse(Int, get(ENV, "SPLIT_ITRS", "1"))
    filtertype = Symbol(get(ENV, "FILTERTYPE", "stacked")) #Options are stacked, dual, or split for use with MVIA
    return merge(p, (; split_ens_size, split_itrs, filtertype))
end

function name_scenario(scenario, parameters)
    @unpack ntimes, ens_size, ρ, statetypes, datatypes, xy_file, new_folder, 
    timeofday, pathset, modelsteps, localization_mask,  h_off, b_off, estimator_name = parameters()
    scenario = scenario * "$(ntimes)itr_$(ens_size)ens_$(ρ)_$(estimator_name)"

    if h_off !=0
        scenario = scenario * "_h$(h_off)"
    end
    if b_off !=0
        intb_off = Int(b_off*100)
        scenario = scenario * "_b$(intb_off)"
    end

    shuffle_check = false
    if :tx in statetypes
        @unpack shuffle_tx, TX_range = parameters()
        if shuffle_tx
            shuffle_check=true
        end

        if TX_range != 500*1000.0
            scenario = scenario * "_tx_constrained_$(TX_range/1000)kW"
        else
            scenario = scenario * "_tx"
        end

        scenario = scenario * "_log10"
    end
    if :rx in statetypes
        @unpack rx_scenario = parameters()
        scenario = scenario * "_rx_scen_$rx_scenario"
        η_rx = get(parameters(), :η, 1.0)
        η_rx == 1.0 || (scenario = scenario * "_eta$(η_rx)")
    end
    if :xy in statetypes
        @unpack shuffle_xy = parameters()
        if shuffle_xy
            shuffle_check=true
        end
    end
    if shuffle_check
        scenario = scenario * "_shuffle"
    end

    if !(:amp in datatypes && :phase in datatypes)
        if :amp in datatypes
            scenario = scenario* "_amp_only"
        elseif :phase in datatypes
            scenario = scenario * "_phase_only"
        end 
    end

    if :tx in statetypes || :rx in statetypes
        @unpack filtertype = parameters()
        if filtertype == :split
            @unpack split_itrs, split_ens_size = parameters()
            scenario = scenario*"_split_$(split_itrs)splititrs_$(split_ens_size)splitens"
        elseif filtertype == :dual
            scenario = scenario*"_dual"
        else
            scenario = scenario*"_stacked"
        end
    end

    if xy_file != "false"
        file_t = parse(Int, get(ENV, "FILE_T","-1"))
        scenario = scenario * "_xy_file_$file_t"
        @info "Background Ionosphere from file"
    end

    scenario = scenario * "_" * timeofday*"1"*pathset

    scenario = scenario * "_" * localization_mask[1]

    if modelsteps[1].dr != 300 * 1e3
        scenario = scenario*"_$(Int(modelsteps[1].dr /1e3))dr"
    end

    @info "Scenario: "*scenario
    return scenario
end

function rx_offsets0(paths)
    return [0.0 for _ in paths]
end

function rx_offsets1(paths)
    offsets = Float64[]
    for l in 0:length(paths)-1
        push!(offsets, mod(l, 4))
    end
    offsets
end

"""
    common_simulation()

Return parameters common to all scenarios.
"""
function common_simulation()
    modelproj = MVIA.esri_102010() # from SubionosphericInversionAlgorithms.jl

    # Simulation boundaries
    west, east = -136.5, -91
    south, north = 46, 64

    return (;modelproj, west, east, south, north)
end


function window_start(x::Int, i::Int, w::Int)
    @assert 1 ≤ i ≤ x "Index i out of bounds"
    @assert 1 ≤ w ≤ x "Window width w invalid"

    half = (w - 1) ÷ 2
    start = i - half

    # Clamp to valid range while preserving window width
    start = max(1, min(start, x - w + 1))

    return start
end


"""
    centered_window(i, split_itrs, ntimes) → UnitRange{Int}
 
Return a range of exactly `split_itrs` observation indices centered on `i`,
clamped to `[1, ntimes]`.  At the boundaries the window is **shifted** (not
truncated) so it always contains exactly `split_itrs` steps.
 
# Examples
```julia
centered_window(5,  5, 20)   # → 3:7   (symmetric)
centered_window(1,  5, 20)   # → 1:5   (shifted right at left boundary)
centered_window(20, 5, 20)   # → 16:20 (shifted left at right boundary)
centered_window(3,  7, 20)   # → 1:7   (left-edge shift)
```
"""
function centered_window(i::Int, split_itrs::Int, ntimes::Int)
    half = split_itrs ÷ 2
    lo   = i - half
    hi   = lo + split_itrs - 1
 
    # Shift window to remain inside [1, ntimes]
    if lo < 1
        lo = 1
        hi = split_itrs
    elseif hi > ntimes
        hi = ntimes
        lo = ntimes - split_itrs + 1
    end
 
    # Hard clamp in case ntimes < split_itrs
    lo = max(lo, 1)
    hi = min(hi, ntimes)
 
    return lo:hi
end

"""
    first_t_with_agreement(A::KeyedArray, frac; atol=0.0)

Return the first `t` key such that for every path, at least `frac`
(0–1) of ensemble values at that t are identical.

If none exists, return `nothing`.

`atol` allows approximate matching for floating-point values.
"""
function first_t_with_agreement(A, frac; atol=0.0)
    data = parent(A)
    npath, nens, nt = size(data)

    tkeys = axiskeys(A, :t)

    frac = float(frac)  # allow Int input like 1

    # maximum multiplicity fraction in vector v
    function agreement_fraction(v)
        counts = Dict{Float64,Int}()
        best = 0

        @inbounds for x in v
            key = atol == 0 ? x : round(x / atol) * atol
            c = get!(counts, key, 0) + 1
            counts[key] = c
            best = max(best, c)
        end

        return best / length(v)
    end

    @inbounds for ti in 1:nt
        for p in 1:npath
            v = @view data[p, :, ti]

            if agreement_fraction(v) < frac
                @goto next_t
            end
        end

        return tkeys[ti]

        @label next_t
    end

    return nothing
end