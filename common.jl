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

function txpower(paths, txname::AbstractString)
    for (tx, _) in paths
        tx.name == txname && return tx.power
    end
    error("No transmitter found with name = $txname")
end


function sample_values(A::KeyedArray, n::Int; rng=Random.default_rng())
    data = parent(A)              # raw Array
    idx = sample(rng, eachindex(data), n, replace=false)
    return data[idx]
end

"""
observations(name,σamp, σphase)

From filename "name" read in a JLD2 file of amp and phase measurments and randomly add noise from provided σamp, σphase.
"""
function observations(name, σamp, σphase)
    paths = buildpaths()

    obs_fname = name*".jld2"
    f = jldopen(obs_fname, "r")
    obsamp, obsphase = f["obsamp"], f["obsphase"]

    rng = StableRNG(1234)

    npaths = length(paths)
    data = KeyedArray(Array{Float64,3}(undef, 4, npaths, DATALENGTH);
        field=[:amp, :phase, :amp_noiseless, :phase_noiseless], path=pathname.(paths), t=1:DATALENGTH)
    data(:amp_noiseless) .= obsamp
    data(:phase_noiseless) .= obsphase
    data(:amp) .= obsamp .+ σamp.*randn(rng, npaths, DATALENGTH)
    data(:phase) .= obsphase .+ σphase.*randn(rng, npaths, DATALENGTH)

    return data
end

"""
observations(name,σamp, σphase)

From filename "name" read in a JLD2 file of amp and phase measurments and randomly add noise from provided σamp, σphase.
"""
function observations(name)
    
    obs_fname = name*".jld2"
    f = jldopen(obs_fname, "r")
    data = f["params"].data.contents

    return data
end


"""
    pathname(p)

Return path name string for (transmitter, receiver) path tuple `p`.

    TODO: Exported by SIA, delete when SIA is properly imported
"""
pathname(p) = p[1].name*"-"*p[2].name

function init_params()

    @unpack west, east, south, north, modelproj = common_simulation()

    dt = DateTime(2020, 3, 1, 20, 00) #using Dates
    pathstep = 50e3 # defined in WGS84

    # Setup Grid
    dr = 300e3
    lengthscale = 600e3
    modelsteps = ((;dr, lengthscale),)

    x_grid, y_grid = MVIA.build_xygrid(west, east, south, north, MVIA.wgs84(), modelproj; dr) #build_xygrid() & wgs84() defined in MVIA

    trans = Proj.Transformation(modelproj, MVIA.wgs84())
    lola = trans.(parent(parent(MVIA.densify(x_grid, y_grid))))

    paths = buildtruepaths()
    localization = MVIA.obs2grid_distance(lola, paths; r=lengthscale)
    filterbounds!(localization, lola, west, east, south, north)

    # Build interpolant (needed for plots)
    itppts = MVIA.build_xygrid(MVIA.anylocal(localization), x_grid, y_grid)
    itp = MVIA.ScatteredInterpolant(SI.GeneralizedPolyharmonic(1,1), modelproj, itppts) #MVIA extends SI.ScatteredInterpolant

    ncells = length(lola)
    hB = fill(2,ncells) # σ_h′
    bB = fill(0.04, ncells) # σ_β

    h_off = parse(Float64, get(ENV, "H_OFF", "0"))
    b_off = parse(Float64, get(ENV, "B_OFF", "0"))
    hb0 = [MVIA.ferguson(ll[2], LMPTools.zenithangle(ll[2], ll[1], dt), dt) for ll in lola]
    h0 = getindex.(hb0, 1) .+ h_off
    b0 = getindex.(hb0, 2) .+ b_off

    @assert length(h0) == length(hB) == ncells

    return(;pathstep, modelsteps, x_grid, y_grid, hB, bB, h0, b0, itp, localization)
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
