reset_rng() = StableRNG(1234)
reset_rng(seed) = StableRNG(seed)

"""
    buildpaths()

Return a vector of `(Transmitter, Receiver)` propagation paths used in the scenarios
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

    pathstep = 50e3 # defined in WGS84

    # Setup Grid
    dr = 300e3
    lengthscale = 600e3
    modelsteps = ((;dr, lengthscale),)

    x_grid, y_grid = MVIA.build_xygrid(west, east, south, north, MVIA.wgs84(), modelproj; dr) #build_xygrid() & wgs84() defined in MVIA

    trans = Proj.Transformation(modelproj, MVIA.wgs84())
    lola = trans.(parent(parent(MVIA.densify(x_grid, y_grid))))

    paths = buildpaths()
    localization = MVIA.obs2grid_distance(lola, paths; r=lengthscale)
    filterbounds!(localization, lola, west, east, south, north)

    # Build interpolant (needed for plots)
    itppts = MVIA.build_xygrid(MVIA.anylocal(localization), x_grid, y_grid)
    itp = MVIA.ScatteredInterpolant(SI.GeneralizedPolyharmonic(1,1), modelproj, itppts) #MVIA extends SI.ScatteredInterpolant

    ncells = length(lola)
    hB = fill(2,ncells) # σ_h′
    bB = fill(0.04, ncells) # σ_β

    hb0 = [MVIA.ferguson(ll[2], LMPTools.zenithangle(ll[2], ll[1], dt), dt) for ll in lola]
    h0 = getindex.(hb0, 1)
    b0 = getindex.(hb0, 2)

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


