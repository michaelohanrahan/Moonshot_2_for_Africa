using Wflow

using Wflow.UnPack
using Wflow.ProgressLogging
using Wflow.NCDatasets
using Wflow.Dates
using Wflow.Base.Threads
using Wflow.Statistics
using Interpolations

"""
Function to run a single timestep, that includes reading coarser resolution forcing data, and
regridding this to the Wflow model grid. Optionally, a lapse rate correction can be applied, to
correct the temperature field with the high resolution model DEM.
"""
function run_timestep_regrid(model::Wflow.Model, correct2sea, correct2dem; update_func=Wflow.update, write_model_output=true)
    Wflow.advance!(model.clock)
    # load_dynamic_input!(model)
    load_dynamic_input_regrid!(model, correct2sea, correct2dem)
    model = update_func(model)
    if write_model_output
        Wflow.write_output(model)
    end
    return model
end

"""
Custom run function of Wflow, adapted from v0.8.0, that includes regridding coarser resolution
input data.
"""
function run_custom_regrid(model::Wflow.Model; close_files=true)
    @unpack network, config, writer, clock = model

    model_type = config.model.type::String

    # determine timesteps to run
    calendar = get(config, "calendar", "standard")::String
    @warn string(
        "The definition of `starttime` has changed (equal to model state time).\n Please",
        " update your settings TOML file by subtracting one model timestep dt from the",
        " `starttime`, if it was used with a Wflow version up to v0.6.3.",
    )
    starttime = clock.time
    dt = clock.dt
    endtime = Wflow.cftime(config.endtime, calendar)
    times = range(starttime + dt, endtime, step=dt)

    #### CUSTOM - get forcing settings
    lapse_correction = get(config.forcing, "lapse_correction", false)::Bool
    @info "Lapse rate correct is set to $lapse_correction"
    if lapse_correction
        lapse_rate = get(config.forcing, "lapse_rate", -0.0065)::Float64
        path_orography = config.forcing.path_orography
        abspath_orography = Wflow.input_path(config, path_orography)
        forcing_name = get(config.forcing, "layer_name", "wflow_dem")::String
        # Read forcing elevation
        @info "Reading forcing dem from `$abspath_orography` with layer `$forcing_name`"
        orography = Wflow.read_standardized(NCDataset(abspath_orography), forcing_name, (x=:, y=:))
        correct2sea = temperature_correction(orography, lapse_rate)
        # Read model elevation
        wflow_name = get(config.input.vertical, "altitude", "wflow_dem")::String
        @info "Reading Wflow DEM from based on layer `$wflow_name` in the staticmaps"
        wflow_dem = Wflow.read_standardized(model.reader.cyclic_dataset, wflow_name, (x=:, y=:))
        # Correct temperature to model elevation
        correct2dem = temperature_correction(wflow_dem, lapse_rate)
    else
        correct2sea = nothing
        correct2dem = nothing
    end


    @info "Run information" model_type starttime dt endtime nthreads()
    runstart_time = now()
    @progress for (i, time) in enumerate(times)
        @debug "Starting timestep." time i now()
        model = run_timestep_regrid(model, correct2sea, correct2dem)
    end
    @info "Simulation duration: $(canonicalize(now() - runstart_time))"

    # write output state NetCDF
    if haskey(config, "state") && haskey(config.state, "path_output")
        @info "Write output states to NetCDF file `$(model.writer.state_nc_path)`."
    end
    Wflow.write_netcdf_timestep(model, writer.state_dataset, writer.state_parameters)

    Wflow.reset_clock!(model.clock, config)

    # option to support running function twice without re-initializing
    # and thus opening the NetCDF files
    if close_files
        custom_close_files(model, delete_output=false)
    end

    # copy TOML to dir_output, to archive what settings were used
    if haskey(config, "dir_output")
        src = normpath(pathof(config))
        dst = Wflow.output_path(config, basename(src))
        if src != dst
            @debug "Copying TOML file." src dst
            cp(src, dst, force=true)
        end
    end
    return model
end


"""
Function to calculate the temperature correction
"""
function temperature_correction(dem, lapse_rate)
    return dem * lapse_rate
end

"""
Function to load forcing and regrid on the fly
"""
function load_dynamic_input_regrid!(model, correct2sea, correct2dem)
    update_forcing_regrid!(model, correct2sea, correct2dem)
    if haskey(model.config.input, "cyclic")
        Wflow.update_cyclic!(model)
    end
end

"""
Read forcing and regrid to model resolution
"""
function update_forcing_regrid!(model, correct2sea, correct2dem)
    @unpack vertical, clock, reader, network, config = model
    @unpack datasets, dataset_times, forcing_parameters = reader

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool

    if do_reservoirs
        sel_reservoirs = network.reservoir.indices_coverage
        param_res = Wflow.get_param_res(model)
    end
    if do_lakes
        sel_lakes = network.lake.indices_coverage
        param_lake = Wflow.get_param_lake(model)
    end

    # get forcing settings
    lapse_correction = get(config.forcing, "lapse_correction", false)::Bool

    # Wflow expects `right` labeling of the forcing time interval, e.g. daily precipitation
    # at 01-02-2000 00:00:00 is the accumulated total precipitation between 01-01-2000
    # 00:00:00 and 01-02-2000 00:00:00.

    # load from NetCDF into the model according to the mapping
    for ((par, ncvar), dataset) in zip(forcing_parameters, datasets)
        # no need to update fixed values
        ncvar.name === nothing && continue

        time = convert(eltype(dataset_times), clock.time)
        data = Wflow.get_at(dataset, ncvar.name, dataset_times, time)

        if ncvar.scale != 1.0 || ncvar.offset != 0.0
            data .= data .* ncvar.scale .+ ncvar.offset
        end

        xy_orig = [
            Wflow.read_x_axis(dataset),
            Wflow.read_y_axis(dataset)
        ]
        xy_target = [
            Wflow.read_x_axis(model.reader.cyclic_dataset),
            Wflow.read_y_axis(model.reader.cyclic_dataset),
        ]

        if par[2] == :temperature
            if lapse_correction

                data = regrid_data(
                    data,
                    xy_orig,
                    xy_target,
                    correct2sea,
                    correct2dem)

            else
                data = regrid_data(data, xy_orig, xy_target)
            end
        else
            data = regrid_data(data, xy_orig, xy_target)
        end

        # calculate the mean precipitation and evaporation over the lakes and reservoirs
        # and put these into the lakes and reservoirs structs
        # and set the precipitation and evaporation to 0 in the vertical model
        if par in Wflow.mover_params
            if do_reservoirs
                for (i, sel_reservoir) in enumerate(sel_reservoirs)
                    avg = mean(data[sel_reservoir])
                    data[sel_reservoir] .= 0
                    param_res[par][i] = avg
                end
            end
            if do_lakes
                for (i, sel_lake) in enumerate(sel_lakes)
                    avg = mean(data[sel_lake])
                    data[sel_lake] .= 0
                    param_lake[par][i] = avg
                end
            end
        end

        param_vector = Wflow.param(model, par)
        sel = Wflow.active_indices(network, par)
        data_sel = data[sel]
        if any(ismissing, data_sel)
            print(par)
            msg = "Forcing data has missing values on active model cells for $(ncvar.name)"
            throw(ArgumentError(msg))
        end
        param_vector .= data_sel
    end

    return model
end


"""
Regrid forcing to model grid, without lapse rate correction
"""
function regrid_data(raw, x_y_orig, x_y_target)
    # Get x-y arrays
    x_orig, y_orig = x_y_orig
    x_target, y_target = x_y_target

    # Create an interpolation function from the temperature data
    interp_func = interpolate((x_orig, y_orig), raw, Gridded(Constant{Nearest}()))

    # Interpolate the temperature to the entire target grid in one step
    target = interp_func[x_target, y_target]
    return target
end

"""
Regrid forcing to model grid, with lapse rate correction
"""
function regrid_data(
    raw,
    x_y_orig,
    x_y_target,
    correct2sea,
    correct2dem)

    raw2sea = raw - correct2sea

    x_orig, y_orig = x_y_orig
    x_target, y_target = x_y_target

    # Create an interpolation function from the temperature data
    interp_func = interpolate((x_orig, y_orig), raw2sea, Gridded(Constant{Nearest}()))

    # Interpolate the temperature to the entire target grid in one step
    target = interp_func[x_target, y_target]

    target_out = target + correct2dem

    return target_out
end

"""
Adapted from v0.8.0, but removed the run(model) option, in order to insert a custom function
"""
function run_custom_regrid(config::Wflow.Config)
    modeltype = config.model.type

    model = if modeltype == "sbm"
        custom_initialize_sbm_model(config)
    else
        error("this script is only compatibel with model type sbm")
    end
    Wflow.load_fixed_forcing(model)
    run_custom_regrid(model)
end

"""
Adapted from v0.8.0, but link to the custom_regrid function
"""
function run(tomlpath::AbstractString; silent=nothing)
    config = Wflow.Config(tomlpath)
    # if the silent kwarg is not set, check if it is set in the TOML
    if silent === nothing
        silent = get(config, "silent", false)::Bool
    end
    fews_run = get(config, "fews_run", false)::Bool
    logger, logfile = Wflow.init_logger(config; silent)
    Wflow.with_logger(logger) do
        @info "Wflow version `v0.8.0_custom`"
        # to catch stacktraces in the log file a try-catch is required
        try
            run_custom_regrid(config)
        catch e
            # avoid logging backtrace for the single line FEWS log format
            # that logger also uses SimpleLogger which doesn't result in a good backtrace
            if fews_run
                @error "Wflow simulation failed" exception = e _id = :wflow_run
            else
                @error "Wflow simulation failed" exception = (e, catch_backtrace()) _id =
                    :wflow_run
            end
            rethrow()
        finally
            close(logfile)
        end
    end
end


function custom_read_dataset(path_forcing, config)
    abspath_forcing = Wflow.input_path(config, path_forcing)

    # absolute paths are not supported, see Glob.jl#2
    # the path separator in a glob pattern is always /
    if isabspath(path_forcing)
        parts = splitpath(path_forcing)
        # use the root/drive as the dir, to support * in directory names as well
        glob_dir = parts[1]
        glob_path = join(parts[2:end], '/')
    else
        tomldir = dirname(config)
        dir_input = get(config, "dir_input", ".")
        glob_dir = normpath(tomldir, dir_input)
        glob_path = replace(path_forcing, '\\' => '/')
    end
    @info "Reading `$abspath_forcing` for forcing parameters."
    
    dynamic_paths = Wflow.glob(glob_path, glob_dir)  # expand "data/forcing-year-*.nc"
    if isempty(dynamic_paths)
        error("No files found with name '$glob_path' in '$glob_dir'")
    end
    return dynamic_paths
end


function custom_prepare_reader(config)
    paths_forcing = config.input.paths_forcing
    datasets = []
    # alphabetical order so first pet, then precip and then temp
    for path_forcing in paths_forcing
        dynamic_paths = custom_read_dataset(path_forcing, config)
        dataset = NCDataset(dynamic_paths, aggdim = "time", deferopen = false)
        push!(datasets, dataset)
    end

    cyclic_path = Wflow.input_path(config, config.input.path_static)
    @info "Cyclic parameters are provided by `$cyclic_path`."

    dataset = datasets[1]
    if haskey(dataset["time"].attrib, "_FillValue")
        @warn "Time dimension contains `_FillValue` attribute, this is not in line with CF conventions."
        nctimes = dataset["time"][:]
        times_dropped = collect(skipmissing(nctimes))
        # check if lenght has changed (missings in time dimension are not allowed), and throw
        # an error if the lenghts are different
        if length(times_dropped) != length(nctimes)
            error("Time dimension in `$abspath_forcing` contains missing values")
        else
            nctimes = times_dropped
            nctimes_type = eltype(nctimes)
        end
    else
        nctimes = dataset["time"][:]
        nctimes_type = eltype(nctimes)
    end

    # check for cyclic parameters
    do_cyclic = haskey(config.input, "cyclic")

    # create map from internal location to netCDF variable name for forcing parameters
    forcing_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
    for par in config.input.forcing
        fields = Wflow.symbols(par)
        ncname, mod = Wflow.ncvar_name_modifier(Wflow.param(config.input, fields))
        forcing_parameters[fields] =
            (name = ncname, scale = mod.scale, offset = mod.offset, value = mod.value)

        @info "Set `$par` using netCDF variable `$ncname` as forcing parameter."
    end

    # create map from internal location to netCDF variable name for cyclic parameters and
    # store cyclic times for each internal location (duplicate cyclic times are possible
    # this way, however it seems not worth to keep track of unique cyclic times for now
    # (memory usage))
    if do_cyclic == true
        cyclic_dataset = NCDataset(cyclic_path)
        cyclic_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
        cyclic_times = Dict{Tuple{Symbol,Vararg{Symbol}},Vector{Tuple{Int,Int}}}()
        for par in config.input.cyclic
            fields = Wflow.symbols(par)
            ncname, mod = Wflow.ncvar_name_modifier(Wflow.param(config.input, fields))
            i = findfirst(x -> startswith(x, "time"), dimnames(cyclic_dataset[ncname]))
            dimname = Wflow.dimnames(cyclic_dataset[ncname])[i]
            cyclic_nc_times = collect(cyclic_dataset[dimname])
            cyclic_times[fields] = Wflow.timecycles(cyclic_nc_times)
            cyclic_parameters[fields] =
                (name = ncname, scale = mod.scale, offset = mod.offset)

            @info "Set `$par` using netCDF variable `$ncname` as cyclic parameter, with `$(length(cyclic_nc_times))` timesteps."
        end
    else
        cyclic_parameters = Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}()
        cyclic_dataset = nothing
        cyclic_times = Dict{Tuple{Symbol,Vararg{Symbol}},Vector{Tuple{Int,Int}}}()
    end
    # check if there is overlap
    if do_cyclic == true
        forcing_keys = keys(forcing_parameters)
        cyclic_keys = keys(cyclic_parameters)
        for forcing_key in forcing_keys
            if forcing_key in cyclic_keys
                error("parameter specified in both forcing and cyclic")
            end
        end
    end

    return NCReader{nctimes_type}(
        datasets,
        nctimes,
        cyclic_dataset,
        cyclic_times,
        forcing_parameters,
        cyclic_parameters,
    )
end

"""
    initialize_sbm_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function custom_initialize_sbm_model(config::Wflow.Config)

    model_type = config.model.type::String
    @info "Initialize model variables for model type `$model_type`."

    # unpack the paths to the netCDF files
    static_path = Wflow.input_path(config, config.input.path_static)
    custom_reader = custom_prepare_reader(config)

    reader = Wflow.prepare_reader(config)
    clock = Wflow.Clock(config, reader)
    dt = clock.dt

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool
    routing_options = ("kinematic-wave", "local-inertial")
    floodplain_1d = get(config.model, "floodplain_1d", false)::Bool
    river_routing = Wflow.get_options(
        config.model,
        "river_routing",
        routing_options,
        "kinematic-wave",
    )::String
    land_routing =
        Wflow.get_options(config.model, "land_routing", routing_options, "kinematic-wave")::String
    do_water_demand = haskey(config.model, "water_demand")

    snow = get(config.model, "snow", false)::Bool
    reservoirs = do_reservoirs
    lakes = do_lakes
    glacier = get(config.model, "glacier", false)::Bool
    masswasting = get(config.model, "masswasting", false)::Bool
    @info "General model settings" reservoirs lakes snow masswasting glacier

    nc = NCDataset(static_path)

    subcatch_2d = Wflow.ncread(nc, config, "subcatchment"; optional = false, allow_missing = true)
    # indices based on catchment
    inds, rev_inds = Wflow.active_indices(subcatch_2d, missing)
    n = length(inds)
    modelsize_2d = size(subcatch_2d)

    river_2d =
        Wflow.ncread(nc, config, "river_location"; optional = false, type = Bool, fill = false)
    river = river_2d[inds]
    riverwidth_2d =
        Wflow.ncread(nc, config, "lateral.river.width"; optional = false, type = Float64, fill = 0)
    riverwidth = riverwidth_2d[inds]
    riverlength_2d =
        Wflow.ncread(nc, config, "lateral.river.length"; optional = false, type = Float64, fill = 0)
    riverlength = riverlength_2d[inds]

    # read x, y coordinates and calculate cell length [m]
    y_nc = Wflow.read_y_axis(nc)
    x_nc = Wflow.read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer = (1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = Wflow.cell_lengths(y, cellength, sizeinmetres)
    riverfrac = Wflow.river_fraction(river, riverlength, riverwidth, xl, yl)

    inds_riv, rev_inds_riv = Wflow.active_indices(river_2d, 0)
    nriv = length(inds_riv)

    sbm = Wflow.initialize_sbm(nc, config, riverfrac, inds)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            Wflow.initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, Wflow.tosecond(dt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            Wflow.initialize_lake(config, nc, inds_riv, nriv, pits, Wflow.tosecond(dt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = Wflow.ncread(nc, config, "ldd"; optional = false, allow_missing = true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = Wflow.ncread(nc, config, "pits"; optional = false, type = Bool, fill = false)
        ldd = Wflow.set_pit_ldd(pits_2d, ldd, inds)
    end

    landslope =
        Wflow.ncread(nc, config, "lateral.land.slope"; optional = false, sel = inds, type = Float64)
    clamp!(landslope, 0.00001, Inf)

    dl = map(Wflow.detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl

    # check if lateral subsurface flow component is defined for the SBM model, when coupled
    # to another groundwater model, this component is not defined in the TOML file.
    subsurface_flow = haskey(config.input.lateral, "subsurface")
    if subsurface_flow
        khfrac = Wflow.ncread(
            nc,
            config,
            "lateral.subsurface.ksathorfrac";
            sel = inds,
            defaults = 1.0,
            type = Float64,
        )

        # unit for lateral subsurface flow component is [m³ d⁻¹], sbm.kv_0 [mm Δt⁻¹]
        kh_0 = khfrac .* sbm.kv_0 .* 0.001 .* (Wflow.basetimestep / dt)
        f = sbm.f .* 1000.0
        zi = sbm.zi .* 0.001
        soilthickness = sbm.soilthickness .* 0.001
        z_exp = sbm.z_exp .* 0.001

        ssf = Wflow.LateralSSF{Float64}(
            kh_0 = kh_0,
            f = f,
            kh = fill(Wflow.mv, n),
            khfrac = khfrac,
            zi = zi,
            z_exp = z_exp,
            soilthickness = soilthickness,
            theta_s = sbm.theta_s,
            theta_r = sbm.theta_r,
            dt = dt / Wflow.basetimestep,
            slope = landslope,
            dl = dl,
            dw = dw,
            exfiltwater = fill(Wflow.mv, n),
            recharge = fill(Wflow.mv, n),
            ssf = fill(Wflow.mv, n),
            ssfin = fill(Wflow.mv, n),
            ssfmax = fill(Wflow.mv, n),
            to_river = zeros(n),
            volume = (sbm.theta_s .- sbm.theta_r) .* (soilthickness .- zi) .* (xl .* yl),
        )
        # update variables `ssf`, `ssfmax` and `kh` (layered profile) based on ksat_profile
        ksat_profile = get(config.input.vertical, "ksat_profile", "exponential")::String
        if ksat_profile == "exponential"
            Wflow.initialize_lateralssf_exp!(ssf::Wflow.LateralSSF)
        elseif ksat_profile == "exponential_constant"
            Wflow.initialize_lateralssf_exp_const!(ssf::Wflow.LateralSSF)
        elseif ksat_profile == "layered" || ksat_profile == "layered_exponential"
            Wflow.initialize_lateralssf_layered!(ssf::Wflow.LateralSSF, sbm::Wflow.SBM, ksat_profile)
        end
    else
        # when the SBM model is coupled (BMI) to a groundwater model, the following
        # variables are expected to be exchanged from the groundwater model.
        ssf = Wflow.GroundwaterExchange{Float64}(
            dt = dt / Wflow.basetimestep,
            exfiltwater = fill(mv, n),
            zi = fill(mv, n),
            to_river = fill(mv, n),
            ssf = zeros(n),
        )
    end

    graph = Wflow.flowgraph(ldd, inds, Wflow.pcr_dir)
    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = Wflow.flowgraph(ldd_riv, inds_riv, Wflow.pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = Wflow.fraction_runoff_toriver(graph, ldd, index_river, landslope, n)

    inds_allocation_areas = Vector{Int}[]
    inds_riv_allocation_areas = Vector{Int}[]
    if do_water_demand
        areas = unique(sbm.allocation.areas)
        for a in areas
            area_index = findall(x -> x == a, sbm.allocation.areas)
            push!(inds_allocation_areas, area_index)
            area_riv_index = findall(x -> x == a, sbm.allocation.areas[index_river])
            push!(inds_riv_allocation_areas, area_riv_index)
        end
    end

    if land_routing == "kinematic-wave"
        olf = Wflow.initialize_surfaceflow_land(
            nc,
            config,
            inds;
            sl = landslope,
            dl,
            width = map(Wflow.det_surfacewidth, dw, riverwidth, river),
            iterate = kinwave_it,
            tstep = kw_land_tstep,
            dt,
        )
    elseif land_routing == "local-inertial"
        index_river_nf = rev_inds_riv[inds] # not filtered (with zeros)
        olf, indices = Wflow.initialize_shallowwater_land(
            nc,
            config,
            inds;
            modelsize_2d,
            indices_reverse = rev_inds,
            xlength = xl,
            ylength = yl,
            riverwidth = riverwidth_2d[inds_riv],
            graph_riv,
            ldd_riv,
            inds_riv,
            river,
            waterbody = !=(0).(resindex + lakeindex),
            dt,
        )
    end

    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")
    if river_routing == "kinematic-wave"
        rf = Wflow.initialize_surfaceflow_river(
            nc,
            config,
            inds_riv;
            dl = riverlength,
            width = riverwidth,
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            iterate = kinwave_it,
            tstep = kw_river_tstep,
            dt = dt,
        )
    elseif river_routing == "local-inertial"
        rf, nodes_at_link = Wflow.initialize_shallowwater_river(
            nc,
            config,
            inds_riv;
            graph = graph_riv,
            ldd = ldd_riv,
            dl = riverlength,
            width = riverwidth,
            reservoir_index = resindex,
            reservoir = reservoirs,
            lake_index = lakeindex,
            lake = lakes,
            dt = dt,
            floodplain = floodplain_1d,
        )
    else
        error(
            """An unknown "river_routing" method is specified in the TOML file ($river_routing).
            This should be "kinematic-wave" or "local-inertial".
            """,
        )
    end

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = Wflow.topological_sort_by_dfs(graph)
    if land_routing == "kinematic-wave" ||
       river_routing == "kinematic-wave" ||
       subsurface_flow
        streamorder = Wflow.stream_order(graph, toposort)
    end
    if land_routing == "kinematic-wave" || subsurface_flow
        toposort = Wflow.topological_sort_by_dfs(graph)
        index_pit_land = findall(x -> x == 5, ldd)
        min_streamorder_land = get(config.model, "min_streamorder_land", 5)
        subbas_order, indices_subbas, topo_subbas = Wflow.kinwave_set_subdomains(
            graph,
            toposort,
            index_pit_land,
            streamorder,
            min_streamorder_land,
        )
    end
    if river_routing == "kinematic-wave"
        min_streamorder_river = get(config.model, "min_streamorder_river", 6)
        toposort_riv = Wflow.topological_sort_by_dfs(graph_riv)
        index_pit_river = findall(x -> x == 5, ldd_riv)
        subriv_order, indices_subriv, topo_subriv = Wflow.kinwave_set_subdomains(
            graph_riv,
            toposort_riv,
            index_pit_river,
            streamorder[index_river],
            min_streamorder_river,
        )
    end

    if nthreads() > 1
        if river_routing == "kinematic-wave"
            @info "Parallel execution of kinematic wave" min_streamorder_land min_streamorder_river
        elseif land_routing == "kinematic-wave" || subsurface_flow
            @info "Parallel execution of kinematic wave" min_streamorder_land
        end
    end

    modelmap = (vertical = sbm, lateral = (subsurface = ssf, land = olf, river = rf))
    indices_reverse = (
        land = rev_inds,
        river = rev_inds_riv,
        reservoir = isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake = isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer = Wflow.prepare_writer(
        config,
        modelmap,
        indices_reverse,
        x_nc,
        y_nc,
        nc,
        extra_dim = (name = "layer", value = Float64.(1:sbm.maxlayers)),
    )
    close(nc)

    # for each domain save:
    # - the directed acyclic graph (graph),
    # - the traversion order (order),
    # - upstream_nodes,
    # - subdomains for the kinematic wave domains for parallel execution (execution order of
    #   subbasins (subdomain_order), traversion order per subbasin (topo_subdomain) and
    #   Vector indices per subbasin matching the traversion order of the complete domain
    #   (indices_subdomain))
    # - the indices that map it back to the two dimensional grid (indices)

    # for the land domain the x and y length [m] of the grid cells are stored
    # for reservoirs and lakes indices information is available from the initialization
    # functions
    land = (
        graph = graph,
        upstream_nodes = Wflow.filter_upsteam_nodes(graph, pits[inds]),
        subdomain_order = subbas_order,
        topo_subdomain = topo_subbas,
        indices_subdomain = indices_subbas,
        order = toposort,
        indices = inds,
        reverse_indices = rev_inds,
        area = xl .* yl,
        slope = landslope,
        indices_allocation_areas = inds_allocation_areas,
    )
    if land_routing == "local-inertial"
        land = merge(land, (index_river = index_river_nf, staggered_indices = indices))
    end
    if do_water_demand
        # exclude waterbodies for local surface and ground water abstraction
        inds_riv_2d = copy(rev_inds_riv)
        inds_2d = ones(Bool, modelsize_2d)
        if !isempty(reservoir)
            inds_cov = collect(Iterators.flatten(reservoir.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        if !isempty(lake)
            inds_cov = collect(Iterators.flatten(lake.indices_coverage))
            inds_riv_2d[inds_cov] .= 0
            inds_2d[inds_cov] .= 0
        end
        land = merge(land, (index_river_wb = inds_riv_2d[inds], index_wb = inds_2d[inds]))
    end
    if river_routing == "kinematic-wave"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for kinematic_wave
            upstream_nodes = Wflow.filter_upsteam_nodes(graph_riv, pits[inds_riv]),
            subdomain_order = subriv_order,
            topo_subdomain = topo_subriv,
            indices_subdomain = indices_subriv,
            order = toposort_riv,
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    elseif river_routing == "local-inertial"
        river = (
            graph = graph_riv,
            indices = inds_riv,
            reverse_indices = rev_inds_riv,
            # reservoir and lake index
            reservoir_index = resindex,
            lake_index = lakeindex,
            reservoir_index_f = filter(x -> x ≠ 0, resindex),
            lake_index_f = filter(x -> x ≠ 0, lakeindex),
            # specific for local-inertial
            nodes_at_link = nodes_at_link,
            links_at_node = Wflow.adjacent_links_at_node(graph_riv, nodes_at_link),
            # water allocation areas
            indices_allocation_areas = inds_riv_allocation_areas,
            area = xl[index_river] .* yl[index_river],
        )
    end

    model = Wflow.Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface = ssf, land = olf, river = rf),
        sbm,
        clock,
        custom_reader,
        writer,
        Wflow.SbmModel(),
    )

    model = Wflow.set_states(model)

    @info "Initialized model with custom script"
    return model
end

struct NCReader{T}
    datasets::Vector{Wflow.CFDataset}
    dataset_times::Vector{T}
    cyclic_dataset::Union{NCDataset,Nothing}
    cyclic_times::Dict{Tuple{Symbol,Vararg{Symbol}},Vector{Tuple{Int,Int}}}
    forcing_parameters::Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}
    cyclic_parameters::Dict{Tuple{Symbol,Vararg{Symbol}},NamedTuple}
end


"Close input and output datasets that are opened on model initialization"
function custom_close_files(model; delete_output::Bool = false)
    @unpack reader, writer, config = model
    for dataset in reader.datasets
        close(dataset)
    end
    if haskey(config.input, "cyclic")
        close(reader.cyclic_dataset)
    end
    writer.dataset === nothing || close(writer.dataset)
    writer.dataset_scalar === nothing || close(writer.dataset_scalar)
    close(writer.csv_io)  # can be an IOBuffer
    writer.state_dataset === nothing || close(writer.state_dataset)

    if delete_output
        if writer.nc_path !== nothing
            isfile(writer.nc_path) && rm(writer.nc_path)
        end
        if writer.csv_path !== nothing
            isfile(writer.csv_path) && rm(writer.csv_path)
        end
        if writer.state_nc_path !== nothing
            isfile(writer.state_nc_path) && rm(writer.state_nc_path)
        end
        if writer.nc_scalar_path !== nothing
            isfile(writer.nc_scalar_path) && rm(writer.nc_scalar_path)
        end
    end
    return nothing
end

run(ARGS[1])