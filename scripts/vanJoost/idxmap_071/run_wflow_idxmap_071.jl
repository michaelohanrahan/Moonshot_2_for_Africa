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
mapping this via the indices to the Wflow model grid. Optionally, a lapse rate correction can be applied, to
correct the temperature field with the high resolution model DEM.
"""
function run_timestep_idxmapping(model::Wflow.Model, correct2sea, correct2dem, lon_idx, lat_idx, indices; update_func=Wflow.update, write_model_output=true)
    Wflow.advance!(model.clock)
    # load_dynamic_input!(model)
    load_dynamic_input_idxmapping!(model, correct2sea, correct2dem, lon_idx, lat_idx, indices)
    model = update_func(model)
    if write_model_output
        Wflow.write_output(model)
    end
    return model
end

"""
Custom run function of Wflow, adapted from v0.7.1, that includes index mapping coarser resolution
input data to the Wflow grid.
"""
function run_custom_idxmapping(model::Wflow.Model; close_files=true)
    @unpack network, config, writer, clock = model

    model_type = config.model.type::String

    # determine timesteps to run
    calendar = get(config, "calendar", "standard")::String
    @warn string(
        "The definition of `starttime` has changed (equal to model state time).\n Please",
        " update your settings TOML file by subtracting one model timestep Δt from the",
        " `starttime`, if it was used with a Wflow version up to v0.6.3.",
    )
    starttime = clock.time
    Δt = clock.Δt
    endtime = Wflow.cftime(config.endtime, calendar)
    times = range(starttime + Δt, endtime, step=Δt)

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
        # orography = Wflow.read_standardized(NCDataset(abspath_orography), forcing_name, (x=:, y=:))
        orography = NCDataset(abspath_orography)[forcing_name]
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

    # get index mappings
    path_idx = Wflow.input_path(config, get(config.forcing, "path_idx", ""))
    ds_idx = NCDataset(path_idx)
    lon_idx = Wflow.read_standardized(ds_idx, get(config.forcing, "lon_idx_name", "lon_idx"), (x=:, y=:))
    lat_idx = Wflow.read_standardized(ds_idx, get(config.forcing, "lat_idx_name", "lat_idx"), (x=:, y=:))
    indices = LinearIndices(size(lat_idx))

    @info "Run information" model_type starttime Δt endtime nthreads()
    runstart_time = now()
    @progress for (i, time) in enumerate(times)
        @debug "Starting timestep." time i now()
        model = run_timestep_idxmapping(model, correct2sea, correct2dem, lon_idx, lat_idx, indices)
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
        Wflow.close_files(model, delete_output=false)
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
Function to load forcing and map to the Wflow grid on the fly
"""
function load_dynamic_input_idxmapping!(model, correct2sea, correct2dem, lon_idx, lat_idx, indices)
    update_forcing_idxmapping!(model, correct2sea, correct2dem, lon_idx, lat_idx, indices)
    if haskey(model.config.input, "cyclic")
        Wflow.update_cyclic!(model)
    end
end

"""
Read forcing and map via indices to model resolution
"""
function update_forcing_idxmapping!(model, correct2sea, correct2dem, lon_idx, lat_idx, indices)
    @unpack vertical, clock, reader, network, config = model
    @unpack dataset, dataset_times, forcing_parameters = reader

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
    for (par, ncvar) in forcing_parameters
        # no need to update fixed values
        ncvar.name === nothing && continue

        time = convert(eltype(dataset_times), clock.time)
        t_index = findfirst(>=(time), dataset_times)
        time < first(dataset_times) && throw(DomainError("time $time before dataset begin $(first(dataset_times))"))
        t_index === nothing && throw(DomainError("time $time after dataset end $(last(dataset_times))"))

        # data = Wflow.get_at(dataset, ncvar.name, dataset_times, time)
        data = dataset[ncvar.name][:, :, t_index]

        if ncvar.scale != 1.0 || ncvar.offset != 0.0
            data .= data .* ncvar.scale .+ ncvar.offset
        end

        if par[2] == :temperature
            if lapse_correction

                data = data - correct2sea
                data = forcing_to_wflow_grid(data, lon_idx, lat_idx, indices)
                data = data + correct2dem

            else
                data = forcing_to_wflow_grid(data, lon_idx, lat_idx, indices)
            end
        else
            data = forcing_to_wflow_grid(data, lon_idx, lat_idx, indices)
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


function forcing_to_wflow_grid(raw_data, lon_idx, lat_idx, indices)
    # Create a temperorary array to store the forcing values
    result = similar(indices, eltype(raw_data))
    # Fill the data array with the corresponding values from the raw data
    for i in eachindex(indices)
        result[i] = raw_data[lon_idx[indices[i]], lat_idx[indices[i]]]
    end
    # Reshape array into a 2D array similar to the wflow grid
    result_2d = reshape(result, size(lat_idx))
    return result_2d
end

"""
Adapted from v0.7.1, but removed the run(model) option, in order to insert a custom function
"""
function run_custom_idxmapping(config::Wflow.Config)
    modeltype = config.model.type

    model = if modeltype == "sbm"
        Wflow.initialize_sbm_model(config)
    elseif modeltype == "sbm_gwf"
        Wflow.initialize_sbm_gwf_model(config)
    elseif modeltype == "hbv"
        Wflow.initialize_hbv_model(config)
    elseif modeltype == "sediment"
        Wflow.initialize_sediment_model(config)
    elseif modeltype == "flextopo"
        Wflow.initialize_flextopo_model(config)
    else
        error("unknown model type")
    end
    Wflow.load_fixed_forcing(model)
    run_custom_idxmapping(model)
end

"""
Adapted from v0.7.1, but link to the custom_idxmapping function
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
        @info "Wflow version `v0.7.1_custom`"
        # to catch stacktraces in the log file a try-catch is required
        try
            run_custom_idxmapping(config)
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

run(ARGS[1])