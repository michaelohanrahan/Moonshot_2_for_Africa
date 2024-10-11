# README
# Summary of the workflow:
# 1) A config is processed as a Config and contains:
# - points of interests (spatial)
# - start date and duration (temporal)
# - model settings such as forcing type
# 2) The Jobs consist of a Forecast per cluster
# 3) Each Forecast requires an initial State
# - Either an existing state is used
# - Or a new state must be created
# 4) Create TOML files and estimate run durations
# - Each Forecast gets its own TOML file
# - Only a new State requires a TOML file

import logging
import datetime
import os
import bisect
import yaml
import glob
import argparse
import traceback
from helper import syscheck
import geopandas as gpd
from hydromt_wflow import WflowModel
from hydromt.config import configread

DRIVE = syscheck()

print(f"{'*'*10}\nDRIVE: {DRIVE}\n{'*'*10}")

logger = logging.getLogger("MS2")

# Set the logging level for the root logger
logging.getLogger().setLevel(logging.INFO)

# Optionally, set logging levels for specific loggers
logging.getLogger('fiona').setLevel(logging.WARNING)
logging.getLogger('geopandas').setLevel(logging.WARNING)
logging.getLogger('rasterio').setLevel(logging.WARNING)
logging.getLogger('gdal').setLevel(logging.WARNING)

_ROOT = os.path.normpath(f"{DRIVE}/moonshot2-casestudy/Wflow/africa/data/")
_TOML_STATE = os.path.normpath(
    f"{DRIVE}/moonshot2-casestudy/Wflow/africa/config/3_warmup_wflow_sbm.toml"
)
_TOML_FORECAST = os.path.normpath(
    f"{DRIVE}/moonshot2-casestudy/Wflow/africa/config/3_forecast_wflow_sbm.toml"
)
_CLUSTERS = os.path.normpath(
    f"{DRIVE}/moonshot2-casestudy/Wflow/africa/data/2-interim/clustered_basins.geojson"
)

_DATE_FORMAT_FNAME = r"%Y%m%d"
_DATE_FORMAT_SHORT = r"%Y-%m-%d"
_DATE_FORMAT_LONG = r"%Y-%m-%dT%H:%M:%S"

_FORCING_FILES = {
    "era5_daily": {
        'precip': ("tp", "p:/wflow_global/hydromt/meteo/era5_daily/tp/era5_tp_*_daily.nc"),
        'temp': ("t2m", "p:/wflow_global/hydromt/meteo/era5_daily/tp/era5_t2m_*_daily.nc"),
        'pet': ("pet", "p:/moonshot2-casestudy/Wflow/africa/data/3-input/global_era5_pet/era5_daily_debruin_PET_daily_*.nc"),
        'temp_in_celsius': False,
    },
    "chirps":  {
        'precip': ("precipitation", "p:/wflow_global/hydromt/meteo/chirps_africa_caily_v2.0/CHIRPS_rainfall*.nc"),
        'temp': ("t2m", "p:/wflow_global/hydromt/meteo/era5_daily/tp/era5_t2m_*_daily.nc"),
        'pet': ("pet", "p:/moonshot2-casestudy/Wflow/africa/data/3-input/global_era5_pet/era5_daily_debruin_PET_daily_*.nc"),
        'temp_in_celsius': False,
    },
    # "era5_hourly": "...", TODO support hourly PET
}

_MAX_RUNTIME = 10  # upper limit for runtime, unit: ms per timestep per km²
# TODO: dict per cluster type / number of cores


def convert_path(path: str) -> str:
    """
    Convert Windows path to Linux path and replace 'p:' or {} with DRIVE.
    """
    
    # self.logger.debug(f"Original path: {path}")
    # self.logger.debug(f"DRIVE value: {DRIVE}")
    
    if os.name != 'nt':  # If not Windows
        path = path.replace('\\', '/')
        if path.startswith('p:/'):
            path = path.replace('p:/', f'{DRIVE}/', 1)
        elif 'p:/' in path:
            path = path.replace('p:/', f'{DRIVE}/')
    
    # Replace {} with DRIVE if present
    path = path.replace('{}', DRIVE)
    
    return path

def time_in_dhms(seconds: float) -> str:
    """
    Converts a given time in seconds to a string representation in
    days, hours, minutes, and seconds.

    Parameters
    ----------
    seconds : float
        The time in seconds to be converted.

    Returns
    -------
    str
        A string representation of the time in the format 'X days, Y hours, Z minutes, W seconds'.
    """
    # Define the time components and their corresponding number of seconds
    components = [("days", 86400), ("hours", 3600), ("minutes", 60), ("seconds", 1)]
    time_values = []

    # Calculate the value for each time component
    for name, secs_per_unit in components:
        value, seconds = divmod(seconds, secs_per_unit)
        if value > 0 or name == "seconds":
            time_values.append(f"{int(value)} {name}")
    return ", ".join(time_values)


class Config:
    """
    A class used to easily access the attributes read from the config file.

    TODO: It might be worh it to replace this Config object by a HydroMT/Hydroflows config?
    """

    def __init__(self, config_path: str, logger=logger) -> None:
        """
        Constructs all the necessary attributes for the Config object.

        Parameters
        ----------
            config_path : str
                The path to the configuration file.
            logger : Logger, optional
                A logger object for logging messages. If not provided, a default logger is used.
        """
        self.logger = logger
        self.logger.info(f"Reading config from {config_path}")
        self.read_config(config_path)

    def read_config(self, config_path: str) -> None:
        """
        Reads the configuration from a YAML file and sets the attributes of the Config object.

        Parameters
        ----------
            config_path : str
                The path to the configuration file.
        """
        with open(config_path, "r") as file:
            config_dict = yaml.safe_load(file)

        self.name = config_dict["name"]
        self.points_of_interest = convert_path(config_dict["points_of_interest"])
        # self.logger.debug(f"Points of interest path: {self.points_of_interest}")
        self.start_date = config_dict["start_date"]
        self.timestepsecs = config_dict["timestepsecs"]
        self.duration_in_timesteps = config_dict["duration_in_timesteps"]
        self.variables = config_dict["variables"]  # not used yet
        self.forcing = config_dict["forcing"]
        self.warmup_forcing = config_dict["warmup_forcing"]

class Jobs:
    """
    A class used to contain the Forecasts for one of more clusters.
    """

    def __init__(self, config_path: str, logger=logger) -> None:
        """
        Constructs all the necessary attributes for the Jobs object.

        Parameters
        ----------
            config : Config
                A Config object containing the configuration for the forecast.
        """
        self.logger = logger
        self.config = Config(config_path, self.logger)
        self.name = self.config.name
        self.logger.info(f"Preparing forecast jobs for {self.name}!")

        self.tstart = datetime.datetime.strptime(
            self.config.start_date, _DATE_FORMAT_SHORT
        )
        self.duration = int(self.config.duration_in_timesteps)
        self.timestepsecs = int(self.config.timestepsecs)
        self.tend = self.tstart + datetime.timedelta(
            seconds=self.timestepsecs * self.config.duration_in_timesteps
        )

        self.forcing = self.config.forcing
        self.warmup_forcing = self.config.warmup_forcing
        for _forcing in [self.forcing, self.warmup_forcing]:
            if _forcing not in _FORCING_FILES.keys():
                raise ValueError(
                    f"Invalid forcing type {self.forcing}, choose from {_FORCING_FILES.keys()}"
                )

        self.points = self.read_points_of_interest()
        self.clusters = gpd.read_file(convert_path(_CLUSTERS))

        self.cluster_ids = []
        self.runtimes = {}
        self.forecasts = {}

    def get_state_file_name(self, forecast, cluster):
        return self.forecasts[cluster].state_file_name

    
    def prepare(self):
        """
        Prepares the Forecast for each cluster.
        """
        self.cluster_ids = self.locate_clusters()
        for cluster_id in self.cluster_ids:
            forecast = Forecast(self, cluster_id)
            forecast.prepare()
            self.forecasts[cluster_id] = forecast

    def locate_clusters(self, column_with_ids: str = "cluster_key") -> list[int]:
        """
        Finds the clusters corresponding to the points of interest using 'geopandas.sjoin'.

        Parameters
        ----------
        column_with_ids : str, optional
            The column in the clusters GeoDataFrame that contains the cluster ids.

        Returns
        -------
        list[int]
            A list of unique cluster ids.
        """
        self.logger.info("Finding clusters corresponding to points of interest")
        cluster_ids = gpd.sjoin(
            self.points, self.clusters, how="inner", predicate="within"
        )[column_with_ids]
        cluster_ids = list(
            map(lambda x: int(x) if str(x).isdigit() else x, set(cluster_ids))
        )
        self.logger.info(f"Found one or more clusters: {cluster_ids}")
        return cluster_ids

    def read_points_of_interest(self):
        try:
            return gpd.read_file(self.config.points_of_interest)
        except Exception as e:
            self.logger.error(f"Error reading points of interest: {e}")
            self.logger.info(f"Attempted to read file: {self.config.points_of_interest}")
            self.logger.info(f"Current working directory: {os.getcwd()}")
            raise


class Run:
    """
    A class used to represent a Wflow model run.

    This class is used as the parent class for the Forecast and State classes.
    A forecast is represented by multiple Run instances.
    """

    def __init__(self, jobs: Jobs, cluster_id: int) -> None:
        """
        Constructs all the necessary attributes for the Run object.

        Parameters
        ----------
            jobs : Jobs
                A Jobs object containing all the jobs for the forecast.
            cluster_id : int
                The cluster id for this specific Wflow model run.
        """
        self.logger = jobs.logger
        self.cluster_id = cluster_id
        self.jobs = jobs

        self.dir_input = os.path.join(_ROOT, "3-input", "wflow_build", str(self.cluster_id))
        self.dir_output = os.path.join(
            _ROOT, "4-output", "wflow_forecast", self.jobs.name, str(self.cluster_id)
        )

        self.forcing = None
        self.toml = None
        self.duration = None
        self.max_runtime = None

        # TOML contents
        self.starttime = None
        self.endtime = None
        self.timestepsecs = None
        self.path_log = None
        self.reinit = None
        self.path_output = None
        self.state_input = None
        self.state_output = None

    def set_toml_forcing(self, forcing_dict=_FORCING_FILES):
        self.var_precip, self.path_precip = forcing_dict[self.forcing]["precip"]
        self.var_pet, self.path_pet = forcing_dict[self.forcing]["pet"]
        self.var_temp, self.path_temp = forcing_dict[self.forcing]["temp"]

    def create_toml(
        self, template: str, hydromt_config_fn: str = None, limit_logging: bool = True
    ) -> None:
        """
        Writes the configuration for the Wflow model to a TOML file.

        This method sets and writes the configuration for the Wflow model using the Run attributes.
        It uses a dummy instance of WflowModel from hydromt_wflow to update and write the TOML file.
        A hydromt_wflow configuration file can be used to further update specific parts of the
        TOML file. The TOML file for both the state and forecast are written to the output directory
        of the forecast.

        Parameters
        ----------
        template : str
            The path to the template TOML file to use for the Wflow model configuration.
        hydromt_config_fn : str, optional
            The path to a hydromt_wflow config file. If provided, the configuration options
            in this file will be read and used to update the Wflow model configuration.
        limit_logging: bool, optional
            Will prevent logging from the hydromt_wflow module to be added to the current log.
        """
        if limit_logging:
            w_logger = logging.getLogger("hydromt_wflow")
            w_logger.setLevel(logging.CRITICAL)
            w_logger.addHandler(logging.NullHandler())
        else:
            w_logger = self.logger

        w = WflowModel(logger=w_logger)
        w.read_config(template)

        w.set_config("starttime", self.starttime.strftime(_DATE_FORMAT_LONG))
        w.set_config("endtime", self.endtime.strftime(_DATE_FORMAT_LONG))
        w.set_config("timestepsec", self.timestepsecs)
        w.set_config("dir_input", self.dir_input)
        w.set_config("dir_output", self.dir_output)
        w.set_config("path_log", self.path_log)

        self.set_toml_forcing()
        w.set_config("input", "path_precip", self.path_precip)
        w.set_config("input", "vertical", "precipitation", self.var_precip)
        w.set_config("input", "path_pet", self.path_pet)
        w.set_config("input", "vertical", "potential_evaporation", self.var_pet)
        w.set_config("input", "path_temp", self.path_temp)
        w.set_config("input", "vertical", "temperature", self.var_temp)

        w.set_config("model", "reinit", self.reinit)
        w.set_config("state", "path_input", self.state_input)
        w.set_config("state", "path_output", self.state_output)
        w.set_config("output", "path", self.path_output)

        if hydromt_config_fn:
            self.logger.info(f"Reading hydromt_wflow config: {hydromt_config_fn}")
            opt = configread(config_fn=hydromt_config_fn)
            w.update(write=False, opt=opt)

        self.logger.info(f"Writing Wflow config file {self.toml} to {self.dir_output}")
        w.write_config(self.toml, self.dir_output)
        w = None  # close model

    def estimate_max_runtime(
        self, clusters: gpd.GeoDataFrame, column_with_ids: str = "cluster_key"
    ) -> float:
        """
        Estimates the maximum runtime for the Wflow model.

        This method calculates the area of the cluster in square kilometers,
        and multiplies it by the duration in timesteps of the model to obtain an
        estimate of the maximum runtime. The global variable _MAX_RUNTIME is
        used which contains a constant for the maximum runtime per km² per timestep in seconds.
        The maximum runtime can be used when the simulation
        is calculated on a computational cluster.

        Parameters
        ----------
        clusters : gpd.GeoDataFrame
            A GeoDataFrame containing the clusters. Each row represents a
            cluster and must have a geometry column with the polygon of the cluster.
        column_with_ids : str, optional
            The name of the column in 'clusters' that contains the cluster IDs.
            The default is 'cluster_key'.

        Returns
        -------
        max_runtime: float
            The maximum runtime of the run in seconds.
        """
        cluster = clusters[clusters[column_with_ids] == self.cluster_id]
        cluster = cluster.to_crs(epsg=3857)
        area_cluster = cluster["geometry"].area / 1e6  # from m² to km²
        max_runtime = (
            area_cluster * self.duration * _MAX_RUNTIME / 1e3
        )  # from ms to seconds
        return max_runtime.values[0]


class State(Run):
    """
    A class used to represent the Run used to create the warm state for a Forecast.

    This class inherits from the Run class and is used in a Forecast object.
    It used the toml defined in the global variable _TOML_STATE as a template for the run settings.
    """

    def __init__(self, jobs: Jobs, cluster_id: int) -> None:
        """
        Constructs all the necessary attributes for the State object.

        Parameters
        ----------
            jobs : Jobs
                A Jobs object containing the jobs for the forecast.
            cluster_id : int
                The cluster id for the Wflow model run.
        """
        super().__init__(jobs, cluster_id)
        self.logger.info(f"Initializing forecast for cluster {cluster_id}")

        self.state_dir = os.path.join(_ROOT, "4-output", "wflow_state", str(self.cluster_id))
        if not os.path.exists(self.state_dir):
            self.logger.info(
                "No states exist yet for this cluster and forcing combination, creating folder"
            )
            os.makedirs(self.state_dir)

        self.create_state = False
        self.forcing = self.jobs.warmup_forcing

        self.timestepsecs = 86400

        self.path_log = f"log_{self.jobs.name}_{self.cluster_id}_warmup.log"
        self.reinit = False
        self.path_forcing = _FORCING_FILES[self.forcing]
        self.toml = "warmup.toml"

    def get_all_states(self) -> list[datetime.datetime]:
        """
        Returns a list of all available states for the Wflow model forecast.

        This method retrieves all state files from the state director
        and extracts the date from each file name. It will also filter these states
        by the selected forcing type when specified in the original config.

        Returns
        -------
        list[datetime.datetime]
            A list of datetime objects representing the dates of the available states.
        """
        fname = f"{self.forcing}_*.nc" if self.forcing else "*.nc"
        files = glob.glob(os.path.join(self.state_dir, fname))
        datestrings = [fname_state.split("_")[1] for fname_state in files]
        states = [
            datetime.datetime.strptime(date, _DATE_FORMAT_FNAME) for date in datestrings
        ]
        return sorted(states, reverse=True)  # from recent to old

    def get_new_state(self, warmup_days: int = 750) -> datetime.datetime:
        """
        Prepares a Wflow model run to create a new state for a forecast.

        This method prepares a new Wflow model run to be used as a warm state for a forecast.
        It will try to use an exisiting state if
        it can be found within the tolerance of 'warmup_days'.

        Parameters
        ----------
        warmup_days : int, optional
            The tolerance to search for an existing state as starting point.
        """
        self.create_state = True
        self.endtime = self.jobs.tstart
        available_states = self.get_all_states()

        index = bisect.bisect_right(available_states, self.endtime)
        if index:
            state = available_states[index - 1]
            self.logger.info(
                f"Preparing new warmup run, starting at closest existing state {state}"
            )
            self.reinit = False
        else:
            state = self.jobs.tstart - datetime.timedelta(days=warmup_days)
            self.reinit = True
            self.logger.info(
                f"Preparing new warmup run,"
                f"starting with cold state {state} (warmup_days = {warmup_days})"
            )
        self.starttime = state
        return state


class Forecast(Run):
    """
    A class used to represent the Run that forms the Forecast for a specific cluster.

    This class inherits from the Run class, and also contains a State object that represents
    the warm state for the Forecast. It used the toml defined in the global variable _TOML_FORECAST
    as a template for the run settings.
    """

    def __init__(self, jobs: Jobs, cluster_id: int) -> None:
        super().__init__(jobs, cluster_id)
        self.state = State(jobs, cluster_id)
        
        self.forcing = jobs.forcing
        self.starttime = jobs.tstart
        self.endtime = jobs.tend
        self.timestepsecs = jobs.timestepsecs
        self.duration = jobs.duration
        self.path_log = f"log_{jobs.name}_{self.cluster_id}_forecast.log"
        self.reinit = False
        self.path_forcing = _FORCING_FILES[self.forcing]
        self.path_output = "output.nc"
        self.toml = "forecast.toml"

    def prepare(self):
        """
        Prepares the Wflow model run by checking for states,
        writing TOML files and estimating the maximum runtime.

        """
        self.logger.info(f"Preparing run for cluster {self.cluster_id}")
        self.find_recent_state()
        forecast_runtime = self.estimate_max_runtime(self.jobs.clusters)
        self.jobs.runtimes[f"{self.cluster_id}_forecast"] = forecast_runtime
        self.logger.info(
            f"Write Wflow TOML file for forecast, estimated runtime: {time_in_dhms(forecast_runtime)}"
        )
        self.create_toml(template=_TOML_FORECAST)

        if self.state.create_state:
            state_runtime = self.estimate_max_runtime(self.jobs.clusters)
            self.logger.info(
                f"Write Wflow TOML file for state, estimated runtime: {time_in_dhms(state_runtime)}"
            )
            self.state.create_toml(template=_TOML_STATE)
            self.jobs.runtimes[f"{self.cluster_id}_state"] = state_runtime

    def find_recent_state(self, recent_days: int = 14) -> None:
        """
        Finds the most recent state for a Forecast.

        This method searches the most recent state for a Forecast.
        When a matching state if found within the tolerance of 'recent_days',
        it will use that state and change the start date of the forecast to match this state.
        If no recent state is found, it will prepare a new warump run
        that will create the state for the forecast.

        Parameters
        ----------
        recent_days : int, optional
            The number of days to consider for the search of the most recent state.
        """
        self.logger.info(
            f"Searching for most recent state (start date: {self.starttime},"
            f"recent_days: {recent_days})"
        )
        available_states = self.state.get_all_states()
        if len(available_states) > 0:
            self.logger.info(
                f"Found {len(available_states)} existing states for this combination"
                " of cluster and forcing"
            )

            for state in available_states:
                if state == self.starttime:
                    self.logger.info(
                        f"Found state with start date {state}, use as warm state"
                    )
                    break

                elif state < self.starttime <= state + datetime.timedelta(days=recent_days):
                    self.logger.info(
                        f"Found state with start date {state}, use as warm state and new start time"
                    )
                    self.state_input = state
                    self.starttime = state
                    self.duration = int((self.endtime - self.starttime) / self.timestepsecs)
                    break
        else:
            self.logger.info("Found no matching states, need to create a new state")
            state = self.state.get_new_state()
        
        self.set_state_input(state)
 
        
    def set_state_input(self, state: datetime.datetime) -> None:
        """
        Set the input path to the state given a state date

        Parameters
        ----------
        state : datetime.dateime
            The date of the input state
        """
        state_dir = os.path.join(_ROOT, "3-input", "wflow_state", str(self.cluster_id))
        self.state_input = os.path.join(state_dir, f"{self.forcing}_{state.strftime(_DATE_FORMAT_FNAME)}.nc")
        self.logger.debug(f"set state_input to {self.state_input}")
        

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt=r"%Y-%m-%d %H:%M:%S",
    )
    try:
        runs = Jobs(r"C:\git\Moonshot_2_for_Africa\forecasts\2024-10-16_MS2_workshop\forecast_mozambique_freddy.yml")
        # args = argparse.ArgumentParser()
        # args.add_argument("--config", type=str, required=True)
        # args = args.parse_args()
        
        # runs = Jobs(
        #     args.config
        # )
        runs.prepare()
    except Exception as e:
        traceback.print_exc()
        logger.error(e) 
        exit(1)