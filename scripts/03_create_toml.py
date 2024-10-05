import logging
import geopandas as gpd
import datetime
import os
import bisect
import yaml
import glob

from hydromt_wflow import WflowModel
from hydromt.config import configread

logger = logging.getLogger("MS2")

_ROOT = os.path.normpath("p:/moonshot2-casestudy/Wflow/africa/src/3-model")
_TOML_STATE = os.path.normpath("p:/moonshot2-casestudy/Wflow/africa/config/03_warmup_wflow_sbm.toml")
_TOML_FORECAST = os.path.normpath("p:/moonshot2-casestudy/Wflow/africa/config/03_forecast_wflow_sbm.toml")

_DATE_FORMAT_FNAME = r"%Y%m%d"
_DATE_FORMAT_SHORT = r"%Y-%m-%d"
_DATE_FORMAT_LONG = r"%Y-%m-%dT%H:%M:%S"

_FORCING_FILES = {
    "era5": "...",
    "era5_hourly": "...",
    "chirps": "...",
    }

_MAX_RUNTIME = 10 # upper limit for runtime in seconds per timestep per km²
#TODO: dict per cluster type / number of cores

'''
TO DO LIST
- Add consistency checks: to start of script, e.g. do all states have a config assigned to them?
'''

class Config:
    def __init__(self, config_path: str, logger=logger) -> None:
        self.logger = logger
        self.logger.info(f"Preparing forecast! Reading config from {config_path}")
        
        self.read_config(config_path)

    def read_config(self, config_path):
        with open(config_path, 'r') as file:
            config_dict = yaml.safe_load(file)
        
        self.name = config_dict['name']
        self.tstart = datetime.datetime.strptime(config_dict['start_date'], _DATE_FORMAT_SHORT)
        self.timestepsecs = config_dict['timestepsecs']
        self.duration = config_dict['duration_in_timesteps']
        self.tend = self.tstart + datetime.timedelta(seconds= self.timestepsecs * self.duration)
        self.vars = config_dict['variables']
        self.forcing = config_dict['forcing']
        self.warmup_forcing = config_dict["warmup_forcing"]
        
        for _forcing in [self.forcing, self.warmup_forcing]:
            if _forcing not in _FORCING_FILES.keys():
                raise ValueError(f"Invalid forcing type {self.forcing}, choose from {_FORCING_FILES.keys()}")

        self.points = gpd.read_file(os.path.abspath(config_dict['points_of_interest']))

        self.clusters = None

    def locate_cluster(self, clusters: gpd.GeoDataFrame, column_with_ids: str = 'cluster_key') -> list[int]:
        self.logger.info("Finding clusters corresponding to points of interest")
        cluster_ids = gpd.sjoin(self.points, clusters, how='inner', predicate='within')[column_with_ids]
        self.clusters = list(map(lambda x: int(x) if str(x).isdigit() else x, set(cluster_ids)))
        self.logger.info(f"Found one or more clusters: {self.clusters}")


class Run():
    def __init__(self, config: Config, cluster_id: int) -> None:
        self.logger = config.logger
        self.logger.info(f"Initializing forecast for cluster {cluster_id}")
        self.config = config
        self.cluster = cluster_id
        self.toml = None
        self.duration = None
        self.dir_input = os.path.join(_ROOT, "wflow_build", str(self.cluster))
        self.dir_output = os.path.join(_ROOT, "wflow_forecast", self.config.name, str(self.cluster))
        

        # TOML contents
        self.starttime = None
        self.endtime = None
        self.timestepsecs = None
        self.path_log = None
        self.reinit = None
        self.path_forcing = None
        self.path_output = None
        self.state_input = None
        self.state_output = None

    def create_toml(self, template: str, hydromt_config_fn: str = None, limit_logging: bool = True) -> None:
        """
        Writes the configuration for the Wflow model to a TOML file.

        This method sets and writes the configuration for the Wflow model using the Run attributes.
        It uses a dummy instance of WflowModel from hydromt_wflow to update and write the TOML file.
        A hydromt_wflow configuration file can be used to further update specific parts of the TOML file.
        The TOML file for both the state and forecast are written to the output directory of the forecast.

        Parameters
        ----------
        template : str
            The path to the template TOML file to use for the Wflow model configuration.
        hydromt_config_fn : str, optional
            The path to a hydromt_wflow config file. If provided, the configuration options in this file 
            will be read and used to update the Wflow model configuration.
        limit_logging: bool, optional
            Will prevent logging from the hydromt_wflow module to be added to the current log.

        Returns
        -------
        None
        """
        if limit_logging:
            w_logger = logging.getLogger('hydromt_wflow')
            w_logger.setLevel(logging.CRITICAL)
            w_logger.addHandler(logging.NullHandler())
        else:
            w_logger = self.logger
    
        w = WflowModel(config_fn=template, logger=w_logger)
        
        w.set_config("starttime", self.starttime.strftime(_DATE_FORMAT_LONG))
        w.set_config("endtime", self.endtime.strftime(_DATE_FORMAT_LONG))
        w.set_config("timestepsec", self.timestepsecs)
        w.set_config("dir_input", self.dir_input)
        w.set_config("dir_output", self.dir_output)
        w.set_config("path_log", self.path_log)

        w.set_config("model", "reinit", self.reinit)
        w.set_config("input", "path_forcing", self.path_forcing)
        w.set_config("state", "path_input", self.state_input)
        w.set_config("state", "path_output", self.state_output)
        w.set_config("output", "path", self.path_output)

        if hydromt_config_fn:
            self.logger.info(f"Reading hydromt_wflow config: {hydromt_config_fn}")
            opt = configread(config_fn=hydromt_config_fn)
            w.update(write=False, opt=opt)

        self.logger.info(f"Writing Wflow config (warmup) to {self.dir_output}")
        w.write_config(self.toml, self.dir_output)
        w = None # close model

    def estimate_max_runtime(self, clusters: gpd.GeoDataFrame, column_with_ids: str = 'cluster_key') -> None:
        gdf_cluster = clusters[clusters[column_with_ids] == self.cluster]
        gdf_cluster = gdf_cluster.to_crs(epsg=3857)
        area_cluster = gdf_cluster['geometry'].area
        self.max_runtime = area_cluster * self.duration * _MAX_RUNTIME

class State(Run):
    def __init__(self, config: Config, cluster_id: int) -> None:
        super().__init__(config, cluster_id)
        
        self.state_dir = os.path.join(_ROOT, "wflow_state", str(self.cluster))
        if not os.path.exists(self.state_dir):
            self.logger.info("No states exist yet for this cluster and forcing combination, creating folder")
            os.makedirs(self.state_dir)
        
        self.create_state = False
        self.forcing = self.config.warmup_forcing

        self.timestepsecs = 86400
        
        self.path_log = f"log_{self.config.name}_{self.cluster}_warmup.log"
        self.reinit = False
        self.path_forcing = _FORCING_FILES[self.config.forcing]
        self.path_output = f"{self.config.name}_{self.cluster}_output.nc"
        self.toml = f"{self.config.name}_{self.cluster}_warmup.toml"

    def get_all_states(self) -> list[datetime.datetime]:
        fname = f"{self.forcing}_*.nc" if self.forcing else "*.nc"
        files = glob.glob(os.path.join(self.state_dir, fname))
        datestrings = [fname_state.split("_")[1] for fname_state in files]
        states = [datetime.datetime.strptime(date, _DATE_FORMAT_FNAME) for date in datestrings]
        return sorted(states, reverse=True) # from recent to old

    def get_new_state(self, warmup_days: int = 750) -> None:
        available_states = self.get_all_states()
        index = bisect.bisect_right(available_states, self.config.tstart)
        if index:
            state = available_states[index - 1]
            self.logger.info(f"Preparing new warmup run, starting at closest existing state {state}")
            self.reinit = False
        else:
            state = self.config.tstart - datetime.timedelta(days=warmup_days)
            self.reinit = True
            self.logger.info(f"Preparing new warmup run, starting with cold state {state} (warmup_days = {warmup_days})")
        self.state_date = state
        


class Forecast(Run):
    def __init__(self, config: Config, cluster_id: int) -> None:
        super().__init__(config, cluster_id)
        self.state = State(config, cluster_id)
    
        self.starttime = self.config.tstart
        self.endtime = self.config.tend
        self.timestepsecs = self.config.timestepsecs
        self.duration = self.config.duration
        self.path_log = f"log_{self.config.name}_{self.cluster}_forecast.log"
        self.reinit = False
        self.path_forcing = _FORCING_FILES[self.config.forcing]
        self.path_output = f"{self.config.name}_{self.cluster}_output.nc"
        self.toml = f"{self.config.name}_{self.cluster}_warmup.toml"

    def prepare(self):
        self.logger.info(f"Preparing run for cluster {self.cluster}")
        self.find_recent_state()
        self.logger.info("Write Wflow TOML file for forecast")
        self.create_toml(template=_TOML_FORECAST)
        if self.state.create_state:
            self.logger.info("Write Wflow TOML file for state")
            self.state.create_toml(template=_TOML_STATE)

    def find_recent_state(self, recent_days: int = 14) -> None:
        self.logger.info(f"Searching for most recent state (start date: {self.starttime}, recent_days: {recent_days})")
        available_states = self.state.get_all_states(self.config.warmup_forcing)
        self.logger.info(f"Found {len(available_states)} existing states for this combintation of cluster and forcing")
        for state in available_states:
            if state == self.starttime:
                self.logger.info(f"Found state with start date {state}, use as warm state")
                self.state_input = state
                
            elif state < self.starttime <= state + datetime.timedelta(days=recent_days):
                self.logger.info(f"Found state with start date {state}, use as warm state and new start time")
                self.state_input = state
                self.starttime = state
                self.duration = int((self.endtime - self.starttime) / self.timestepsecs)
            else:
                self.logger.info("Found no matching states, need to create a new state")
                self.state.get_new_state()



    

def main(config_path: str, cluster_path: str):
    config = Config(config_path)

    clusters = gpd.read_file(cluster_path)
    config.locate_cluster(clusters)
    for cluster_id in config.clusters:
        forecast = Forecast(config=config, cluster_id=cluster_id)
        forecast.prepare()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO) 

    cnfg_path = "p:/moonshot2-casestudy/Wflow/africa/src/0-setup/forecast_mozambique_freddy.yml"
    clstr_path = "p:/moonshot2-casestudy/Wflow/africa/data/2-interim/clustered_basins.geojson"

    main(cnfg_path, clstr_path)
