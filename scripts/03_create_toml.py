import logging
import geopandas as gpd
import datetime
import os
import bisect
import yaml
import glob

from hydromt_wflow import WflowModel

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


class Forecast:
    def __init__(self, config: Config, cluster_id: int) -> None:
        self.logger = config.logger
        self.logger.info(f"Initializing forecast for cluster {cluster_id}")
        self.config = config
        self.cluster = cluster_id
        self.model_dir = os.path.join(_ROOT, "wflow_build", str(self.cluster))
        
        self.state_dir = os.path.join(_ROOT, "wflow_state", str(self.cluster), "state")
        self.state_config = os.path.join(_ROOT, "wflow_state", str(self.cluster), "config")
        self.cold_state = False
        self.create_state = False

        self.forecast_dir = os.path.join(_ROOT, "wflow_forecast", config.name, str(self.cluster))
  
    def get_recent_state(self, recent_days: int = 14) -> None:
        self.logger.info(f"Searching for most recent state (start date: {self.config.tstart})")
        if not os.path.exists(self.state_dir):
            self.logger.info("No states exist for this cluster and forcing combination, creating folder")
            os.makedirs(self.state_dir)
            os.makedirs(self.state_config)
            available_states = []
        else:
            available_states = self.get_all_states(self.config.warmup_forcing)
            self.logger.info(f"Searching through existing states (found {len(available_states)})")
        for state in available_states:
            if state == self.config.tstart:
                self.state_date = state
                self.logger.info(f"Found state with start date {state}, use as warm state")
                return
            if state < self.config.tstart <= state + datetime.timedelta(days=recent_days):
                self.state_date = state
                self.config.tstart = state
                #TODO: update forecast duration as it gets extended

                self.logger.info(f"Found state with start date {state}, use as warm state and new start time")
                return
        self.create_state = True # no recent state? create new state
        self.get_new_state()

        
    def get_all_states(self, filter_forcing: str = None) -> list[datetime.datetime]:
        fname = f"{filter_forcing}_*.nc" if filter_forcing else "*.nc"
        files = glob.glob(os.path.join(self.state_dir, fname))
        datestrings = [fname_state.split("_")[1] for fname_state in files]
        states = [datetime.datetime.strptime(date, _DATE_FORMAT_FNAME) for date in datestrings]
        return sorted(states, reverse=True) # from recent to old

    def get_new_state(self, warmup_days: int = 750) -> None:
        available_states = self.get_all_states(self.config.warmup_forcing)
        index = bisect.bisect_right(available_states, self.config.tstart)
        if index:
            state = available_states[index - 1]
            self.logger.info(f"Found no matching states, preparing warmup run, starting at existing state {state}")
        else:
            state = self.config.tstart - datetime.timedelta(days=warmup_days)
            self.cold_state = True
            self.logger.info(f"Found no matching states, preparing warmup run, starting with cold state {state} (warmup_days = {warmup_days})")
        self.state_date = state


    def create_toml_state(self):
        # dummy Wflow model to access config methods and limit logging from hydromt_wflow
        w_logger = logging.getLogger('hydromt_wflow')
        w_logger.setLevel(logging.CRITICAL)
        w_logger.addHandler(logging.NullHandler())
        w = WflowModel(logger=w_logger)
        
        w.read_config(_TOML_STATE)

        w.set_config("starttime", self.state_date.strftime(_DATE_FORMAT_LONG))
        w.set_config("endtime", self.config.tstart.strftime(_DATE_FORMAT_LONG))
        w.set_config("timestepsec", self.config.timestepsecs)
        w.set_config("dir_input", self.model_dir)
        w.set_config("dir_output", self.forecast_dir)
        w.set_config("path_log", f"log_{self.config.name}_{self.cluster}_warmup.log")

        w.set_config("model", "reinit", self.cold_state)
        w.set_config("input", "path_forcing", self.config.forcing)
        w.set_config("state", "path_input", os.path.join(f"{self.state_dir}", f"{self.state_date.strftime(_DATE_FORMAT_FNAME)}_{self.config.warmup_forcing}.nc"))
        w.set_config("state", "path_output", os.path.join(f"{self.state_dir}", f"{self.config.tstart.strftime(_DATE_FORMAT_FNAME)}_{self.config.warmup_forcing}.nc"))
   
        self.logger.info(f"Writing Wflow config (warmup) to {self.state_config}")
        w.write_config(f"{self.config.tstart.strftime(_DATE_FORMAT_FNAME)}_{self.config.warmup_forcing}.toml", self.state_config)
        w = None

    def create_toml_forecast(self):
        # dummy Wflow model to access config methods and limit logging from hydromt_wflow
        w_logger = logging.getLogger('hydromt_wflow')
        w_logger.setLevel(logging.CRITICAL)
        w_logger.addHandler(logging.NullHandler())
        w = WflowModel(logger=w_logger)
        
        w.read_config(_TOML_FORECAST)

        w.set_config("starttime", self.config.tstart.strftime(_DATE_FORMAT_LONG))
        w.set_config("endtime", self.config.tend.strftime(_DATE_FORMAT_LONG))
        w.set_config("timestepsec", self.config.timestepsecs)
        w.set_config("dir_input", self.model_dir)
        w.set_config("dir_output", self.forecast_dir)
        w.set_config("path_log", f"log_{self.config.name}_{self.cluster}_forecast.log")

        w.set_config("model", "reinit", self.cold_state)
        w.set_config("input", "path_forcing", self.config.forcing)
        w.set_config("state", "path_input", os.path.join(f"{self.state_dir}", f"{self.config.tstart.strftime(_DATE_FORMAT_FNAME)}_{self.config.warmup_forcing}.nc"))
        w.set_config("output", "path", f"{self.config.name}_{self.cluster}_output.nc")
        
        self.logger.info(f"Writing Wflow config (forecast) to {self.forecast_dir}")
        w.write_config(f"{self.config.name}_{self.cluster}.toml", self.forecast_dir)
        w = None

    def estimate_runtime(self, clusters: gpd.GeoDataFrame, column_with_ids: str = 'cluster_key') -> None:
        gdf_cluster = clusters[clusters[column_with_ids] == self.cluster]
        gdf_cluster = gdf_cluster.to_crs(epsg=3857)
        area_cluster = gdf_cluster['geometry'].area
        if self.create_state:
            self.state_runtime = area_cluster * (self.config.tstart - self.state_date) / datetime.timedelta(seconds=86400)
        self.forecast_runtime = area_cluster * self.config.duration

def main(config_path: str, cluster_path: str):
    config = Config(config_path)

    clusters = gpd.read_file(cluster_path)
    config.locate_cluster(clusters)
    for cluster_id in config.clusters:
        forecast = Forecast(config=config, cluster_id=cluster_id)
        forecast.get_recent_state()
        if forecast.create_state:
            forecast.create_toml_state()
        forecast.create_toml_forecast()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO) 

    cnfg_path = "p:/moonshot2-casestudy/Wflow/africa/src/0-setup/forecast_mozambique_freddy.yml"
    clstr_path = "p:/moonshot2-casestudy/Wflow/africa/data/2-interim/clustered_basins.geojson"

    main(cnfg_path, clstr_path)
