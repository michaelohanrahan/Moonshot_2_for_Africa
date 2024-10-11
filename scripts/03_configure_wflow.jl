using Wflow

config_template = "c:/Users/hartgrin/OneDrive - Stichting Deltares/Projecten/Moonshot/Moonshot_2_for_Africa/config/03_wflow_sbm.toml"

config = Wflow.Config(config_template)

config.starttime = ""
config.endtime = ""
config.state.path_input = ""

config.path_forcing = ""