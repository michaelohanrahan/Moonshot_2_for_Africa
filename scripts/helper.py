import sys 
from pathlib import Path

def syscheck():
    if sys.platform.startswith('win'):
        DRIVE='P:'
    else:
        DRIVE='/p'
    return DRIVE

def create_directories(config):
    dirs = {}
    for key, dir_path in config['data']['directories'].items():
        dir_path = Path(dir_path)
        dir_path.mkdir(exist_ok=True, parents=True)
        dirs[key] = dir_path
    return dirs