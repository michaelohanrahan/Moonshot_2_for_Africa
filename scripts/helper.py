import sys 
import logging
import os
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

def setup_logging(dir, logname):
    # Set up logging
    log_file_path = Path(Path(os.getcwd()).as_posix(), dir, logname)
    # print('LOGGING TO:', log_file_path)
    os.makedirs(Path(os.path.dirname(log_file_path)).as_posix(), exist_ok=True)

    # Create a custom logger
    logger = logging.getLogger(logname)
    logger.setLevel(logging.INFO)  # Set the root logger level

    # Create handlers
    c_handler = logging.StreamHandler()  # Console handler
    f_handler = logging.FileHandler(log_file_path)  # File handler
    c_handler.setLevel(logging.INFO)  # Set level for console handler
    f_handler.setLevel(logging.INFO)  # Set level for file handler

    # Create formatters and add it to handlers
    format_str = '%(name)s - %(asctime)s - %(levelname)s - %(message)s'  # Modify as needed
    c_format = logging.Formatter(format_str)
    f_format = logging.Formatter(format_str)
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)
    assert os.path.exists(log_file_path)
    return logger

if __name__ == "__main__":
    os.chdir('..')
    # Test the logging setup
    l = setup_logging('data/0-log', 'test_logging.log')
    l.info("This is a test log message.")
    p = syscheck()
    l.info(f"System check: {p}")