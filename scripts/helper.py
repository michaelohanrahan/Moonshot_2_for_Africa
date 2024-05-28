import sys 

def syscheck():
    if sys.platform.startswith('win'):
        DRIVE='P:'
    else:
        DRIVE='/p'
    return DRIVE