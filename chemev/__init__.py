"""
chemev
=======

A modular Chemical Evolution code that hopes to get contributions from
the entire astrophysical community.

For more information, either build the latest documentation included
in our git repository, or view the online version here:
http://chemev.github.io/chemev/

"""

# This file is based on a copy from pynbody. 
# pynbody was written mostly by Andrew Pontzen and Rok Roskar 
# with contributions from Greg Stinson

from . import backcompat

# Import basic dependencies
import ConfigParser
import os
import imp
import numpy
import warnings
import sys
import logging
import multiprocessing


# Create config dictionaries which will be required by subpackages
# We use the OrderedDict, which is default in 2.7, but provided here for 2.6/2.5 by
# the backcompat module. This keeps things in the order they were parsed (important
# for units module, for instance).
config_parser = ConfigParser.ConfigParser(dict_type=backcompat.OrderedDict)
config = {}


# Process configuration options
config_parser.optionxform = str
config_parser.read(
    os.path.join(os.path.dirname(__file__), "default_config.ini"))
config_parser.read(os.path.join(os.path.dirname(__file__), "config.ini"))
config_parser.read(os.path.expanduser("~/.chemevrc"))
config_parser.read("config.ini")


config = {'verbose': config_parser.getboolean('general', 'verbose')}


config['threading'] = config_parser.get('general', 'threading')
config['number_of_threads'] = int(
    config_parser.get('general', 'number_of_threads'))

if config['number_of_threads']<0:
    config['number_of_threads']=multiprocessing.cpu_count()

# Create the logger for chemev
logger = logging.getLogger('chemev')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(name)s : %(message)s')
for existing_handler in list(logger.handlers):
    logger.removeHandler(existing_handler)

ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)

if config['verbose']:
    ch.setLevel(logging.INFO)
    logger.info("Verbose mode is on")
else:
    ch.setLevel(logging.WARNING)

    warning = """
Welcome to chemev v0.30. Note this new version by default is much quieter than old versions.
To get back the verbose output, edit your config.ini or .chemevrc file and insert the following
section

[general]
verbose: True

The information is now parsed through python's standard logging module; using logging.getLogger('chemev')
you can customize the behaviour. See here https://docs.python.org/2/howto/logging-cookbook.html#logging-cookbook."""

    if not os.path.exists(os.path.expanduser("~/.chemev_v03_touched")):
        print warning
        with open(os.path.expanduser("~/.chemev_v03_touched"), "w") as f:
            print>>f, "This file tells chemev not to reprint the welcome-to-v-0.3 warning"



# Import subpackages
from . import engine, snii, starlifetime, zones, agb, enrich, snia, star, imf, data

try:
    from . import plot
except:
    warnings.warn(
        "Unable to import plotting package (missing matplotlib or running from a text-only terminal? Plotting is disabled.", RuntimeWarning)

# The following code resolves inter-dependencies when reloading
imp.reload(engine)
imp.reload(imf)
imp.reload(snii)
# imp.reload(family) # reloading this causes problems for active snapshots
imp.reload(snia)
imp.reload(agb)
imp.reload(star)
imp.reload(starlifetime)
imp.reload(enrich)
imp.reload(zones)
imp.reload(data)

try:
    imp.reload(plot)
except:
    pass


#def load(filename, *args, **kwargs):
#    """Loads a file using the appropriate class, returning a SimSnap
#    instance."""

#    for c in config['snap-class-priority']:
#        if c._can_load(filename):
#            logger.info("Loading using backend %s" % str(c))
#            return c(filename, *args, **kwargs)

#    raise IOError(
#        "File %r: format not understood or does not exist" % filename)


#from snapshot import _new as new

#derived_array = snapshot.SimSnap.derived_quantity

#__all__ = ['load', 'new', 'derived_array']
