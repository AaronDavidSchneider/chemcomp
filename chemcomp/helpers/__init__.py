import os
from ast import literal_eval

CONFIG_PATH = os.path.join(os.getcwd(), "config")  # set default path of config directory
JOB_PATH = os.path.join(os.getcwd(), "jobs")  # set default path of config directory
OUTPUT_PATH = os.path.join(os.getcwd(), "output")  # set default path of output directory
ZIP_OUTPUT_PATH = os.path.join(os.getcwd(), "zip_output")  # set default path of output directory

def eval_kwargs(inp):
    try:
        if isinstance(inp, str):
            return literal_eval(inp)
        else:
            return inp
    except ValueError:
        return inp


from .import_config import import_config
from .main_helpers import run_setup, run_pipeline
from .solvediffonedee import solvediffonedee_components, getfluxonedee
from .analysis_helper import GrowthPlot, molecules, elements

