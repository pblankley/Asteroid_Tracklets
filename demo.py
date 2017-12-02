# Imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# os.path.join(BASE_DIR, 'data/plist_df.json')
""" Here we will have a demo for the code that gets a subset of data from the ITF
    and runs our clustering process over the subset. It will then give the
    resulting clusters and the realted orbit-fit validated results, along with
    the total number of asteroids we accurately matched, based on the orbit fitting.
"""

## TODO: add __init__.py file to the itf directory and make the whole package more module-like
## TODO: copy process in the "itf mjh trimmed ipynb"
