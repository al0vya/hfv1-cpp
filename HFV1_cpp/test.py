import os
import sys
import imageio
import subprocess
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab  as pylab

from mpl_toolkits.mplot3d import Axes3D

####################################
# NATURAL SORTING, READ UP ON THIS #
####################################

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def get_filenames_natural_order():
    path = os.path.dirname(__file__)
    
    filenames_natural_order = os.listdir(path)
    
    filenames_natural_order.sort(key=natural_keys)
    
    return filenames_natural_order

####################################
####################################
####################################

def EXIT_HELP():
    help_message = (
        "Use as:\n\n" +
        "    python test.py test <MODE> <NUM_CELLS> <REFINEMENT_LEVEL> <EPSILON> <SAVE_INT>\n\n" +
        "        MODE     : [debug,release]\n"
        "        SAVE_INT : interval in seconds that the solution data are saved\n\n"
        "    python test.py run <MODE> <TEST_CASE> <NUM_CELLS> <REFINEMENT_LEVEL> <EPSILON> <SAVE_INT>\n\n" +
        "        MODE      : [debug,release]\n" +
        "        TEST_CASE : [1,2,3,4,5,6]\n" +
        "        SAVE_INT  : interval in seconds that the solution data are saved\n\n"
        "    Available test cases:\n" +
        "        1. Wet dam break\n" +
        "        2. Dry dam break\n" +
        "        3. Dry dam break with friction\n" +
        "        4. Wet c property\n" +
        "        5. Wet dry c property\n" +
        "        6. Building overtopping"
    )
    
    sys.exit(help_message)

def clear_jpg_files():
    path = os.path.dirname(__file__)
    
    [ os.remove(filename) for filename in os.listdir(path) if filename.endswith(".jpg") ]

test_names = [
    "wet-dam-break",
    "dry-dam-break",
    "dry-dam-break-fric",
    "wet-c-prop",
    "wet-dry-c-prop",
    "building-overtopping",
]

sim_times = [
    2.5,
    1.3,
    1.3,
    0.5,
    0.5,
    10
]

def plot_soln(
        x1,
        x2,
        y,
        quantity,
        interval,
        length,
        ylabel,
        test_num,
        test_name
    ):
        filename = test_name + "-" + str(interval) + "-" + quantity + ".jpg"
        
        xlim = [ np.amin(x1), np.amax(x2) ]
        
        ylim = [np.amin(y) * 0.98, np.amax(y) * 1.02]
        
        for i in range(length):
            plt.plot( [ x1[i], x2[i] ], [ y[i], y[i] ], color='b', linewidth="1.5" )    
            
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel("$x \, (m)$")
        plt.ylabel(ylabel)
        plt.savefig(filename, bbox_inches="tight")
        plt.clf()

class Solution:
    def __init__(
        self,
        mode,
        interval,
        test_num=0,
        test_name="ad-hoc"
    ):
        print("Initialising solution...")
        
        if mode != "debug" and mode != "release":
            EXIT_HELP()
        
        self.interval  = interval
        self.test_num  = test_num
        self.test_name = test_name
        
        filename = "solution_data-" + str(self.interval) + ".csv"
        
        dataframe = pd.read_csv(filename)
        
        self.x1  = dataframe["x1"].values
        self.x2  = dataframe["x2"].values
        self.q   = dataframe["q"].values
        self.z   = dataframe["z"].values
        self.eta = dataframe["eta"].values
        
        self.length = self.x1.size
        
    def plot_soln(
        self
    ):
        print("Plotting solution " + str(self.interval) + "...")
        
        plot_soln(self.x1, self.x2, self.q,   "q",   self.interval, self.length, "$q \, (m^2s^{-1})$", self.test_num, self.test_name)
        plot_soln(self.x1, self.x2, self.eta, "eta", self.interval, self.length, "$\eta \, (m)$", self.test_num, self.test_name)
        plot_soln(self.x1, self.x2, self.z,   "z",   self.interval, self.length, "$z \, (m)$", self.test_num, self.test_name)

def animate():
    images = []
    
    filenames_natural_order = get_filenames_natural_order()
    
    vars = ["eta", "q"]
    
    for test_name in test_names:
        for var in vars:
            suffix = var + ".jpg"
            for filename in filenames_natural_order:
                if filename.startswith(test_name) and filename.endswith(suffix):
                    image = imageio.imread(filename)
                    
                    images.append(image)
            
            if images:
                imageio.mimsave(test_name + "-" + var + ".gif", images)
                images = []

def run():
    print("Attempting to run...")
    
    if len(sys.argv) > 7:
        dummy, action, mode, test, num_cells, max_ref_lvl, epsilon, saveint = sys.argv
        
        if   mode == "debug":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "HFV1_cpp.exe")
        elif mode == "release":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "HFV1_cpp.exe")
        else:
            EXIT_HELP()
        
        subprocess.run( [solver_file, str(test), num_cells, max_ref_lvl, epsilon, saveint] )
        
        intervals = int( sim_times[int(test) - 1] / float(saveint) )
        
        [ Solution(mode, interval, test, test_names[int(test) - 1]).plot_soln() for interval in range(intervals) ]
    else:
        EXIT_HELP()
        
def run_tests():
    print("Attempting to run all tests...")
    
    if len(sys.argv) > 6:
        dummy, action, mode, num_cells, max_ref_lvl, epsilon, saveint = sys.argv
        
        if   mode == "debug":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "HFV1_cpp.exe")
        elif mode == "release":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "HFV1_cpp.exe")
        else:
            EXIT_HELP()
        
        tests = [1, 2, 3, 4, 5, 6]
        
        for i, test in enumerate(tests):
            subprocess.run( [solver_file, str(test), num_cells, max_ref_lvl, epsilon] )
            Solution(mode, test, test_names[i]).plot_soln()
    else:
        EXIT_HELP()

params = {
    "legend.fontsize" : "xx-large",
    "axes.labelsize"  : "xx-large",
    "axes.titlesize"  : "xx-large",
    "xtick.labelsize" : "xx-large",
    "ytick.labelsize" : "xx-large"
}

pylab.rcParams.update(params)

if len(sys.argv) > 1:
    clear_jpg_files()
    
    action = sys.argv[1]
    
    if action == "run":
        run()
    elif action == "test":
        run_tests()
    else:
        EXIT_HELP()
        
    animate()
else:
    EXIT_HELP()