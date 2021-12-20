import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab  as pylab

from mpl_toolkits.mplot3d import Axes3D

def EXIT_HELP():
    help_message = (
        "Use as:\n\n" +
        "    python test.py test <MODE> <NUM_CELLS> <REFINEMENT_LEVEL> <EPSILON>\n\n" +
        "        MODE : [debug,release]\n\n"
        "    python test.py run <MODE> <TEST_CASE> <NUM_CELLS> <REFINEMENT_LEVEL> <EPSILON>\n\n" +
        "        MODE      : [debug,release]\n" +
        "        TEST_CASE : [1,2,3,4,5,6]\n\n" +
        "    Available test cases:\n" +
        "        1. Wet dam break\n" +
        "        2. Dry dam break\n" +
        "        3. Dry dam break with friction\n" +
        "        4. Wet c property\n" +
        "        5. Wet dry c property\n" +
        "        6. Building overtopping"
    )
    
    sys.exit(help_message)

test_names = [
    "wet-dam-break",
    "dry-dam-break",
    "dry-dam-break-fric",
    "wet-c-prop",
    "wet-dry-c-prop",
    "building-overtopping",
]

def plot_soln(
        x1,
        x2,
        y,
        quantity,
        length,
        ylabel,
        test_num,
        test_name
    ):
        filename = str(test_num) + "-" + test_name + "-" + quantity
        
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
        test_num=0,
        test_name="ad-hoc"
    ):
        print("Initialising solution...")
        
        if mode != "debug" and mode != "release":
            EXIT_HELP()
        
        self.test_num  = test_num
        self.test_name = test_name
        
        dataframe = pd.read_csv("solution_data.csv")
        
        self.x1  = dataframe["x1"].values
        self.x2  = dataframe["x2"].values
        self.q   = dataframe["q"].values
        self.z   = dataframe["z"].values
        self.eta = dataframe["eta"].values
        
        self.length = self.x1.size
        
    def plot_soln(
        self
    ):
        print("Plotting solution...")
        
        plot_soln(self.x1, self.x2, self.q,   "q",   self.length, "$q \, (m^2s^{-1})$", self.test_num, self.test_name)
        plot_soln(self.x1, self.x2, self.eta, "eta", self.length, "$\eta \, (m)$", self.test_num, self.test_name)
        plot_soln(self.x1, self.x2, self.z,   "z",   self.length, "$z \, (m)$", self.test_num, self.test_name)
    
def run():
    print("Attempting to run...")
    
    if len(sys.argv) > 6:
        dummy, action, mode, test, num_cells, max_ref_lvl, epsilon = sys.argv
        
        if   mode == "debug":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "FV1_cpp.exe")
        elif mode == "release":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "FV1_cpp.exe")
        else:
            EXIT_HELP()
        
        subprocess.run( [solver_file, str(test), num_cells, max_ref_lvl, epsilon] )
        Solution(mode, test, test_names[int(test) - 1]).plot_soln()
    else:
        EXIT_HELP()
        
def run_tests():
    print("Attempting to run all tests...")
    
    if len(sys.argv) > 5:
        dummy, action, mode, num_cells, max_ref_lvl, epsilon = sys.argv
        
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
    action = sys.argv[1]
    
    if action == "run":
        run()
    elif action == "test":
        run_tests()
    else:
        EXIT_HELP()
else:
    EXIT_HELP()