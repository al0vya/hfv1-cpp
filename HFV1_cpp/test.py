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
        "python test.py plot <MODE>\n" +
        "MODE : [debug,release]"
    )
    
    sys.exit(help_message)

def plot_soln(
        x1,
        x2,
        y,
        quantity,
        length,
        ylabel,
        test_num=0,
        test_name="ad-hoc",
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
        mode
    ):
        if   mode == "debug":
            path = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "solution_data.csv")
        elif mode == "release":
            path = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "solution_data.csv")
        else:
            EXIT_HELP()
        
        print("Searching for solution data in", path)
        
        dataframe = pd.read_csv(path)
        
        self.x1  = dataframe["x1"].values
        self.x2  = dataframe["x2"].values
        self.q   = dataframe["q"].values
        self.z   = dataframe["z"].values
        self.eta = dataframe["eta"].values
        
        self.length = self.x1.size
        
    def plot_soln(
        self
    ):
        plot_soln(self.x1, self.x2, self.q,   "q",   self.length, "$q \, (m^2s^{-1})$")
        plot_soln(self.x1, self.x2, self.eta, "eta", self.length, "$\eta \, (m)$")
        plot_soln(self.x1, self.x2, self.z,   "z",   self.length, "$z \, (m)$")
    
def plot():
    if len(sys.argv) > 2:
        dummy, action, mode = sys.argv
        
        Solution(mode).plot_soln()
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
    
    if action == "plot":
        plot()
    else:
        EXIT_HELP()
else:
    EXIT_HELP()
   


