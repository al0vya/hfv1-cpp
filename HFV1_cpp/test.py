import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab  as pylab

from mpl_toolkits.mplot3d import Axes3D

params = {
    "legend.fontsize" : "xx-large",
    "axes.labelsize"  : "xx-large",
    "axes.titlesize"  : "xx-large",
    "xtick.labelsize" : "xx-large",
    "ytick.labelsize" : "xx-large"
}

pylab.rcParams.update(params)

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
        self
    ):
        self.x1  = pd.read_csv("solution_data.csv")["x1"].values
        self.x2  = pd.read_csv("solution_data.csv")["x2"].values
        self.q   = pd.read_csv("solution_data.csv")["q"].values
        self.z   = pd.read_csv("solution_data.csv")["z"].values
        self.eta = pd.read_csv("solution_data.csv")["eta"].values
        
        self.length = self.x1.size
        
    def plot_soln(
        self
    ):
        plot_soln(self.x1, self.x2, self.q,   "q",   self.length, "$q \, (m^2s^{-1})$")
        plot_soln(self.x1, self.x2, self.eta, "eta", self.length, "$\eta \, (m)$")
        plot_soln(self.x1, self.x2, self.z,   "z",   self.length, "$z \, (m)$")
        
Solution().plot_soln()