#!/usr/bin/env python3
'''
Simple file to plot the output.
'''

import re
import numpy as np
import matplotlib.pyplot as plt

def get_data(filename):
    '''Read in the files and return data for plotting'''
    x_values = np.array([])
    flux = np.array([])
    pattern = re.compile(r'\s*(\S+)\s*(\S+)\s*')
    with open(filename, mode='r') as datafile:
        for line in datafile:
            if re.match(pattern, line):
                x_values = np.append(x_values,
                                     float(re.match(pattern, line).group(1)))
                flux = np.append(flux,
                                 float(re.match(pattern, line).group(2)))
    return (x_values, flux)

def realization_plot(realname):
    '''Plot of a geometry realization'''
    plotname = realname.replace("out", "plots").replace(".plots", ".png")
    x_values, flux = get_data(realname)
    plt.plot(x_values, flux)
    plt.grid(which='major', axis='both')
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.savefig(plotname)
    plt.clf()
    plt.cla()

def hist_plot(histname):
    '''Plot of a reflection or transmission PDF histogram'''
    if "refl" in histname:
        label = "Reflection"
    else:
        label = "Transmission"
    plotname = histname.replace("out", "plots").replace(".plots", ".png")
    prob, num = get_data(histname)
    # Trim zero values here
    prob = np.array([p for i, p in enumerate(prob) if num[i] != 0])
    num = np.array([n for n in num if n != 0])
    plt.plot(prob, num)
    plt.xlim(xmin=0.0)
    plt.xlabel(f"Percentage of {label}")
    plt.ylabel(f"Probability of {label}")
    plt.savefig(plotname)
    plt.clf()
    plt.cla()

def main():
    '''Main wrapper'''
    plotname = "steady_state_slab.png"
    filenames = ["./out/steady_state_slab.out",
                 "./out/steady_state_slab_1.out",
                 "./out/steady_state_slab_2.out"]
    for i, filename in enumerate(filenames):
        x_values, flux = get_data(filename)
        if i == 0:
            legend_title = "Total Material"
        else:
            legend_title = "Material " + str(i)
        plt.plot(x_values, flux, label=legend_title)
    plt.grid(which='major', axis='both')
    plt.legend()
    #plt.title("Steady State Slab")
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")
    plt.clf()
    plt.cla()
    realnames = ["./out/realization_1.out",
                 "./out/realization_2.out",
                 "./out/realization_3.out",
                 "./out/last_realization.out"]
    for realname in realnames:
        realization_plot(realname)
    histnames = ["./out/trans_histogram.out",
                 "./out/refl_histogram.out"]
    for histname in histnames:
        hist_plot(histname)

if __name__ == '__main__':
    main()
