#!/usr/bin/env python3
'''
Simple file to plot the output.
'''

import re
import os
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

def main():
    '''Main wrapper'''
    mat_1_filenames = ["./out/steady_state_slab_1.out",
        "./out/steady_state_slab_closure_1.out",
        "./out/steady_state_slab_non_mark_closure_1.out",
        "./out/hold/steady_state_slab_non_mark_closure_1.out"]
    mat_2_filenames = ["./out/steady_state_slab_2.out",
        "./out/steady_state_slab_closure_2.out",
        "./out/steady_state_slab_non_mark_closure_2.out",
        "./out/hold/steady_state_slab_non_mark_closure_2.out"]
    filenames = ["./out/steady_state_slab.out",
        "./out/steady_state_slab_closure.out",
        "./out/steady_state_slab_non_mark_closure.out",
        "./out/hold/steady_state_slab_non_mark_closure.out"]
    for i, filename in enumerate(mat_1_filenames):
        if os.path.isfile(filename):
            x_values, flux = get_data(filename)
            if i == 0:
                legend_title = "Material 1"
                color = 'g'
                linestyle = '-'
            elif i == 1:
                legend_title = "Material 1 LP"
                color = 'r'
                linestyle = '-'
            elif i == 2:
                legend_title = "Material 1 NM, Eta = 0.6"
                color = 'b'
                linestyle = ':'
            else:
                legend_title = "Material 1 NM, Eta = 1.0"
                color = 'k'
                linestyle = '-'
            plt.plot(x_values, flux, color=color, linestyle=linestyle, label=legend_title)
    plt.grid(which='both', axis='both')
    plt.legend()
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.yscale('log')
    #plt.xscale('log')
    plotname = "steady_state_slab_non_mark_closure_overlap_m1.png"
    plt.tight_layout()
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")
    plt.cla()
    plt.clf()
    for i, filename in enumerate(mat_2_filenames):
        if os.path.isfile(filename):
            x_values, flux = get_data(filename)
            if i == 0:
                legend_title = "Material 2"
                color = 'g'
                linestyle = '-'
            elif i == 1:
                legend_title = "Material 2 LP"
                color = 'r'
                linestyle = '-'
            elif i == 2:
                legend_title = "Material 2 NM, Eta = 0.6"
                color = 'b'
                linestyle = ':'
            else:
                legend_title = "Material 2 NM, Eta = 1.0"
                color = 'k'
                linestyle = '-'
            plt.plot(x_values, flux, color=color, linestyle=linestyle, label=legend_title)
    plt.grid(which='both', axis='both')
    plt.legend()
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.yscale('log')
    #plt.xscale('log')
    plotname = "steady_state_slab_non_mark_closure_overlap_m2.png"
    plt.tight_layout()
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")
    plt.cla()
    plt.clf()
    for i, filename in enumerate(filenames):
        if os.path.isfile(filename):
            x_values, flux = get_data(filename)
            if i == 0:
                legend_title = "Benchmark"
                color = 'g'
                linestyle = '-'
            elif i == 1:
                legend_title = "LP"
                color = 'r'
                linestyle = '-'
            elif i == 2:
                legend_title = "NM, Eta = 0.6"
                color = 'b'
                linestyle = ':'
            else:
                legend_title = "NM, Eta = 1.0"
                color = 'k'
                linestyle = '-'
            plt.plot(x_values, flux, color=color, linestyle=linestyle, label=legend_title)
    plt.grid(which='both', axis='both')
    plt.legend()
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.yscale('log')
    #plt.xscale('log')
    plotname = "steady_state_slab_non_mark_closure_overlap.png"
    plt.tight_layout()
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")
    plt.cla()
    plt.clf()

if __name__ == '__main__':
    main()
