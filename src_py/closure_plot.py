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

def main():
    '''Main wrapper'''
    plotname = "steady_state_slab_closure.png"
    filenames = ["./out/steady_state_slab_closure.out",
                 "./out/steady_state_slab_closure_1.out",
                 "./out/steady_state_slab_closure_2.out"]
    for i, filename in enumerate(filenames):
        x_values, flux = get_data(filename)
        if i == 0:
            legend_title = "Total Material"
        else:
            legend_title = "Material " + str(i)
        plt.plot(x_values, flux, label=legend_title)
    plt.grid(which='major', axis='both')
    plt.legend()
    plt.title("Steady State Slab")
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")

if __name__ == '__main__':
    main()
