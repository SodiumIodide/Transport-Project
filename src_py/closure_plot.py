#!/usr/bin/env python3
'''
Simple file to plot the output.
'''

import re
import numpy as np
import matplotlib.pyplot as plt

def main():
    '''Main wrapper'''
    plotname = "steady_state_slab_closure.png"
    x_values = np.array([])
    flux = np.array([])
    pattern = re.compile(r'\s+(\S+)\s+(\S+)\s+')
    with open('./out/steady_state_slab_closure.out', mode='r') as datafile:
        for line in datafile:
            if re.match(pattern, line):
                x_values = np.append(x_values,
                                     float(re.match(pattern, line).group(1)))
                flux = np.append(flux,
                                 float(re.match(pattern, line).group(2)))
    plt.plot(x_values, flux)
    plt.grid(which='major', axis='both')
    plt.title("Steady State Slab")
    plt.xlabel("x (cm)")
    plt.ylabel("Flux (1/cm^2-s-MeV)")
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")

if __name__ == '__main__':
    main()
