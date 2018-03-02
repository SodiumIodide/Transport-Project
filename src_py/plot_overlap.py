#!/usr/bin/env python3
'''
Python script to loop through benchmark data and produce presentation-ready
graphics involving a comparison of solution methods.
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

def plot_overlap(first_x, first_y, first_label, second_x, second_y, second_label, destination):
    '''Produce a figure of overlapping data'''
    length = 3  # Guaranteed length of array: total, mat_1, mat_2
    colors = ['b', 'c', 'r']
    for index in range(length):
        if index == 0:
            legend_entry = "Total Material"
        else:
            legend_entry = f"Material {index}"
        plt.plot(first_x[index], first_y[index], label=f"{first_label} {legend_entry}", color=colors[index])
        plt.plot(second_x[index], second_y[index], label=f"{second_label} {legend_entry}", linestyle=':', color=colors[index])
    plt.grid(which='major', axis='both')
    plt.legend()
    plt.xlabel('x (cm)')
    plt.ylabel('Flux (1/cm^2-s-MeV)')
    plt.xlim(xmin=0.0)
    plt.ylim(ymin=0.0)
    plt.savefig(f"{destination}.png")
    plt.clf()
    plt.cla()

def main():
    '''Main wrapper'''
    directory = "./benchmarks"
    print("Creating files in:")
    for item in os.listdir(directory):
        datafolder = f"{directory}/{item}/data"
        if os.path.isdir(f"{directory}/{item}") and os.path.exists(datafolder):
            exact_x = []
            exact_y = []
            model_x = []
            model_y = []
            alpha_x = []
            alpha_y = []
            print(f"{directory}/{item}")
            for datafile in os.listdir(datafolder):
                filename = f"{datafolder}/{datafile}"
                if "alpha" in datafile:
                    # Order is: total, mat_1, mat_2
                    temp_x, temp_y = get_data(filename)
                    alpha_x.append(temp_x)
                    alpha_y.append(temp_y)
                elif "closure" in datafile:
                    # Order is: total, mat_1, mat_2
                    temp_x, temp_y = get_data(filename)
                    model_x.append(temp_x)
                    model_y.append(temp_y)
                elif "slab" in datafile:
                    # Order is: total, mat_1, 
                    temp_x, temp_y = get_data(filename)
                    exact_x.append(temp_x)
                    exact_y.append(temp_y)
            destination_1 = f"{directory}/{item}/closure_overlap"
            destination_2 = f"{directory}/{item}/alpha_overlap"
            if exact_x and model_x:
                plot_overlap(exact_x, exact_y, "Exact", model_x, model_y, "Closure", destination_1)
            if exact_x and alpha_x:
                plot_overlap(exact_x, exact_y, "Exact", alpha_x, alpha_y, "Alpha", destination_2)

if __name__ == '__main__':
    main()
