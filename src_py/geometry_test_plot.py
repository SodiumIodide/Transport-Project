#!/usr/bin/env python3
'''
Simple file to plot the output of a geometry test.
'''

import re
import numpy as np
import matplotlib.pyplot as plt

def main():
    '''Main wrapper'''
    plotname = "geometry_test.png"
    x_values = np.array([])
    mat_1_x = np.array([])
    mat_2_x = np.array([])
    material = np.array([])
    mat_1 = np.array([])
    mat_2 = np.array([])
    pattern = re.compile(r'\s+(\S+)\s+(\S+)\s+')
    switch_flag_1 = False
    switch_flag_2 = False
    with open('./out/geometry_test.out', mode='r') as datafile:
        for line in datafile:
            if re.match(pattern, line):
                x_values = np.append(x_values,
                                     float(re.match(pattern, line).group(1)))
                material = np.append(material,
                                     int(re.match(pattern, line).group(2)))
                if material[-1] == 1:
                    if switch_flag_2:
                        mat_2_x = np.append(mat_2_x, x_values[-1])
                        mat_2 = np.append(mat_2, 2)
                        switch_flag_2 = False
                    switch_flag_1 = True
                    mat_1_x = np.append(mat_1_x, x_values[-1])
                    mat_1 = np.append(mat_1, 1)
                if material[-1] == 2:
                    if switch_flag_1:
                        mat_1_x = np.append(mat_1_x, x_values[-1])
                        mat_1 = np.append(mat_1, 1)
                        switch_flag_1 = False
                    switch_flag_2 = True
                    mat_2_x = np.append(mat_2_x, x_values[-1])
                    mat_2 = np.append(mat_2, 2)
    if material[-1] == 1:
        mat_1_x = mat_1_x[0:-1]
    else:
        mat_2_x = mat_2_x[0:-1]
    counter = 0
    while counter < len(mat_1_x):
        plt.plot(mat_1_x[counter:counter + 2], mat_1[counter:counter + 2], color='b')
        plt.axvline(x=mat_1_x[counter + 1], linestyle=':', color='k')
        counter += 2
    counter = 0
    while counter < len(mat_2_x):
        plt.plot(mat_2_x[counter:counter + 2], mat_2[counter:counter + 2], color='r')
        plt.axvline(x=mat_2_x[counter + 1], linestyle=':', color='k')
        counter += 2
    plt.title("Sample Markovian Geometry")
    plt.xlabel("x")
    plt.ylabel("Material Number")
    plt.grid(which='major', axis='x')
    plt.ylim(0, 5)
    plt.savefig(f"./plots/{plotname}")
    print(f"Plot saved as {plotname}")
    with open("./out/geometry_1.out", "w") as out_1:
        for i, point in enumerate(mat_1_x):
            out_1.write("{0}, {1}\n".format(point, mat_1[i]))
    with open("./out/geometry_2.out", "w") as out_2:
        for i, point in enumerate(mat_2_x):
            out_2.write("{0}, {1}\n".format(point, mat_2[i]))

if __name__ == '__main__':
    main()
