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
    pattern = re.compile(r'\s*(\S+)\s+(\S+)\s*')
    with open('./out/geometry_test.out', mode='r') as datafile:
        for line in datafile:
            if re.match(pattern, line):
                x_value = float(re.match(pattern, line).group(1))
                mat = int(re.match(pattern, line).group(2))
                x_values = np.append(x_values, x_value)
                material = np.append(material, mat)
                if mat == 1:
                    mat_1_x = np.append(mat_1_x, x_value)
                    mat_1 = np.append(mat_1, mat)
                if mat == 2:
                    mat_2_x = np.append(mat_2_x, x_value)
                    mat_2 = np.append(mat_2, mat)
    plt.scatter(mat_1_x, mat_1, marker='_', color='b', s=5)
    plt.scatter(mat_2_x, mat_2, marker='_', color='r', s=5)
    print(mat_1_x)
    plt.title("Sample Markovian Geometry")
    plt.xlabel("x")
    plt.ylabel("Material Number")
    plt.ylim(0, 5)
    plt.xlim(0, 10)
    plt.tight_layout()
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
