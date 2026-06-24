import matplotlib.pyplot as plt
import numpy as np
import math

def plot_data(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    data = []
    for line in lines:
        values = list(map(float, line.split()))
        # The entries with an imaginary part too small correspond to the bandgap
        # region, where there is no phase shift.
        if abs(values[3]) > 1e-14:
            data.append(values)

    omega = [row[0] for row in data]
    re_data = [row[2] for row in data]
    im_data = [row[3] for row in data]
    wvect = [math.atan2(y, x)/2/math.pi for x, y in zip(re_data, im_data)]
    plt.scatter(wvect, omega, marker=".")
    plt.title('Photonic band structure GaAs/air multilayer')
    plt.xlabel('Wave vector ka/2π')
    plt.ylabel('Frequency ωa/2πc')
    plt.show()



plot_data("bande.dat")