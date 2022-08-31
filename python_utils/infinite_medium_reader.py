import os
import struct

import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import physical_constants as constants

k = constants["Boltzmann constant in eV/K"][0]


def read_psi(file_path):
    file_path = os.path.abspath(file_path)
    if not os.path.isfile(file_path):
        raise FileNotFoundError

    def read_unsigned_int(f):
        return struct.unpack('I', f.read(4))[0]

    def read_double(f):
        return struct.unpack('d', f.read(8))[0]

    with open(file_path, mode='rb') as file:
        file.read(249)  # skip the header

        n_angles = read_unsigned_int(file)
        n_groups = read_unsigned_int(file)

        x = np.zeros((n_angles, n_groups))

        for n in range(n_angles):
            for g in range(n_groups):
                x[n][g] = read_double(file)
    return x


def read_phi(file_path):
    file_path = os.path.abspath(file_path)
    if not os.path.isfile(file_path):
        raise FileNotFoundError

    def read_unsigned_int(f):
        return struct.unpack('I', f.read(4))[0]

    def read_double(f):
        return struct.unpack('d', f.read(8))[0]

    with open(file_path, mode='rb') as file:
        file.read(249)  # skip the header

        n_moments = read_unsigned_int(file)
        n_groups = read_unsigned_int(file)

        x = np.zeros((n_moments, n_groups))

        for ell in range(n_moments):
            for g in range(n_groups):
                x[ell][g] = read_double(file)
    return x


def maxwell(E, T):
    assert isinstance(E, np.ndarray)
    print(k, T)
    return (2.0 * np.sqrt(E / np.pi) *
            (1.0 / (k * T)) ** (3 / 2) *
            np.exp(-E / (k * T)))


if __name__ == "__main__":
    path = os.path.abspath(os.path.dirname(__file__))

    psi = read_psi(os.path.join(path, "..", "psi.data"))
    phi = read_phi(os.path.join(path, "..", "phi.data"))
    E_bounds = np.loadtxt(os.path.join(path, "..", "e_bounds.txt"))

    T_therm = 293.0  # in K or 2.530064236e-08 in MeV
    plt.semilogx(E_bounds, maxwell(E_bounds, T_therm), "")
    plt.show()
