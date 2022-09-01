import os
import struct
import sys

import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import physical_constants as constants

k = constants["Boltzmann constant in eV/K"][0]/1.0e6


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
    return (2.0 * np.pi / (np.pi * k * T)**(3/2) *
            np.sqrt(E) * np.exp(- E / (k * T)))


if __name__ == "__main__":
    assert len(sys.argv)
    assert os.path.isdir(sys.argv[1])
    path = os.path.abspath(sys.argv[1])

    phi, E_bounds = [], []
    for entry in sorted(os.listdir(sys.argv[1])):
        entry_path = os.path.join(path, entry)

        if "e_bounds.txt" in entry_path:
            E_bounds = np.loadtxt(os.path.join(entry_path))[:, 1]

        elif os.path.isdir(entry_path):
            for file in sorted(os.listdir(entry_path)):
                fpath = os.path.join(entry_path, file)
                if "sflux" in fpath:
                    phi.append(read_phi(fpath))

    phi = np.array(phi)
    E_avg = np.zeros(len(E_bounds) - 1)
    dE = np.zeros(len(E_bounds) - 1)
    for g in range(len(E_bounds) - 1):
        E_avg[g] = 0.5 * (E_bounds[g] + E_bounds[g + 1])
        dE[g] = E_bounds[g] - E_bounds[g + 1]
        phi[:, :, g] *= E_avg[g]/dE[g]

    # Initialize the figure
    plt.figure()
    plt.xlabel("Energy (MeV)", fontsize=12)
    plt.ylabel(rf"$\phi(E) \bar{{E}} / \Delta E$", fontsize=12)
    plt.grid(True)

    # Plot the infinite medium results at various times\
    for i in range(0, len(phi), len(phi)//5):
        vals = phi[i][0]
        plt.semilogx(E_avg, phi[i][0], label=f"n = {int(i)}")
    plt.legend()

    savepath = os.path.abspath(os.path.dirname(__file__))
    plt.savefig(f"{savepath}/maxwellians.pdf")
    plt.show()
