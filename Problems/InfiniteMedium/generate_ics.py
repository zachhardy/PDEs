import os
import sys

import numpy as np


def maxwellian(E, kT):
    """
    Evaluate a Maxwellian distribution in energy.

    Parameters
    ----------
    E : float, The energy in MeV.
    kT : float, The temperature parameter in MeV.
    """
    return 2.0*np.pi / (np.pi * kT)**1.5 * \
        np.sqrt(E) * np.exp(-E / kT)


def read_energy_bins(p):
    """
    Read the energy bins from a cross-section file.

    Parameters
    ----------
    p : The path the the cross-section file.
    """
    x = []
    with open(p, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if "e_bounds" in line:
                i += 1
                line = lines[i]
                entries = line.split()

                while not entries[0].isalpha():
                    for entry in entries:
                        x.append(float(entry))

                    i += 1
                    line = lines[i]
                    entries = line.split()
                break
    return x


def generate_ics(outfile, E_b, kT):
    """
    Generate an initial condition file for the given energy
    groups and temperature parameter 'kT'

    Parameters
    ----------
    outfile : A path to the file to write to.
    E_b : The energy group boundaries in MeV
    kT : The temperature parameter in MeV
    """
    # Print summary
    print()
    print("##################################################")
    print(f"# Output file set to: {outfile}")
    print(f"# kT set to {kT} MeV")
    print("##################################################")

    with open(outfile, 'w') as f:
        f.write(f"##################################################\n"
                f"# Maxwellian Initial Condition\n"
                f"# Number of Groups: {len(E_bounds) - 1}\n"
                f"# kT: {kT} MeV\n"
                f"##################################################\n")

        for g in range(len(E_bounds) - 1):
            E = 0.5*(E_bounds[g] + E_bounds[g + 1])
            M = maxwellian(E, kT)
            if M > np.finfo(float).eps:
                f.write(f"{g:<5}  {M:.10e}\n")


if __name__ == "__main__":
    assert len(sys.argv) > 2

    # Get the cross-sections path
    xs_path = os.path.abspath(sys.argv[1])
    assert os.path.isfile(xs_path)
    E_bounds = read_energy_bins(xs_path)

    # Define output directors
    filepath = os.path.abspath(os.path.dirname(__file__))
    outdir = os.path.join(filepath, "ics")
    outdir = os.path.join(outdir, f"{len(E_bounds) - 1}g")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Setup the kT parameters
    kTs = [2.524878645766e-08, 0.75, 1.0, 1.5]
    names = ["room", "0", "1", "2"]
    assert len(kTs) == len(names)

    # Construct initial conditions
    for i in range(len(kTs)):
        filename = "maxwell_" + names[i] + ".txt"
        outfile = os.path.join(outdir, filename)
        generate_ics(outfile, E_bounds, kTs[i])
