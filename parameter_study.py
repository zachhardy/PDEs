import os
import sys
import itertools
import time
import copy

import numpy as np


def setup_directory(path):
    """
    Prepare a directory for simulation outputs.

    Parameters
    ----------
    path : str, The path to the directory.
    """
    if os.path.isdir(path):
        os.system(f"rm -r {path}")
    os.makedirs(path)


def define_range(reference, variance, n):
    """
    Define sample points for a parameter given the nominal value,
    the variance, and the number of samples desired. This routine
    generates uniformly spaced samples within plus or minus the
    specified variance about the specified reference value.

    Parameters
    ----------
    reference : float, The nominal parameter value.
    variance : float, The variance of the parameter.
    n : int, The number of samples to generate.

    Returns
    -------
    numpy.ndarray : The samples to use in the parameter study.
    """
    samples = np.linspace(-1.0, 1.0, n)
    return reference*(1.0 + variance * samples)


def parameter_study(problem_name, study):
    """
    Define and run a parameter study.

    Parameters
    ----------
    problem_name : str {'Sphere3g', 'TWIGL', 'LRA'}
        The problem to run a parameter study for.
    study : int, The pre-defined parameter study to run.
    """

    path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(path, "Problems", problem_name)
    path_to_exe = os.path.join(path, "bin", problem_name)

    study = int(study)

    ##################################################
    # Define the parameter study
    ##################################################

    # Three group sphere problem
    if problem_name == "Sphere3g":

        # Default value
        radius = 6.0
        density = 0.05
        sig_s01 = 1.46

        # Define the parameter space
        parameters = {}
        if study == 0:
            parameters['radius'] = define_range(6.1612, 0.025, 21)
        elif study == 1:
            parameters['density'] = define_range(0.05134325, 0.025, 21)
        elif study == 2:
            parameters['scatter'] = define_range(sig_s01, 0.2, 21)
        elif study == 3:
            parameters['radius'] = define_range(radius, 0.01, 6)
            parameters['density'] = define_range(density, 0.01, 6)
        elif study == 4:
            parameters['radius'] = define_range(radius, 0.01, 6)
            parameters['scatter'] = define_range(sig_s01, 0.1, 6)
        elif study == 5:
            parameters['density'] = define_range(density, 0.01, 6)
            parameters['scatter'] = define_range(sig_s01, 0.1, 6)
        elif study == 6:
            parameters['radius'] = define_range(radius, 0.02, 4)
            parameters['density'] = define_range(density, 0.005, 4)
            parameters['scatter'] = define_range(sig_s01, 0.1, 4)

        else:
            raise NotImplementedError(
                f"Invalid study number. For the {problem_name} problem,"
                f"the maximum input value is 3.")

    # Infinite slab problem
    elif problem_name == "InfiniteSlab":

        # Default values
        magnitude = -0.01
        duration = 1.0
        interface = 40.0

        # Define the parameter space
        parameters = {}
        if study == 0:
            parameters['magnitude'] = define_range(magnitude, 0.2, 21)
        elif study == 1:
            parameters['duration'] = define_range(duration, 0.2, 21)
        elif study == 2:
            parameters['interface'] = define_range(interface, 0.05, 21)
        elif study == 3:
            parameters['magnitude'] = define_range(magnitude, 0.2, 6)
            parameters['duration'] = define_range(duration, 0.2, 6)
        elif study == 4:
            parameters['magnitude'] = define_range(magnitude, 0.1, 6)
            parameters['interface'] = define_range(interface, 0.025, 6)
        elif study == 5:
            parameters['duration'] = define_range(duration, 0.2, 6)
            parameters['interface'] = define_range(interface, 0.025, 6)
        elif study == 6:
            parameters['magnitude'] = define_range(magnitude, 0.05, 4)
            parameters['duration'] = define_range(duration, 0.05, 4)
            parameters['interface'] = define_range(40.0, 0.025, 4)
        else:
            raise NotImplementedError(
                f"Invalid study number. For the {problem_name} problem,"
                f"the maximum input value is 6.")

    # TWIGL problem
    elif problem_name == "TWIGL":

        # Default values
        magnitude = 0.97667 - 1.0
        duration = 0.2
        scatter = 0.01

        # Define the parameter space
        parameters = {}
        if study == 0:
            parameters['magnitude'] = define_range(magnitude, 0.2, 21)
        elif study == 1:
            parameters['duration'] = define_range(duration, 0.2, 21)
        elif study == 2:
            parameters['scatter'] = define_range(scatter, 0.25, 21)
        elif study == 3:
            parameters['magnitude'] = define_range(magnitude, 0.02, 6)
            parameters['duration'] = define_range(duration, 0.25, 6)
        elif study == 4:
            parameters['magnitude'] = define_range(magnitude, 0.2, 6)
            parameters['scatter'] = define_range(scatter, 0.25, 6)
        elif study == 5:
            parameters['duration'] = define_range(duration, 0.2, 6)
            parameters['scatter'] = define_range(scatter, 0.25, 6)
        elif study == 6:
            parameters['magnitude'] = define_range(magnitude, 0.2, 4)
            parameters['duration'] = define_range(duration, 0.2, 4)
            parameters['scatter'] = define_range(scatter, 0.25, 4)
        else:
            raise NotImplementedError(
                f"Invalid study number. For the {problem_name} problem,"
                f"the maximum input value is 6.")

    # LRA benchmark problem
    elif problem_name == "LRA":

        # Default values
        magnitude = 0.8787631 - 1.0
        duration = 2.0
        feedback = 3.034e-3

        # Define parameter space
        parameters = {}
        if study == 0:
            parameters['magnitude'] = define_range(magnitude, 0.025, 21)
        elif study == 1:
            parameters['duration'] = define_range(duration, 0.05, 21)
        elif study == 2:
            parameters['feedback'] = define_range(feedback, 0.05, 21)
        elif study == 3:
            parameters['magnitude'] = define_range(magnitude, 0.025, 6)
            parameters['duration'] = define_range(duration, 0.05, 6)
        elif study == 4:
            parameters['magnitude'] = define_range(magnitude, 0.025, 6)
            parameters['feedback'] = define_range(feedback, 0.05, 6)
        elif study == 5:
            parameters['duration'] = define_range(duration, 0.05, 6)
            parameters['feedback'] = define_range(feedback, 0.05, 6)
        elif study == 6:
            parameters['magnitude'] = define_range(magnitude, 0.025, 4)
            parameters['duration'] = define_range(duration, 0.05, 4)
            parameters['feedback'] = define_range(feedback, 0.05, 4)
        else:
            raise NotImplementedError(
                f"Invalid study number. For the {problem_name} problem,"
                f"the maximum input value is 6.")

    else:
        raise NotImplementedError("Invalid problem name.")

    keys = list(parameters.keys())
    max_len = np.max([len(key) for key in keys])

    values = np.array(list(itertools.product(*parameters.values())))
    values = np.round(values, 10)

    ##################################################
    # Setup the output paths
    ##################################################

    # Define the path to the output directory
    output_path = f"{path}/parametric/{keys[0]}"
    for k, key in enumerate(keys[1:]):
        output_path += f"_{key}"
    setup_directory(output_path)

    # Save the parameters to a file
    param_path = f"{output_path}/params.txt"
    header = " ".join([f"{key:<13} " for key in keys])
    np.savetxt(param_path, values, fmt='%.8e', header=header)

    ##################################################
    # Run the reference problem
    ##################################################

    # sim_path = os.path.join(output_path, "reference")
    # setup_directory(sim_path)
    #
    # cmd = f"{path_to_exe} output_directory={sim_path} "
    # cmd += f"xs_directory={path}/xs >> {sim_path}/log.txt"
    # os.system(cmd)

    ##################################################
    # Run the parameter study
    ##################################################

    total_time = 0.0
    for n, params in enumerate(values):

        # Setup output path
        sim_path = os.path.join(output_path, str(n).zfill(3))
        setup_directory(sim_path)

        cmd = f"{path_to_exe} "
        for k, key in enumerate(keys):
            cmd += f"{key}={params[k]} "
        cmd += f"output_directory={sim_path} "
        cmd += f"xs_directory={path}/xs >> {sim_path}/log.txt"

        msg = f"="*50 + f"\nRunning Simulation {n}\n" + f"="*50
        for k, key in enumerate(keys):
            s = " ".join([w.capitalize() for w in key.split("_")])
            msg += f"\n{s:<{max_len}}:\t{params[k]:<5.3e}"
        print()
        print(msg)

        t_start = time.time()
        os.system(cmd)
        sim_time = time.time() - t_start
        total_time += sim_time

        print()
        print(f"Simulation Time = {sim_time:.3f} s")

    print()
    print(f"Average Simulation Time = {total_time/len(values):.3e} s")
    print(f"Total Parameter Study Time = {total_time:.3e} s")


if __name__ == "__main__":

    ############################################################
    # Maximum Parameter Study Number:
    #   Sphere3g     - 6
    #   InfiniteSlab - 6
    #   TWIGL        - 6
    #   LRA          - 6
    ############################################################

    # If only
    if len(sys.argv) == 2:
        name = sys.argv[1]
        if name == "Sphere3g":
            max_study = 6
        elif name == "InfiniteSlab":
            max_study = 6
        elif name == "TWIGL":
            max_study = 6
        elif name == "LRA":
            max_study = 6
        else:
            raise NotImplementedError("Invalid problem name.")

        check = input("Are you sure you want to run every parameter "
                      f"study for {name}? [y/n] ")
        if "y" in check:
            for i in range(max_study + 1):
                parameter_study(name, i)
        if "n" in check:
            print("Terminating program.")
            exit(0)

    elif len(sys.argv) == 3:
        parameter_study(*sys.argv[1:])

    else:
        raise AssertionError("Invalid command line inputs.")


