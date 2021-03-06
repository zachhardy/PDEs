import os
import sys
import itertools
import time
import copy

import numpy as np


def setup_directory(path):
    """
    Prepare a directory for simulation outputs.

    :param path: The path to the directory.
    :type path: str
    """
    if os.path.isdir(path):
        os.system(f"rm -r {path}")
    os.makedirs(path)


def define_range(reference, variance, n):
    """
    Define sample points for a parameter given the nominal value,
    the variance, and the number of samples desired.

    :param reference: The nominal parameter value.
    :type reference: float
    :param variance: The variance of the parameter.
    :type variance: float
    :param n: The number of samples to define within the range.
    :return: The samples to use in the parameter study. This routine
        defines a range which spans one times the variance in either
        direction about the nominal value.
    :rtype: numpy.ndarray
    """
    samples = np.linspace(-1.0, 1.0, n)
    return reference*(1.0 + variance * samples)


def parameter_study(problem_name, study):
    """
    Define and run a parameter study.

    :param problem_name: The name of the problem to run.
    :type problem_name: str {'Sphere3g', 'TWIGL', 'LRA'}
    :param study: The study number
    :return:
    """

    path = os.path.dirname(os.path.abspath(__file__))
    path = "/".join(path.split("/")[:-1])
    path = os.path.join(path, "Problems")
    study = int(study)

    ##################################################
    # Define the parameter study
    ##################################################

    # Three group sphere problem
    if problem_name == "Sphere3g":
        path = os.path.join(path, problem_name)
        path_to_exe = os.path.join(path, "bin", "Sphere3g_Transient")

        # Define the parameter
        parameters = {}
        if study == 0:
            parameters['radius'] = define_range(6.1612, 0.025, 31)
        elif study == 1:
            parameters['density'] = define_range(0.05134325, 0.025, 31)
        elif study == 2:
            parameters['size'] = define_range(6.0, 0.01, 6)
            parameters['density'] = define_range(0.05, 0.01, 6)
        elif study == 3:
            parameters['radius'] = define_range(6.0, 0.02, 5)
            parameters['density'] = define_range(0.05, 0.005, 5)
            parameters['scatter'] = define_range(1.46, 0.1, 5)
        else:
            raise NotImplementedError(
                f"Invalid study number. For the {problem_name} problem,"
                f"the maximum input value is 3.")

    # TWIGL problem
    elif problem_name == "TWIGL":
        path = os.path.join(path, problem_name)
        path_to_exe = os.path.join(path, "bin", "TWIGL")

        magnitude = 0.97667 - 1.0
        duration = 0.2
        scatter = 0.01

        # Define the parameter
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

        path = os.path.join(path, problem_name)
        path_to_exe = os.path.join(path, "bin", "LRA")

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
    output_path = f"{path}/studies/{keys[0]}"
    for k, key in enumerate(keys[1:]):
        output_path += f"_{key}"
    setup_directory(output_path)

    # Save the parameters to a file
    param_path = f"{output_path}/params.txt"
    np.savetxt(param_path, values, fmt='%.8e')

    ##################################################
    # Run the reference problem
    ##################################################
    sim_path = os.path.join(output_path, "reference")
    setup_directory(sim_path)

    cmd = f"{path_to_exe} output_directory={sim_path} >> {sim_path}/log.txt"
    os.system(cmd)

    ##################################################
    # Run the parameter study
    ##################################################
    t_avg = 0.0
    for n, params in enumerate(values):

        # Setup output path
        sim_path = os.path.join(output_path, str(n).zfill(3))
        setup_directory(sim_path)

        cmd = f"{path_to_exe} "
        for k, key in enumerate(keys):
            cmd += f"{key}={params[k]} "
        cmd += f"output_directory={sim_path} >> {sim_path}/log.txt"

        msg = f"="*50 + f"\nRunning Simulation {n}\n" + f"="*50
        for k, key in enumerate(keys):
            s = " ".join([w.capitalize() for w in key.split("_")])
            msg += f"\n{s:<{max_len}}:\t{params[k]:<5.3e}"
        print()
        print(msg)

        t_start = time.time()
        os.system(cmd)
        sim_time = time.time() - t_start

        print()
        print(f"Simulation Time = {sim_time:.3f} s")


if __name__ == "__main__":
    for i in range(7):
        parameter_study("LRA", i)

    # parameter_study(*sys.argv[1:])


