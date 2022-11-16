import os
import sys
import matplotlib.pyplot as plt

from readers import NeutronicsSimulationReader


if len(sys.argv) < 2:
    raise AssertionError("No problem specified.")

problem = sys.argv[1]
if problem not in ["Sphere3g", "InfiniteSlab", "TWIGL", "LRA"]:
    raise ValueError(f"{problem} is not a valid problem.")

path = os.path.abspath(os.path.dirname(__file__))
path = f"{path}/Problems/{problem}/outputs"
if not os.path.isdir(path):
    raise NotADirectoryError(f"{path} is not a valid directory.")


reader = NeutronicsSimulationReader(path).read()
reader.plot_power(mode="BOTH")
