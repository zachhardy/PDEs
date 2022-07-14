# from readers import NeutronicsSimulationReader
from python_utils import NeutronicsSimulationReader

import matplotlib.pyplot as plt


sim = NeutronicsSimulationReader("outputs")
sim.read()

sim.plot_flux_moment()
# sim.plot_precursors_species(0)
