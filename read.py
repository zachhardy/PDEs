# from readers import NeutronicsSimulationReader
from python_utils import KEigenvalueNeutronicsReader

import matplotlib.pyplot as plt


sim = KEigenvalueNeutronicsReader("Test/TWIGL/result.data")
sim.read()

sim.plot_flux_moment()
# sim.plot_precursors_species(0)
