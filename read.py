# from readers import NeutronicsSimulationReader
from python_utils import KEigenvalueNeutronicsReader

import matplotlib.pyplot as plt


sim = KEigenvalueNeutronicsReader("Test/LRA/outputs/result.data")
sim.read()

sim.plot_material_ids()
# sim.plot_precursors_species(0)
