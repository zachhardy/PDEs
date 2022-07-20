# from readers import NeutronicsSimulationReader
from python_utils import TransientNeutronicsReader

import matplotlib.pyplot as plt


sim = TransientNeutronicsReader("Test/LRA/outputs")
sim.read()

# sim.plot_flux_moment(groups=0, times=[0.0, 1.44])
sim.plot_temperature('peak', False)
plt.show()
