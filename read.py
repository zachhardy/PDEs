# from readers import NeutronicsSimulationReader
from python_utils.readers import TransientNeutronicsReader

import matplotlib.pyplot as plt


sim = TransientNeutronicsReader("Problems/Sphere3g/studies/radius/reference")
sim.read()

sim.plot_power('total')
plt.show()
