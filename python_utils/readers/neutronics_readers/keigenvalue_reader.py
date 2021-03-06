import os

import matplotlib.pyplot as plt

from .steadystate_reader import SteadyStateNeutronicsReader

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 12


class KEigenvalueNeutronicsReader(SteadyStateNeutronicsReader):

    def __init__(self, path):

        super().__init__(path)

        self.k_eff = 0.0

    def read(self):
        """
        Read a neutronics data files to populate simulation data.
        """
        SteadyStateNeutronicsReader.read(self)

        with open(f"{os.path.dirname(self.path)}/k_eff.txt", "r") as file:
            self.k_eff = float(file.readline())



