import os
import numpy as np


class SimulationReader:

    def __init__(self, path):
        self.path = os.path.abspath(path)

    def read(self):
        raise NotImplementedError

    def init_storage(self):
        raise NotImplementedError

    def default_keys(self):
        raise NotImplementedError

    def get(self, key):
        raise NotImplementedError

    def clear(self):
        self.__init__(self.path)

    def create_matrix(self, keys=None):
        if isinstance(keys, str):
            keys = [keys]
        elif keys is None:
            keys = self.default_keys()

        matrix = self.get(keys[0])
        for i, key in enumerate(keys[1:]):
            matrix = np.hstack((matrix, self.get(key)))
        return matrix

    def create_vector(self, keys=None):
        data = self.create_matrix(keys)
        return data.flatten()
