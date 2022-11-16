import os
import numpy as np

from .simulation_reader import NeutronicsSimulationReader


class NeutronicsDatasetReader:
    """
    A utility class for reading collections of neutronics simulation data.
    """

    def __init__(self, path):
        """
        Create a neutronics dataset readers.

        Parameters
        ----------
        path : str, Path to collection of simulation output data.
        """
        if not os.path.isdir(path):
            raise NotADirectoryError(f"{path} is not a valid directory.")

        self._path = os.path.abspath(path)

        self.simulations: list[NeutronicsSimulationReader] = []
        self.parameters: np.ndarray = None

    @property
    def n_simulations(self):
        """
        Return the number of simulations in the dataset.

        Returns
        -------
        int
        """
        return len(self.simulations)

    @property
    def n_parameters(self):
        """
        Return the number of parameters which describe the dataset.

        Returns
        -------
        int
        """
        return len(self.parameters)

    @property
    def n_snapshots(self):
        """
        Return the number of temporal snapshots.

        Returns
        -------
        int
        """
        return self.simulations[0].n_snapshots

    @property
    def times(self):
        """
        Return the snapshot times.

        Returns
        -------
        list[float]
        """
        return self.simulations[0].times

    @property
    def dimension(self):
        """
        Return the spatial dimension.

        Returns
        -------
        int
        """
        return self.simulations[0].dimension

    @property
    def n_materials(self):
        """
        Return the number of materials in the problem.

        Returns
        -------
        int
        """
        return self.simulations[0].n_materials

    @property
    def n_cells(self):
        """
        Return the number of cells.

        Returns
        -------
        int
        """
        return self.simulations[0].n_cells

    @property
    def centroids(self):
        """
        Return the centroid locations.

        Returns
        -------
        numpy.ndarray
        """
        return self.simulations[0].centroids

    @property
    def n_nodes(self):
        """
        Return the number of nodes.

        Returns
        -------
        int
        """
        return self.simulations[0].n_nodes

    @property
    def nodes(self):
        """
        Return the ndoe locations.

        Returns
        -------
        numpy.ndarray
        """
        return self.simulations[0].nodes

    @property
    def n_moments(self):
        """
        Return the number of flux moments.

        Returns
        -------
        int
        """
        return self.simulations[0].n_moments

    @property
    def n_groups(self):
        """
        Return the number of energy groups.

        Returns
        -------
        int
        """
        return self.simulations[0].n_groups

    @property
    def max_precursors(self):
        """
        Return the maximum number of precursors per material.

        Returns
        -------
        int
        """
        return self.simulations[0].max_precursors

    def read(self):
        """
        Read in the dataset.

        Returns
        -------
        NeutronicsDatasetReader
        """
        self.__init__(self._path)

        # Loop through simulations
        for i, entry in enumerate(sorted(os.listdir(self._path))):
            path = os.path.join(self._path, entry)
            if entry == "params.txt":
                params = np.loadtxt(path)
                if params.ndim == 1:
                    params = np.atleast_2d(params).T
                self.parameters = params

            elif os.path.isdir(path) and "reference" not in entry:
                reader = NeutronicsSimulationReader(path).read()
                if self.n_simulations > 1:
                    self._check_compatibility(reader)
                self.simulations.append(reader)
        return self

    def _check_compatibility(self, reader):
        """
        Ensure that the simulations are identical.

        Parameters
        ----------
        reader : NeutronicsSimulationReader
        """
        err_msg = "Simulation setup data does not agree."
        if (reader.n_snapshots != self.n_snapshots or
                reader.dimension != self.dimension or
                reader.n_materials != self.n_materials or
                reader.n_cells != self.n_cells or
                reader.n_nodes != self.n_nodes or
                reader.n_moments != self.n_moments or
                reader.n_groups != self.n_groups or
                reader.max_precursors != self.max_precursors):
            raise AssertionError(err_msg)

