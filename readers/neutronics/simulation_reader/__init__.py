import os
import struct
import numpy as np
import matplotlib.pyplot as plt


class NeutronicsSimulationReader:
    """
    A utility class for reading neutronics simulation data.
    """

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 12

    # Energy released per fission event in J
    energy_per_fission: float = 3.204e-11

    from ._plot import plot_flux_moment
    from ._plot import plot_power_profile
    from ._plot import plot_temperature_profile

    from ._plot import plot_power
    from ._plot import plot_fuel_temperature

    def __init__(self, path):
        """
        Create a neutronics simulation readers.

        This constructor sets the path to the output directory
        and specifies which files to read.

        Parameters
        ----------
        path : str, Path to the output binaries.
        """

        if not os.path.isdir(path):
            raise NotADirectoryError(f"{path} is not a valid directory.")

        self._path = os.path.abspath(path)

        self.nodes: np.ndarray = None
        self.centroids: np.ndarray = None
        self.material_ids: np.ndarray = None

        self.times: np.ndarray = None

        self.powers: np.ndarray = None
        self.peak_power_densities: np.ndarray = None
        self.average_power_densities: np.ndarray = None

        self.average_fuel_temperatures: np.ndarray = None
        self.peak_fuel_temperatures: np.ndarray = None

        self.flux_moments: np.ndarray = None
        self.precursors: np.ndarray = None
        self.fission_rates: np.ndarray = None
        self.temperatures: np.ndarray = None

    @property
    def n_snapshots(self):
        """
        Return the number of time steps.

        Returns
        -------
        int
        """
        return len(self.times)

    @property
    def n_cells(self):
        """
        Return the number of cells.

        Returns
        -------
        int
        """
        return len(self.centroids)

    @property
    def n_nodes(self):
        """
        Return the number of nodes.

        Returns
        -------
        int
        """
        return len(self.nodes)

    @property
    def dimension(self):
        """
        Return the dimension of the spatial domain.

        Returns
        -------
        int
        """
        if all(self.centroids[:, :2].ravel() == 0.0):
            return 1
        elif all(self.centroids[:, 2].ravel() == 0.0):
            return 2
        else:
            return 3

    @property
    def n_materials(self):
        """
        Return the number of unique materials.

        Returns
        -------
        int
        """
        return len(np.unique(self.material_ids))

    @property
    def n_moments(self):
        """
        Return the number of flux moments.

        Returns
        -------
        int
        """
        return len(self.flux_moments[0])

    @property
    def n_groups(self):
        """
        Return the number of energy groups.

        Returns
        -------
        int
        """
        return len(self.flux_moments[0][0])

    @property
    def max_precursors(self):
        """
        Return the maximum number of precursors per material.

        Returns
        -------
        int
        """
        return len(self.precursors[0])

    def create_flux_moment_matrix(self):
        """
        Create a matrix containing the flux moments ordered by moment, then
        group, then node.

        Returns
        -------
        numpy.ndarray
        """
        n = self.n_snapshots
        m = self.n_nodes * self.n_moments * self.n_groups
        phi = self.flux_moments.reshape((n, m))
        return phi

    def create_flux_moment_vector(self):
        """
        Create a vector containing the flux moments ordered by time, then
        moment, then group, then node.

        Returns
        -------
        numpy.ndarray
        """
        return self.flux_moments.reshape((-1, 1))

    def create_power_density_matrix(self):
        """
        Create a matrix containing the power density profiles.

        Returns
        -------
        numpy.ndarray
        """
        Ef = self.energy_per_fission
        n, m = self.n_snapshots, self.n_cells
        return Ef * self.fission_rates.reshape((n, m))

    def create_power_density_vector(self):
        """
        Create a vector containing the power density profiles ordered
        by time then node.

        Returns
        -------
        numpy.ndarray
        """
        Ef = self.energy_per_fission
        return Ef * self.fission_rates.reshape((-1, 1))

    def read(self,
             parse_flux_moments=True,
             parse_precursors=True,
             parse_fission_rate=True,
             parse_temperature=True):
        """
        Read transient simulation data.

        Parameters
        ----------
        parse_flux_moments : bool, Flag for reading multi-group scalar fluxes.
        parse_precursors : bool, Flag for reading precursor concentrations.
        parse_fission_rate : bool, Flag for reading fission rates.
        parse_temperature : bool, Flag for reading temperatures.

        Returns
        -------
        NeutronicsSimulationReader
        """
        self.__init__(self._path)

        # Parse summary file
        tmp = np.loadtxt(f"{self._path}/summary.txt")
        self.times, self.powers = tmp[:, 0], tmp[:, 1]
        self.peak_power_densities = tmp[:, 2]
        self.average_power_densities = tmp[:, 3]
        self.peak_fuel_temperatures = tmp[:, 4]
        self.average_fuel_temperatures = tmp[:, 5]

        # Parse geometry file
        tmp = self.read_geometry_file(f"{self._path}/geom.data")
        self.centroids, self.nodes, self.material_ids = tmp

        # Loop through outputs
        for i, entry in enumerate(sorted(os.listdir(self._path))):
            path = os.path.join(self._path, entry)

            # Go into time step directory
            if os.path.isdir(path):
                for file in os.listdir(path):
                    filepath = os.path.join(path, file)

                    if file == "sflux.data" and parse_flux_moments:
                        phi = self.read_sflux_file(filepath)
                        if self.flux_moments is None:
                            shape = (self.n_snapshots, *phi.shape)
                            self.flux_moments = np.zeros(shape)
                        self.flux_moments[int(entry)] = phi

                    elif file == "precursors.data" and parse_precursors:
                        Cj = self.read_precursors_file(filepath)
                        if self.precursors is None:
                            shape = (self.n_snapshots, *Cj.shape)
                            self.precursors = np.zeros(shape)
                        self.precursors[int(entry)] = Cj

                    elif file == "fission_rate.data" and parse_fission_rate:
                        Sf = self.read_fission_rate_file(filepath)
                        if self.fission_rates is None:
                            shape = (self.n_snapshots, *Sf.shape)
                            self.fission_rates = np.zeros(shape)
                        self.fission_rates[int(entry)] = Sf

                    elif file == "temperature.data" and parse_temperature:
                        T = self.read_temperature_file(filepath)
                        if self.temperatures is None:
                            shape = (self.n_snapshots, *T.shape)
                            self.temperatures = np.zeros(shape)
                        self.temperatures[int(entry)] = T
        return self

    def read_geometry_file(self, file_path):
        """
        Read the geometry file for a time step.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        numpy.ndarray
            The (x, y, z) coordinates of the cell centroids.
        numpy.ndarray
            The (x, y, z) coordinates of the nodes per cell.
        numpy.ndarray
            The material IDs.
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Cannot find file {file_path}.")

        with open(file_path, 'rb') as file:
            file.read(399)

            n_cells = self.read_uint64_t(file)
            n_nodes = self.read_uint64_t(file)

            centroids = np.zeros((n_cells, 3))
            nodes = np.zeros((n_nodes, 3))
            matids = np.zeros(n_cells)

            n = 0
            for c in range(n_cells):
                self.read_uint64_t(file)
                matids[c] = self.read_unsigned_int(file)
                nodes_per_cell = self.read_unsigned_int(file)

                # Parse centroids
                for d in range(3):
                    centroids[c, d] = self.read_double(file)

                # Parse nodes on the cell
                for i in range(nodes_per_cell):
                    for d in range(3):
                        nodes[n, d] = self.read_double(file)
                    n += 1
        return centroids, nodes, matids

    def read_sflux_file(self, file_path):
        """
        Read the scalar flux moments file for a time step.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        numpy.ndarray
            The flux moments per moment per group per node.
            The shape of the output is (n_moments, n_groups, n_nodes).
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Cannot find file {file_path}.")

        with open(file_path, 'rb') as file:
            file.read(249)

            n_nodes = self.read_uint64_t(file)
            n_moments = self.read_unsigned_int(file)
            n_groups = self.read_unsigned_int(file)

            # Parse this time step
            phi = np.zeros((n_moments, n_groups, n_nodes))
            for i in range(n_nodes):
                for m in range(n_moments):
                    for g in range(n_groups):
                        node = self.read_uint64_t(file)
                        moment = self.read_unsigned_int(file)
                        group = self.read_unsigned_int(file)
                        value = self.read_double(file)
                        phi[moment][group][node] = value
        return phi

    def read_precursors_file(self, file_path):
        """
        Read the precursors file for a time step.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        numpy.ndarray
            The precursors per species per cell.
             The shape of the output is (max_precursors, n_cells).
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Cannot find file {file_path}.")

        with open(file_path, 'rb') as file:
            file.read(249)

            n_cells = self.read_uint64_t(file)
            max_precursors = self.read_unsigned_int(file)

            # Parse this time step
            Cj = np.zeros((max_precursors, n_cells))
            for i in range(n_cells):
                for j in range(max_precursors):
                    cell = self.read_uint64_t(file)
                    precursor = self.read_unsigned_int(file)
                    value = self.read_double(file)
                    Cj[precursor][cell] = value
        return Cj

    def read_fission_rate_file(self, file_path):
        """
        Read the fission rate file for a time step.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        numpy.ndarray
            The fission rate per cell.
            The shape of the output is (n_cells,).
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Cannot find file {file_path}.")

        with open(file_path, 'rb') as file:
            file.read(199)

            n_cells = self.read_uint64_t(file)

            # Parse this time step
            Sf = np.zeros(n_cells)
            for i in range(n_cells):
                cell = self.read_uint64_t(file)
                value = self.read_double(file)
                Sf[cell] = value
        return Sf

    def read_temperature_file(self, file_path):
        """
        Read the temperature file for a time step.

        Parameters
        ----------
        file_path : str

        Returns
        -------
        numpy.ndarray
            The temperature per cell.
            The shape of the output is (n_cells,).
        """
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"Cannot find file {file_path}.")

        with open(file_path, 'rb') as file:
            file.read(199)

            n_cells = self.read_uint64_t(file)

            # Parse this time step
            T = np.zeros(n_cells)
            for i in range(n_cells):
                cell = self.read_uint64_t(file)
                value = self.read_double(file)
                T[cell] = value
        return T

    def _interpolate(self, times, data):
        """
        Return interpolated data for a specified time.

        Parameters
        ----------
        self : NeutronicsSimulationReader
        times : list[float]
        data : numpy.ndarray

        Returns
        -------
        numpy.ndarray
            The values of the specified data at the specified times.
        """
        dt = self.times[1] - self.times[0]
        vals = np.zeros((len(times), *data[0].shape))
        for t, time in enumerate(times):
            i = [int(np.floor(time / dt)), int(np.ceil(time / dt))]
            w = [i[1] - time / dt, time / dt - i[0]]
            if i[0] == i[1]:
                w = [1.0, 0.0]
            vals[t] = w[0] * data[i[0]] + w[1] * data[i[1]]
        return vals

    @staticmethod
    def read_double(file) -> float:
        return struct.unpack('d', file.read(8))[0]

    @staticmethod
    def read_uint64_t(file) -> int:
        return struct.unpack('Q', file.read(8))[0]

    @staticmethod
    def read_unsigned_int(file) -> int:
        return struct.unpack('I', file.read(4))[0]
