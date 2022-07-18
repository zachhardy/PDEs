import os
import struct
import bisect

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import Figure
from matplotlib.pyplot import Axes

from python_utils import Point
from python_utils import SimulationReader

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 12


class TransientNeutronicsReader(SimulationReader):

    def __init__(self, path):
        if not os.path.isdir(path):
            raise NotADirectoryError
        super().__init__(path)

        # Temporal Information
        self.n_snapshots = 0
        self.times = []

        # Spatial Information
        self.dimension = 0
        self.n_cells = 0
        self.n_nodes = 0
        self.nodes_per_cell = 0
        self.n_materials = 0

        self.nodes = np.empty(0, dtype=Point)
        self.centroids = np.empty(0, dtype=Point)
        self.material_ids = []

        # Neutronics Information
        self.n_moments = 0
        self.n_groups = 0
        self.n_precursors = 0
        self.max_precursors = 0
        self.precursors_per_material = []

        # Profiles
        self.flux_moments = np.zeros(0)
        self.precursors = np.zeros(0)
        self.power_density = np.zeros(0)
        self.temperature = np.zeros(0)

        # Scalars
        self.power = np.zeros(0)
        self.peak_power_density = np.zeros(0)
        self.average_power_density = np.zeros(0)

        self.peak_fuel_temperature = np.zeros(0)
        self.average_fuel_temperature = np.zeros(0)

    def read(self, skip=1):
        """
        Read a neutronics data files to populate simulation data.

        :param skip: The interval to parse output files. If skip is 1,
            every output file is parsed, otherwise every @p skip
            file is parsed. Default is 1.
        :type skip: int, optional
        """
        self.clear()

        def read_uint64_t(f):
            return struct.unpack('Q', f.read(8))[0]

        def read_unsigned_int(f):
            return struct.unpack('I', f.read(4))[0]

        def read_double(f):
            return struct.unpack('d', f.read(8))[0]

        files = sorted(os.listdir(self.path))[::skip]
        self.n_snapshots = len(files)

        # Loop over files
        for n, snapshot in enumerate(files):
            path = os.path.join(self.path, snapshot)

            # Open the snapshot file
            with open(path, mode='rb') as file:
                file.read(1399)  # skip over header

                step = read_uint64_t(file)
                n_data_blocks = read_unsigned_int(file)

                # Set attributes if this is the first snapshot
                if n == 0:
                    self.n_cells = read_uint64_t(file)
                    self.n_nodes = read_uint64_t(file)

                    self.n_moments = read_unsigned_int(file)
                    self.n_groups = read_unsigned_int(file)
                    self.n_precursors = read_unsigned_int(file)
                    self.max_precursors = read_unsigned_int(file)

                    self.init_storage()

                # Otherwise, ensure compatibility
                else:
                    assert read_uint64_t(file) == self.n_cells
                    assert read_uint64_t(file) == self.n_nodes

                    assert read_unsigned_int(file) == self.n_moments
                    assert read_unsigned_int(file) == self.n_groups
                    assert read_unsigned_int(file) == self.n_precursors
                    assert read_unsigned_int(file) == self.max_precursors

                self.times.append(read_double(file))

                self.power[step] = read_double(file)
                self.peak_power_density[step] = read_double(file)
                self.average_power_density[step] = read_double(file)

                self.peak_fuel_temperature[step] = read_double(file)
                self.average_fuel_temperature[step] = read_double(file)

                # Parse discretization information
                node_count = 0
                for c in range(self.n_cells):
                    cell_id = read_uint64_t(file)
                    material_id = read_uint64_t(file)
                    nodes_per_cell = read_unsigned_int(file)

                    # If on the first snapshot, save material IDs
                    if n == 0:
                        self.material_ids.append(material_id)
                        if c == 0:
                            self.nodes_per_cell = nodes_per_cell
                    else:
                        assert material_id == self.material_ids[cell_id]
                        assert nodes_per_cell == self.nodes_per_cell

                    # Parse the centroid
                    p = Point(*[read_double(file) for _ in range(3)])
                    if n == 0:
                        self.centroids = np.hstack((self.centroids, p))

                    # Parse the nodes
                    for i in range(nodes_per_cell):
                        p = Point(*[read_double(file) for _ in range(3)])
                        if n == 0:
                            self.nodes = np.hstack((self.nodes, p))

                # Set the number of materials
                self.n_materials = len(np.unique(self.material_ids))

                # Check centroids and nodes
                assert len(self.centroids) == self.n_cells
                assert len(self.nodes) == self.n_nodes

                self.centroids = np.array(self.centroids, Point)
                self.nodes = np.array(self.nodes, Point)

                # Parse flux moments
                assert read_unsigned_int(file) == 0
                phi = self.flux_moments[step]

                n_records = read_uint64_t(file)
                for i in range(n_records):
                    dof = self.map_phi_dof(
                        read_uint64_t(file),  # cell
                        read_unsigned_int(file),  # node
                        read_unsigned_int(file),  # moment
                        read_unsigned_int(file))  # group

                    phi[dof] = read_double(file)  # value

                # Parse precursors
                assert read_unsigned_int(file) == 1
                precursors = self.precursors[step]
                ppm = [set() for _ in range(self.n_materials)]

                n_records = read_uint64_t(file)
                for i in range(n_records):
                    cell = read_uint64_t(file)
                    material_id = read_uint64_t(file)
                    precursor = read_unsigned_int(file)

                    ppm[material_id].add(precursor)

                    dof = self.map_precursor_dof(cell, precursor)
                    precursors[dof] = read_double(file)

                # Define the number of precursors per material
                if n == 0:
                    for m in range(self.n_materials):
                        self.precursors_per_material.append(len(ppm[m]))

                # Parse power density
                if read_unsigned_int(file) == 2:
                    power_density = self.power_density[step]

                    n_records = read_uint64_t(file)
                    for i in range(n_records):
                        dof = read_uint64_t(file)
                        power_density[dof] = read_double(file)

                # Parse temperature
                if read_unsigned_int(file) == 3:
                    temperature = self.temperature[step]

                    n_records = read_uint64_t(file)
                    for i in range(n_records):
                        dof = read_uint64_t(file)
                        temperature[dof] = read_double(file)

        assert len(self.times) == self.n_snapshots

        # Determine spatial dimension
        if sum([p.x + p.y for p in self.nodes]) == 0.0:
            self.dimension = 1
        elif sum([p.z for p in self.nodes]) == 0.0:
            self.dimension = 2
        else:
            self.dimension = 3

    def init_storage(self):
        """
        Initialize data storage for the simulation data.
        """
        n_t, n_c, n_n = self.n_snapshots, self.n_cells, self.n_nodes
        n_m, n_g, n_p = self.n_moments, self.n_groups, self.max_precursors

        self.flux_moments = np.zeros((n_t, n_m * n_g * n_n))
        self.precursors = np.zeros((n_t, n_p * n_c))
        self.power_density = np.zeros((n_t, n_c))
        self.temperature = np.zeros((n_t, n_c))

        self.power = np.zeros(n_t)
        self.peak_power_density = np.zeros(n_t)
        self.average_power_density = np.zeros(n_t)
        self.peak_fuel_temperature = np.zeros(n_t)
        self.average_fuel_temperature = np.zeros(n_t)

    def map_phi_dof(self, cell, node, moment, group):
        """
        Return the degree of freedom in a flux moments vector corresponding
        to the provided @p cell, @p node, @p moment, and @p group.

        :param cell: The cell index
        :type cell: int
        :param node: The node index on the cell
        :type node: int
        :param moment: The flux moment index
        :type moment: int
        :param group: The group index
        :type group: int
        :rtype: int
        """
        return (cell * self.nodes_per_cell * self.n_moments * self.n_groups +
                node * self.n_moments * self.n_groups +
                moment * self.n_groups + group)

    def map_precursor_dof(self, cell, precursor):
        """
        Return the degree of freedom in a precursor vector corresponding
        to the provided @p cell and @p precursor.

        :param cell: The cell index
        :type cell: int
        :param precursor: The precursor index on the material.
        :type precursor: int
        :rtype: int
        """
        return cell * self.max_precursors + precursor

    def get(self, key=None):
        """
        Return the specified quantity based on @p key

        :param key: A string identifier for the desired quantity.
        :rtype: numpy.ndarray
        """
        # Get flux moments or subsets
        if "flux" in key:
            assert key.find('flux') == 0

            # Return full flux moments
            if key == "flux_moments":
                return self.flux_moments

            # Find subset of flux moments
            assert "_" in key
            args = key.split('_')[1:]

            # Get the moment index
            assert args[0].find('m') == 0
            m = int(args[0][1:])
            assert m < self.n_moments

            # Return multi-group flux moment
            if len(args) == 1:
                return self.get_flux_moment(m)

            # Get the group index
            assert args[1].find('g') == 0
            g = int(args[1][1:])
            assert g < self.n_groups

            # Return the moment group flux
            return self.get_flux_moment(m, g)

        # Get precursors or subsets
        elif "precursor" in key:
            assert key.find('precursor') == 0

            # Return full precursors
            if key == 'precursors':
                return self.precursors

            # Find subset of precursors
            assert "_" in key
            args = key.split('_')[1:]

            # Get the material index
            assert args[0].find('m') == 0
            m = int(args[0][1:])
            assert m < self.n_materials

            # Return all precursors for a material
            if len(args) == 1:
                return self.get_precursors(m)

            # Get the precursor index
            assert args[1].find('j') == 0
            j = int(args[1][1:])
            assert j < self.precursors_per_material[m]

            # Return a particular precursor
            return self.get_precursors(m, j)

        else:
            assert hasattr(self, key)
            return getattr(self, key)

    def default_keys(self):
        return ['flux_m0']

    def get_flux_moment(self, moment, groups=None, times=None):
        """
        Return a multi-group flux moment vector at all times.

        :param moment: The moment index
        :type moment: int
        :param groups: The group index.
            Default is None, which returns all groups.
        :type groups: int or list of int, optional
        :param times: The times to get the flux moment at.
            Default is None, which returns all available times.
        :type times: float or list of float, optional
        :rtype: numpy.ndarray:
        """
        assert moment < self.n_moments

        if groups is None:
            groups = list(range(self.n_groups))
        elif isinstance(groups, int):
            groups = [groups]
        assert all([g < self.n_groups for g in groups])

        if times is None:
            times = list(self.times)
        elif isinstance(times, float):
            times = [times]
        t0, tf = self.times[0], self.times[-1]
        assert all([t0 <= t <= tf for t in times])

        # Storage for the multigroup flux moment
        out = np.empty((len(times), len(groups), self.n_nodes))

        # Get the flux moments at the specified times
        phi = np.empty((len(times), self.flux_moments.shape[1]))
        for t, time in enumerate(times):
            phi[t] = self._interpolate_time(time, self.flux_moments)

        # Loop over cells, nodes, groups, and times
        for cell in range(self.n_cells):
            for node in range(self.nodes_per_cell):

                # The global node index
                pos = cell * self.nodes_per_cell + node

                for g, group in enumerate(groups):

                    # Position in the full flux moments vector
                    dof = self.map_phi_dof(cell, node, moment, group)

                    for t, time in enumerate(times):
                        out[t][g][pos] = phi[t][dof]
        return out

    def get_precursors(self, material, precursors=None, times=None):
        """
        Return all precursors that live on the specified @p material.

        :param material: The material ID
        :type material: int
        :param precursors: The precursor species index on the material.
            Default is None, which corresponds to all precursors.
        :type precursors: int or list of int, optional
        :param times: The times to get the precursors at.
            Default is None, which returns all available times.
        :type times: float or list of float, optional
        :rtype: numpy.ndarray
        """
        assert material < self.n_materials

        if precursors is None:
            precursors = list(range(self.max_precursors))
        elif isinstance(precursors, int):
            precursors = [precursors]
        n_precursors = self.precursors_per_material[material]
        assert all([j < n_precursors for j in precursors])

        if times is None:
            times = list(self.times)
        elif isinstance(times, float):
            times = [times]
        t0, tf = self.times[0], self.times[-1]
        assert all([t0 <= t <= tf for t in times])

        # Storage for the precursor species
        out = np.zeros((len(times), len(precursors), self.n_cells))

        # Get the flux moments at the specified times
        c_j = np.empty((len(times), self.precursors.shape[1]))
        for t, time in enumerate(times):
            c_j[t] = self._interpolate_time(time, self.precursors)

        # Loop over cells, precursors, and times
        for cell in range(self.n_cells):

            # If not the correct material, skip
            if self.material_ids[cell] != material:
                continue

            for j, precursor in enumerate(precursors):

                # Position in the full precursor vector
                dof = self.map_precursor_dof(cell, precursor)

                for t in range(len(times)):
                    out[t][j][cell] = c_j[t][dof]
        return out

    def _interpolate_time(self, time, vec):
        """
        Interpolate a time-dependent vector @p vec at the specifiec @p time.

        :param time: The time to evaluate the time-dependent vector.
        :type time: float
        :param vec: The time-dependent vector to interpolate on.
        :type vec: numpy.ndarray
        :return: The interpolated vector at time @p time.
        :rtype: numpy.ndarray
        """
        assert self.times[0] <= time <= self.times[-1]

        # Find lower and upper bounds
        try:
            ind = self.times.index(time)
            return vec[ind]

        except ValueError:
            ind = [bisect.bisect_left(self.times, time),
                   bisect.bisect_right(self.times, time)]

            dt = self.times[ind[1]] - self.times[ind[0]]

            # Compute interpolation weight factors

            w = [ind[1] - time / dt, time / dt - ind[0]]
            if ind[0] == ind[1]:
                w = [1.0, 0.0]

            return w[0] * vec[ind[0]] + w[1] * vec[ind[1]]

    def plot_flux_moment(self, moment=0, groups=None,
                         times=None, plot_type='group'):
        """
        Plot a multi-group flux moment at various times.

        :param moment: The moment index. Default is 0.
        :type moment: int, optional
        :param groups: The group indices. Default is all groups.
        :type: groups: int or list of int, optional
        :param times: The simulation times to plot.
            Default is the initial condition and end time.
        :type times: float or list of floats, optional
        :param plot_type: The plotting method. his parameter determines
            whether the each subplot contains all the specified groups at
            a fixed time or vice verse. The parameter specifies what is
            plotted on an individual subplot. Default is `group`.
        :type plot_type: str {`group`, `time`}, optional
        """
        assert plot_type in ['group', 'time']

        # Get the groups to plot
        if groups is None:
            groups = list(range(self.n_groups))
        elif isinstance(groups, int):
            groups = [groups]

        # Get times to plot
        if times is None:
            times = [self.times[0], self.times[-1]]
        elif isinstance(times, float):
            times = [times]
        t0, tf = self.times[0], self.times[-1]
        assert all([t0 <= t <= tf for t in times])

        # Plot 1D profiles
        if self.dimension == 1:

            # Get the coordinates
            z = [p.z for p in self.nodes]

            if plot_type == 'group':

                phi = self.get_flux_moment(moment, groups, times)

                # Create a figure for each time
                for t, time in enumerate(times):
                    fig: Figure = plt.figure(t)
                    ax: Axes = fig.add_subplot(1, 1, 1)
                    ax.set_title(fr'Time = {time:.2e} $\mu$s')
                    ax.set_xlabel(f'z (cm)')
                    ax.set_ylabel(fr"$\phi_{{{moment},g}}$(z)")
                    ax.grid(True)

                    # Plot each group flux
                    for g, group in enumerate(groups):
                        ax.plot(z, phi[t][g], label=f'Group {group}')
                    ax.legend()
                    fig.tight_layout()
        plt.show()

    def plot_precursors_species(self, material=0, precursors=None,
                                times=None, plot_type='precursors'):
        """
        Plot the precursor species at various times.

        :param material: The material ID
        :type material: int
        :param precursors: The precursor species index on the material.
            Default is None, which corresponds to all precursors.
        :type precursors: int or list of int, optional
        :param times: The times to get the precursors at.
            Default is None, which returns all available times.
        :type times: float or list of float, optional
        :param plot_type: The plotting method. his parameter determines
            whether the each subplot contains all the specified species at
            a fixed time or vice verse. The parameter specifies what is
            plotted on an individual subplot. Default is `precursors`.
        :type plot_type: str {`precursors`, `time`}, optional
        """
        assert material < self.n_materials

        # Get precursors to plot
        n_precursors = self.precursors_per_material[material]
        if precursors is None:
            precursors = list(range(n_precursors))
        elif isinstance(precursors, int):
            precursors = [precursors]
        assert all([j < n_precursors for j in precursors])

        # Get times to plot
        if times is None:
            times = [self.times[0], self.times[-1]]
        elif isinstance(times, float):
            times = [times]
        t0, tf = self.times[0], self.times[-1]
        assert all([t0 <= t <= tf for t in times])

        # Plot 1D profiles
        if self.dimension == 1:

            # Get the coordinates
            z = [p.z for p in self.nodes]

            if plot_type == 'precursors':

                c_j = self.get_precursors(material, precursors, times)

                # Create a figure for each time
                for t, time in enumerate(times):
                    fig: Figure = plt.figure(t)
                    ax: Axes = fig.add_subplot(1, 1, 1)
                    ax.set_title(fr'Time = {time:.2e} $\mu$s')
                    ax.set_xlabel(f'z (cm)')
                    ax.set_ylabel(fr"$C^{{{material}}}_{{j}}$(z)")
                    ax.grid(True)

                    # Plot each group flux
                    for j, precursor in enumerate(precursors):
                        ax.plot(z, c_j[t][j], label=f'Species {precursor}')
                    ax.legend()
                    fig.tight_layout()
        plt.show()

