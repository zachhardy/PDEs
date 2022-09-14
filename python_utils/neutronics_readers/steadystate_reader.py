import os
import struct

import numpy as np
import matplotlib.pyplot as plt

from .. import CartesianVector
from .. import SimulationReader

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 12


class SteadyStateNeutronicsReader(SimulationReader):

    def __init__(self, path):
        if not os.path.isfile(path):
            raise FileNotFoundError
        super().__init__(path)

        # Spatial Information
        self.dimension = 0
        self.n_cells = 0
        self.n_nodes = 0
        self.nodes_per_cell = 0
        self.n_materials = 0

        self.nodes = np.empty(0, dtype=CartesianVector)
        self.centroids = np.empty(0, dtype=CartesianVector)
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

    def read(self):
        """
        Read a neutronics data files to populate simulation data.
        """
        self.clear()

        def read_uint64_t(f):
            return struct.unpack('Q', f.read(8))[0]

        def read_unsigned_int(f):
            return struct.unpack('I', f.read(4))[0]

        def read_double(f):
            return struct.unpack('d', f.read(8))[0]

        # Open the snapshot file
        with open(self.path, mode='rb') as file:
            file.read(999)  # skip over header

            n_data_blocks = read_unsigned_int(file)

            self.n_cells = read_uint64_t(file)
            self.n_nodes = read_uint64_t(file)

            self.n_moments = read_unsigned_int(file)
            self.n_groups = read_unsigned_int(file)
            self.n_precursors = read_unsigned_int(file)
            self.max_precursors = read_unsigned_int(file)

            self.init_storage()

            # Parse discretization information
            node_count = 0
            for c in range(self.n_cells):
                cell_id = read_uint64_t(file)
                material_id = read_uint64_t(file)
                nodes_per_cell = read_unsigned_int(file)

                # If on the first snapshot, save material IDs
                self.material_ids.append(material_id)
                if c == 0:
                    self.nodes_per_cell = nodes_per_cell

                # Parse the centroid
                p = CartesianVector(*[read_double(file) for _ in range(3)])
                self.centroids = np.hstack((self.centroids, p))

                # Parse the nodes
                for i in range(nodes_per_cell):
                    p = CartesianVector(*[read_double(file) for _ in range(3)])
                    self.nodes = np.hstack((self.nodes, p))

            # Set the number of materials
            self.n_materials = len(np.unique(self.material_ids))

            # Check centroids and nodes
            assert len(self.centroids) == self.n_cells
            assert len(self.nodes) == self.n_nodes

            self.centroids = np.array(self.centroids, CartesianVector)
            self.nodes = np.array(self.nodes, CartesianVector)

            # Parse flux moments
            assert read_unsigned_int(file) == 0
            phi = self.flux_moments

            n_records = read_uint64_t(file)
            for i in range(n_records):
                dof = self.map_phi_dof(
                    read_uint64_t(file),      # cell
                    read_unsigned_int(file),  # node
                    read_unsigned_int(file),  # moment
                    read_unsigned_int(file))  # group

                phi[dof] = read_double(file)  # value

            # Parse precursors
            assert read_unsigned_int(file) == 1
            precursors = self.precursors
            ppm = [set() for _ in range(self.n_materials)]

            n_records = read_uint64_t(file)
            for i in range(n_records):
                cell = read_uint64_t(file)
                material_id = read_uint64_t(file)
                precursor = read_unsigned_int(file)

                ppm[material_id].add(precursor)

                dof = self.map_precursor_dof(cell, precursor)
                precursors[dof] = read_double(file)

                # print(cell, material_id, precursor, precursors[dof])

            # Define the number of precursors per material
            for m in range(self.n_materials):
                self.precursors_per_material.append(len(ppm[m]))

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
        n_c, n_n = self.n_cells, self.n_nodes
        n_m, n_g, n_p = self.n_moments, self.n_groups, self.max_precursors

        self.flux_moments = np.zeros(n_m * n_g * n_n)
        self.precursors = np.zeros(n_p * n_c)
        self.power_density = np.zeros(n_c)

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

    def get_flux_moment(self, moment, groups=None):
        """
        Return a multi-group flux moment vector at all times.

        :param moment: The moment index
        :type moment: int
        :param groups: The group index.
            Default is None, which returns all groups.
        :type groups: int or list of int, optional
        :rtype: numpy.ndarray:
        """
        assert moment < self.n_moments

        if groups is None:
            groups = list(range(self.n_groups))
        elif isinstance(groups, int):
            groups = [groups]
        assert all([g < self.n_groups for g in groups])

        # Storage for the multigroup flux moment
        out = np.empty((len(groups), self.n_nodes))

        # Loop over cells, nodes, groups, and times
        for cell in range(self.n_cells):
            for node in range(self.nodes_per_cell):

                # The global node index
                pos = cell * self.nodes_per_cell + node

                for g, group in enumerate(groups):
                    dof = self.map_phi_dof(cell, node, moment, group)
                    out[g][pos] = self.flux_moments[dof]
        return out

    def get_precursors(self, material, precursors=None):
        """
        Return all precursors that live on the specified @p material.

        :param material: The material ID
        :type material: int
        :param precursors: The precursor species index on the material.
            Default is None, which corresponds to all precursors.
        :type precursors: int or list of int, optional
        :rtype: numpy.ndarray
        """
        assert material < self.n_materials

        if precursors is None:
            precursors = list(range(self.max_precursors))
        elif isinstance(precursors, int):
            precursors = [precursors]
        n_precursors = self.precursors_per_material[material]
        assert all([j < n_precursors for j in precursors])

        # Storage for the precursor species
        out = np.zeros((len(precursors), self.n_cells))

        # Loop over cells, precursors, and times
        for cell in range(self.n_cells):

            # If not the correct material, skip
            if self.material_ids[cell] != material:
                continue

            for j, precursor in enumerate(precursors):
                dof = self.map_precursor_dof(cell, precursor)
                out[j][cell] = self.precursors[dof]
        return out

    def plot_flux_moment(self, moment=0, groups=None):
        """
        Plot a multi-group flux moment at various times.

        :param moment: The moment index. Default is 0.
        :type moment: int, optional
        :param groups: The group indices. Default is all groups.
        :type: groups: int or list of int, optional
        """
        # Get the groups to plot
        if groups is None:
            groups = list(range(self.n_groups))
        elif isinstance(groups, int):
            groups = [groups]
        assert all([g < self.n_groups for g in groups])

        phi = self.get_flux_moment(moment, groups)

        # Plot 1D profiles
        if self.dimension == 1:
            z = [p.z for p in self.nodes]

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel(f'z (cm)')
            ax.set_ylabel(fr"$\phi_{{{moment},g}}$(z)")
            ax.grid(True)

            # Plot each group flux
            for g, group in enumerate(groups):
                ax.plot(z, phi[g], label=f'Group {group}')
            ax.legend()
            fig.tight_layout()

        # Plot 2D profiles
        elif self.dimension == 2:
            x = np.array([p.x for p in self.nodes])
            y = np.array([p.y for p in self.nodes])
            X, Y = np.meshgrid(np.unique(x), np.unique(y))

            for g, group in enumerate(groups):
                fig = plt.figure()
                fig.suptitle(f"Moment {moment}, Group {group} Flux Moment")
                ax = fig.add_subplot(1, 1, 1)
                ax.set_xlabel("X (cm)")
                ax.set_ylabel("Y (cm)")

                phi_fmtd = phi[g].reshape(X.shape)
                im = ax.pcolor(X, Y, phi_fmtd, cmap="jet", shading="auto",
                               vmin=0.0, vmax=phi_fmtd.max())
                fig.colorbar(im)
                fig.tight_layout()
        plt.show()

    def plot_precursors_species(self, material=0, precursors=None):
        """
        Plot the precursor species at various times.

        :param material: The material ID
        :type material: int
        :param precursors: The precursor species index on the material.
            Default is None, which corresponds to all precursors.
        :type precursors: int or list of int, optional
        """
        assert material < self.n_materials

        # Get precursors to plot
        n_precursors = self.precursors_per_material[material]
        if precursors is None:
            precursors = list(range(n_precursors))
        elif isinstance(precursors, int):
            precursors = [precursors]
        assert all([j < n_precursors for j in precursors])

        # Plot 1D profiles
        if self.dimension == 1:

            # Get the coordinates
            z = [p.z for p in self.nodes]

            c_j = self.get_precursors(material, precursors)

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.set_xlabel(f"z (cm)")
            ax.set_ylabel(fr"$C^{{{material}}}_{{j}}$(z)")
            ax.grid(True)

            # Plot each group flux
            for j, precursor in enumerate(precursors):
                ax.plot(z, c_j[j], label=f'Species {precursor}')
            ax.legend()
            fig.tight_layout()
        plt.show()

    def plot_material_ids(self):
        plt.figure()

        if self.dimension == 1:
            z = [centroid.z for centroid in self.centroids]

            plt.xlabel("z (cm)")
            plt.plot(z, self.material_ids, 'ob')

        elif self.dimension == 2:
            x = [centroid.x for centroid in self.centroids]
            y = [centroid.y for centroid in self.centroids]
            X, Y = np.meshgrid(np.unique(x), np.unique(y))
            matids = np.array(self.material_ids).reshape(X.shape)
            plt.pcolormesh(X, Y, matids, cmap="jet", shading="auto")
            plt.colorbar()

        plt.tight_layout()
        plt.show()

