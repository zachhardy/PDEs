import os
import struct
import bisect

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.pyplot import Figure
from matplotlib.pyplot import Axes

from python_utils import Point
from python_utils import SimulationReader
from .steadystate_reader import SteadyStateNeutronicsReader

plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 12


class KEigenvalueNeutronicsReader(SteadyStateNeutronicsReader):

    def __init__(self, path):

        super().__init__(path)

        self.k_eff = 0.0

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

        # Open the snapshot file
        with open(self.path, mode='rb') as file:
            file.read(799)  # skip over header

            self.k_eff = read_double(file)

            n_data_blocks = read_unsigned_int(file)

            self.n_cells = read_uint64_t(file)
            print(self.n_cells)

            self.n_nodes = read_uint64_t(file)
            print(self.n_nodes)

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
                p = Point(*[read_double(file) for _ in range(3)])
                self.centroids = np.hstack((self.centroids, p))

                # Parse the nodes
                for i in range(nodes_per_cell):
                    p = Point(*[read_double(file) for _ in range(3)])
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
