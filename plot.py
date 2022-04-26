import os
import struct
import numpy as np
import matplotlib.pyplot as plt


def read(filename):
    with open(filename, 'rb') as file:
        file.read(499)
        n_cells = struct.unpack('Q', file.read(8))[0]
        n_nodes = struct.unpack('Q', file.read(8))[0]
        n_groups = struct.unpack('Q', file.read(8))[0]

        grid = []
        solution = [[] for _ in range(n_groups)]

        node_count = 0
        for c in range(n_cells):
            cell_id = struct.unpack('Q', file.read(8))[0]
            material_id = struct.unpack('Q', file.read(8))[0]
            nodes_per_cell = struct.unpack('Q', file.read(8))[0]

            for n in range(nodes_per_cell):
                node = []
                for d in range(3):
                    coord = struct.unpack('d', file.read(8))[0]
                    node.append(coord)
                grid.append(node)

        n_records = n_nodes * n_groups
        for r in range(n_records):
            cell_id = struct.unpack('Q', file.read(8))[0]
            node = struct.unpack('Q', file.read(8))[0]
            group = struct.unpack('Q', file.read(8))[0]
            value = struct.unpack('d', file.read(8))[0]

            solution[group].append(value)
    return np.array(grid), np.array(solution)


if __name__ == "__main__":
    nodes, phi = read("solution.data")
    z = nodes[:, 2]

    plt.figure()
    for g in range(phi.shape[0]):
        plt.plot(z, phi[g], label=f'Group {g}')
    plt.legend()
    plt.show()
