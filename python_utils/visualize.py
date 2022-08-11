import os
import sys

from readers import TransientNeutronicsReader

import matplotlib.pyplot as plt

# Get path to stuff to read
if len(sys.argv) > 1:
    path = sys.argv[1]
else:
    path = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(path, "../Problems/LRA/outputs")
path = os.path.abspath(path)
assert os.path.isdir(path)

# Define defaults
logscale = False

# Parse auxiliary arguments
for arg in sys.argv[1:]:
    if "=" in arg:
        if "log" in arg:
            arg = arg.split("=")[1]
            assert arg in ["0", "1"]
            logscale = bool(int(arg))

sim = TransientNeutronicsReader(path)
sim.read()

sim.plot_power('total', logscale)
plt.show()
