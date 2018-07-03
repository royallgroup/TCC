import sys
import os
sys.path.append(os.path.abspath("../../python_scripts/"))

from tcc import wrapper
from file_readers import xyz

# Open a TCCWrapper object - this holds information about the simulation we want to run
TCC_setup = wrapper.TCCWrapper()

# Get the box size. This can be read from a file or input manually
box = [10, 10, 10]

# Get the coordinates. The file_readers scripts are a good way to read in coordinates from a file.
particle_coordinates = list(xyz.read("../../test/integration_tests/basic_voronoi/sample.xyz", num_frames=1))[0].particle_coordinates

# If a results directory is not specified the TCC will run in a temporary folder returning only
# the average cluster results for the whole run.
TCC_setup.input_parameters['Run']['Frames'] = 1
results = TCC_setup.run(box, particle_coordinates, silent=False)

print(results['Number of clusters'])
print(results['Mean Pop Per Frame'])
