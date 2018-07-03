import sys
import os
sys.path.append(os.path.abspath("../../python_scripts/"))

from file_readers import xyz

particle_coordinates = xyz.read("../../test/integration_tests/basic_voronoi/sample.xyz", num_frames=5)

