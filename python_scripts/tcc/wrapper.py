"""Python wrapper for running the TCC."""

import os
import tempfile
import shutil
import pandas
import subprocess
import platform
from glob import glob

from python_scripts.file_readers import xyz
from python_scripts.tcc import structures


class TCCWrapper:
    """Python interface to the TCC executable.

    The TCC accepts input parameters through a file IO system, so
    this wrapper acts as an intermediate layer to streamline the
    process of running the TCC from within python. All file operations
    are hidden from the user to create a more pythonic interface to
    the TCC.

    On destruction of a wrapper object all of the file input/outputs
    are destroyed, so the user must be careful to extract all of the
    data they need (for e.g. postprocessing) and store this somewhere.
    """

    def __init__(self):
        """On initialisation we have to create a temporary directory
        where file operations will be performed behind the scenes."""
        self.working_directory = tempfile.mkdtemp(prefix='TCC_')
        self.input_parameters = dict()
        self.input_parameters['Box'] = dict()
        self.input_parameters['Run'] = dict()
        self.input_parameters['Simulation'] = dict()
        self.input_parameters['Output'] = dict()

    def __del__(self):
        """Upon deletion we can remove the temporary working folder
        to free up disk space."""
        shutil.rmtree(self.working_directory)

    def run(self, box, particle_coordinates, particle_types='A', silent=True):
        """Run the TCC

        Args:
            box: box size for boundary conditions, list of [len_x, len_y, len_z]
            particle_coordinates: coordinates of atoms
            particle_types: species of atoms individually (if given container) or collectively. This must be either length 1 (if specifying species of all atoms) or the same length as the number of particles.
            silent: if set TCC executable console output will be suppressed
        Returns:
            pandas table giving the static cluster information
        """

        tcc_path = self.get_tcc_executable_path()

        # Create the INI file.
        self.serialise_input_parameters('{}/inputparameters.ini'.format(self.working_directory))

        # Create the box and configuration files.
        self.write_box_file([box], '{}/box.txt'.format(self.working_directory))
        xyz.write(particle_coordinates, '{}/run.xyz'.format(self.working_directory), species=particle_types)

        # Run the TCC executable.
        if silent:
            subprocess_result = subprocess.run([tcc_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=self.working_directory)
        else:
            subprocess_result = subprocess.run([tcc_path], cwd=self.working_directory)

        if subprocess_result.returncode == 0:
            return self.parse_static_clusters()
        else:
            self.__del__()
            print("Error: TCC was not able to run.")
            raise Exception

    @staticmethod
    def get_tcc_executable_path():
        """Find the full path for the tcc executable. It is expected to be in ../../bin relative this this script
        Returns:
            String describing full path of TCC executable.
        """
        bin_directory = os.path.abspath(os.path.dirname(__file__) + '/../../bin/')
        if platform.system() == "Windows":
            tcc_exe = bin_directory + "\\tcc.exe"
        else:
            tcc_exe = bin_directory + "/tcc"
        if os.path.exists(tcc_exe):
            return tcc_exe
        else:
            print("TCC executable not found in bin directory.")
            raise FileNotFoundError

    def serialise_input_parameters(self, output_filename):
        """Serialise the parameters in INI format.

        Args:
            output_filename: file path to write output parameters to
        """
        with open(output_filename, 'w') as output_file:
            for section_heading, section_values in self.input_parameters.items():

                output_file.write('[{}]\n'.format(section_heading))
                # Write the key-value pairs
                for key, value in section_values.items():
                    if type(value) is bool:
                        value = int(value)
                    if type(value) is str:
                        output_file.write('{}={}\n'.format(key, value))
                    else:
                        output_file.write('{}={}\n'.format(key, value))

                output_file.write('\n')

    @staticmethod
    def write_box_file(box, box_filename):
        """Serialise the box box size in the TCC format.

        Args:
            box: Box dimensions are given as a list of the format [len_x, len_y, len_z].
                If NVT this should be a single set of box dimensions.
                If NPT this sould be a list of sets of box dimensions, one for each timestep.
            box_filename: file path to write box to
        """
        with open(box_filename, 'w') as output_file:
            output_file.write('#iter Lx Ly Lz\n')
            for frame, lengths in enumerate(box):
                output_file.write('{} {}\n'.format(frame, ' '.join(map(str, lengths))))

    def parse_static_clusters(self):
        """Retrive the static cluster information after running the TCC.
        Returns:
            pandas table giving the static cluster information
        """
        summary_file = glob('%s/*.static_clust' % self.working_directory)[0]
        table = pandas.read_table(summary_file, index_col='Cluster type', skiprows=1, nrows=len(structures.cluster_list))
        table.fillna(0., inplace=True)
        return table
