"""Python interface to the TCC executable."""

import os
import tempfile
import shutil
import pandas
import subprocess
import platform
from glob import glob

from tcc_python_scripts.file_readers import xyz
from tcc_python_scripts.tcc import structures


class TCCWrapper:
    """Python interface to the TCC executable. Runs the TCC on a
    single configuration.

    The TCC accepts input parameters through a file IO system, so
    this wrapper acts as an intermediate layer to streamline the
    process of running the TCC from within python. All file operations
    are hidden from the user to create a more pythonic interface to
    the TCC.

    On destruction of a wrapper object all of the file input/outputs
    are destroyed, so the user must be careful to extract all of the
    data they need (for e.g. postprocessing) and store this somewhere.

    Attributes:
        input_parameters['Box']: TCC box paramaters used for TCC run
        input_parameters['Run']: TCC run paramaters used for TCC run
        input_parameters['Simulation']: TCC simulation paramaters used for TCC run
        input_parameters['Output']: TCC output paramaters used for TCC run
    """

    def __init__(self):
        """On initialisation we have to create a temporary directory
        where file operations will be performed behind the scenes."""
        self.working_directory = None
        self.input_parameters = dict()
        self.input_parameters['Box'] = dict()
        self.input_parameters['Run'] = dict()
        self.input_parameters['Simulation'] = dict()
        self.input_parameters['Output'] = dict()

    def __del__(self):
        """Upon deletion we can remove the temporary working folder
        to free up disk space."""
        shutil.rmtree(self.working_directory)

    def run(self, box, particle_coordinates, output_directory=None, particle_types='A', silent=True):
        """Invoke the TCC using the provided coordinates and parameters.

        Args:
            box: box size for boundary conditions, list of [len_x, len_y, len_z]
            particle_coordinates: a list of lists, one list for each frame containing coordinates of atoms
            output_directory: If you want to save the output of the TCC specify a directory to store the output
            particle_types: species of atoms individually (if given container) or collectively. This must be either length 1 (if specifying species of all atoms) or the same length as the number of particles.
            silent: if set TCC executable console output will be suppressed
        Returns:
            pandas table containing the static cluster information
        """

        tcc_path = self.get_tcc_executable_path()
        self._set_up_working_directory(output_directory)

        # Create the INI file.
        self._serialise_input_parameters('{}/inputparameters.ini'.format(self.working_directory))

        # Create the box and configuration files.
        self._write_box_file(box, '{}/box.txt'.format(self.working_directory))
        xyz.write('{}/sample.xyz'.format(self.working_directory), particle_coordinates, species=particle_types)

        # Run the TCC executable.
        if silent:
            subprocess_result = subprocess.run([tcc_path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=self.working_directory)
        else:
            subprocess_result = subprocess.run([tcc_path], cwd=self.working_directory)

        if subprocess_result.returncode == 0:
            return self._parse_static_clusters()
        else:
            self.__del__()
            print("Error: TCC was not able to run.")
            raise Exception

    def _set_up_working_directory(self, output_directory):
        """
        Work out where to run the TCC. Create a directory if necessary.

        Args:
            output_directory (str):  If None run the TCC in a temporary directory.
                If a path then run the TCC there
        """
        if output_directory is None:
            self.working_directory = tempfile.mkdtemp(prefix='TCC_')
        else:
            if not os.path.exists(output_directory):
                try:
                    os.makedirs(output_directory)
                except os.error:
                    print("Unable to create output directory: {}. "
                          "Check location is not write protected.".format(output_directory))
                    raise os.error
            self.working_directory = output_directory

    @staticmethod
    def get_tcc_executable_path():
        """Find the full path for the tcc executable. It is expected to be in ../../bin relative this this script
        
        Returns:
            (str) Full path of TCC executable.
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

    def _serialise_input_parameters(self, output_filename):
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
    def _write_box_file(box, box_filename):
        """Serialise the box box size in the TCC format.

        Args:
            box: Box dimensions are given as a list of the format [len_x, len_y, len_z].
            box_filename: file path to write box to
        """
        with open(box_filename, 'w') as output_file:
            output_file.write('#iter Lx Ly Lz\n')
            output_file.write('1\t')
            for dimension in box:
                output_file.write('{}\t'.format(dimension))

    def _parse_static_clusters(self):
        """Retrive the static cluster information after running the TCC.
        Returns:
            Pandas dataframe containing the static cluster information
        """
        summary_file = glob('%s/*.static_clust' % self.working_directory)[0]
        table = pandas.read_table(summary_file, index_col='Cluster type', skiprows=1, nrows=len(structures.cluster_list))
        table.fillna(0., inplace=True)
        return table
