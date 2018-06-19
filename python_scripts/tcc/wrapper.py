"""Python wrapper for running the TCC.
The main class is TCCWrapper, but there are a number of support classes
to handle the parameters needed to interface with the TCC.
"""

import os
import sys
import tempfile
import shutil
import pandas
import subprocess
import platform
from glob import glob

from python_scripts.file_readers import xyz
from python_scripts.tcc import structures


class BoxType:
    """Enumerations for the simulation box options."""
    NVT = 1
    NPT = 2
    triclinic = 3


class BondType:
    """Enumerations for the bond detection algorithm."""
    simple = 0
    voronoi = 1


class TCCDefaults:
    @staticmethod
    def parameters(sigma=1.):
        """Default input parameters in the TCC input INI file.
        Refer to the example INI files for more detailed descriptions of
        the various parameters.

        Args:
            sigma: diameter of particles for rescaling the bond cutoffs
        Returns:
            dictionary of different sections of input script (as dictionaries of key,value pairs).
        """
        parameters = dict()
        parameters['Box'] = TCCDefaults.box()
        parameters['Run'] = TCCDefaults.run()
        parameters['Simulation'] = TCCDefaults.simulation(sigma)
        parameters['Output'] = TCCDefaults.output()
        return parameters

    @staticmethod
    def box():
        """Defaults passed to the [Box] section of the TCC input file.

        Returns:
            dictionary of the default parameters for this input script
        """

        section = dict()

        section['box_type'] = BoxType.NVT
        section['box_name'] = 'box.txt'

        return section

    @staticmethod
    def run():
        """Defaults passed to the [Run] section of the TCC input file.

        Returns:
            dictionary of the default parameters for this input script
        """

        section = dict()

        section['xyzfilename'] = 'run.xyz'
        section['frames'] = 1
        section['sample_frequency'] = 1

        return section

    @staticmethod
    def simulation(sigma=1.):
        """Defaults passed to the [Simulation] section of the TCC input file.

        Args:
            sigma: diameter of particles for rescaling the bond cutoffs
        Returns:
            dictionary of the default parameters for this input script
        """

        section = dict()

        # Distance cutoff for bond detection between different species
        section['rcutAA'] = 1.4*sigma
        section['rcutAB'] = 1.4*sigma
        section['rcutBB'] = 1.4*sigma
        section['min_cutAA'] = 0.

        section['bond_type'] = BondType.voronoi  # Bond detection algorithm
        section['PBCs'] = True                   # Periodic boundary conditions
        section['voronoi_parameter'] = 0.82      # Mysterious Fc parameter (citation?)
        section['num_bonds'] = 30                # Maximum number of bonds for one particle
        section['cell_list'] = True              # Optimise bond detection using cell lists

        return section

    @staticmethod
    def output():
        """Defaults passed to the [Output] section of the TCC input file.

        Returns:
            dictionary of the default parameters for this input script
        """

        section = dict()

        section['bonds'] = False
        section['clusts'] = False
        section['raw'] = False
        section['do_XYZ'] = False
        section['11a'] = False
        section['13a'] = False
        section['pop_per_frame'] = False

        return section


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
        self.input_parameters = TCCDefaults.parameters()

    def __del__(self):
        """Upon deletion we can remove the temporary working folder
        to free up disk space."""
        shutil.rmtree(self.working_directory)

    def serialise_input_parameters(self, out=sys.stdout):
        """Serialise the parameters in INI format.
        It's not strictly necessary to output these with section headings,
        but it makes the resulting input script more legible and thus
        easier for debugging.

        Args:
            out: file stream or path to write to
        """
        if type(out) is str:
            out = open(out, 'w')

        for section_heading, section_values in self.input_parameters.items():
            out.write('[%s]\n' % section_heading)

            # Write the key-value pairs
            for key, value in section_values.items():
                if type(value) is bool:
                    value = int(value)
                if type(value) is str:
                    out.write('%s=%s\n' % (key, value))
                else:
                    out.write('%s=%r\n' % (key, value))

            out.write('\n')

    def serialise_box(self, box, out=sys.stdout):
        """Serialise the box box size in the TCC format.

        Args:
            out: file stream or path to write to
        """
        if type(out) is str:
            out = open(out, 'w')
        out.write('#iter Lx Ly Lz\n')
        for frame, lengths in enumerate(box):
            out.write('%d %s\n' % (frame, ' '.join(map(str, lengths))))

    def parse_static_clusters(self):
        """Retrive the static cluster information after running the TCC.
        Returns:
            pandas table giving the static cluster information
        """
        summary_file = glob('%s/*.static_clust' % self.working_directory)[0]
        table = pandas.read_table(summary_file, index_col='Cluster type', skiprows=1, nrows=len(structures.cluster_list))
        table.fillna(0., inplace=True)
        return table

    @staticmethod
    def get_tcc_executable_path():
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

    def run(self, box, x, atoms='A', silent=True, **output):
        """Run the TCC

        Args:
            box: box size for boundary conditions
            x: coordinates of atoms
            atoms: species of atoms individually (if given container) or collectively. This must be either length 1 (if specifying species of all atoms) or the same length as the number of particles.
            silent: if set TCC executable console output will be suppressed
            output: dictionary of flags setting output options (see example INI files for possible options, or defaults)
        Returns:
            pandas table giving the static cluster information
        """

        tcc_path = self.get_tcc_executable_path()

        # Create the INI file.
        self.input_parameters['Output'] = TCCDefaults.output()
        for key, value in output.items():
            self.input_parameters['Output'][key] = value
        self.serialise_input_parameters('%s/inputparameters.ini' % self.working_directory)

        # Create the box and configuration files.
        self.serialise_box([box], '%s/box.txt' % self.working_directory)
        xyz.write(x, '%s/run.xyz' % self.working_directory, species=atoms)

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
