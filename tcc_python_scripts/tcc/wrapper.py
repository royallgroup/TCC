"""Python interface to the TCC executable."""

import os
import re
import tempfile
import shutil
import numpy
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
        working_directory: The directory in which the TCC will run
        tcc_executable_directory: The directory containing the TCC executable
        tcc_executable_path: The full path of the TCC executable
        input_parameters['Box']: TCC box paramaters used for TCC run
        input_parameters['Run']: TCC run paramaters used for TCC run
        input_parameters['Simulation']: TCC simulation paramaters used for TCC run
        input_parameters['Output']: TCC output paramaters used for TCC run
        input_parameters['Clusters_to_analyse']: List of clusters to include in the analysis, all are detected if list is empty
    """

    def __init__(self, clusters_to_analyse=None):
        """On initialisation we have to create a temporary directory
        where file operations will be performed behind the scenes.

        Args:
            clusters_to_analyse: list of which structures to perform a structural analysis on.
                 If None (default) then all of them will be used."""
        self.working_directory = None
        self.tcc_executable_directory = None
        if platform.system() == "Windows":
            self.tcc_executable_path = shutil.which('tcc.exe')
        else:
            self.tcc_executable_path = shutil.which('tcc')
        self.input_parameters = dict()
        self.input_parameters['Box'] = dict()
        self.input_parameters['Run'] = dict()
        self.input_parameters['Simulation'] = dict()
        self.input_parameters['Output'] = dict()

        self.clusters_to_analyse = clusters_to_analyse
        self.input_parameters['Simulation']['analyse_all_clusters'] = not clusters_to_analyse

    def __del__(self):
        """Upon deletion we can remove the temporary working folder
        to free up disk space."""
        if self.cleanup:
            shutil.rmtree(self.working_directory)

    def save(self, destination, verbatim=True):
        """Save TCC working directory, with all its output files, to a specified path.

        The working directory is deleted when the wrapper falls out of scope, so this must
        be called to retain all the files.

        Args:
            destination: new directory to save data
            verbatim: if True, then the tree is copied verbatim from the output of the TCC,
                otherwise the xyz/box files are ignored, and the directory is flattened so
                all files within subdirectories appear at the root level.
        """
        if verbatim:
            shutil.copytree(self.working_directory, destination)
        else:
            os.makedirs(destination)
            for name in os.listdir(self.working_directory):
                if name.split('.')[-1] in ['xyz', 'txt']: continue

                full_path = os.path.join(self.working_directory, name)
                if os.path.isdir(full_path):
                    for p in os.listdir(full_path):
                        shutil.copy2(os.path.join(full_path, p), destination)
                else:
                    shutil.copy2(full_path, destination)

    def run(self, box, particle_coordinates, output_directory=None, output_clusters=False, particle_types='A', silent=True):
        """Invoke the TCC using the provided coordinates and parameters.

        Args:
            box: box size for boundary conditions, list of [len_x, len_y, len_z]
            particle_coordinates: a list of lists, one list for each frame containing coordinates of atoms
            output_directory: If you want to save the output of the TCC specify a directory to store the output
            particle_types: species of atoms individually (if given container) or collectively. This must be either length 1
                            (if specifying species of all atoms) or the same length as the number of particles.
            silent: if set TCC executable console output will be suppressed
        Returns:
            pandas table containing the static cluster information
        """

        self._check_tcc_executable_path()
        if output_directory is None:
            output_directory = self.working_directory
            self.cleanup = True
        else:
            self.cleanup = False
        self._set_up_working_directory(output_directory)
        self.input_parameters['Output']['clusts'] = output_clusters

        # Create the box and configuration files.
        self._write_box_file(box, self.working_directory)

        # Create the INI file.
        n_frame = self.input_parameters['Run'].get('frames')
        if not n_frame:
            self.input_parameters['Run']['frames'] = xyz.get_frame_number(
                particle_coordinates
            )
        self._serialise_input_parameters(
            '{}/inputparameters.ini'.format(self.working_directory)
        )

        # Create the box and configuration files.
        self._write_box_file(box, self.working_directory)
        xyz.write_multiple(
            '{}/sample.xyz'.format(self.working_directory),
            particle_coordinates, species=particle_types
        )

        if self.clusters_to_analyse:
            self._write_clusters_to_analyse(self.clusters_to_analyse, self.working_directory)

        # Run the TCC executable.
        if silent:
            subprocess_result = subprocess.run(self.tcc_executable_path, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=self.working_directory)
        else:
            subprocess_result = subprocess.run(self.tcc_executable_path, cwd=self.working_directory)

        if subprocess_result.returncode == 0:
            return self._parse_static_clusters()
        else:
            self.__del__()
            print("Error: TCC was not able to run.")
            raise Exception

    def set_tcc_executable_directory(self, path):
        """A method for setting the directory containing the compiled TCC executable.

        Args:
            path: The directory containing the compiled TCC executable.
        """
        self.tcc_executable_directory = path

    def _set_up_working_directory(self, output_directory):
        """
        Work out where to run the TCC. Create a directory if necessary.

        Args:
            output_directory (str):  If None run the TCC in a temporary directory.
                If a path then run the TCC there.
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

    def _check_tcc_executable_path(self):
        """Check the provided path for the tcc executable is valid.

        Returns:
            If provided executable path is valid, returns full path, else raises FileNotFoundError.
        """
        if self.tcc_executable_path is not None and os.path.exists(self.tcc_executable_path): return

        bin_directory = os.path.abspath(self.tcc_executable_directory)
        if platform.system() == "Windows":
            tcc_exe = bin_directory + "\\tcc.exe"
        else:
            tcc_exe = bin_directory + "/tcc"
        if os.path.exists(tcc_exe):
            self.tcc_executable_path = tcc_exe
        else:
            print("TCC executable not found in provided directory.")
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
    def _write_box_file(box, folder_path):
        """Serialise the box box size in the TCC format.

        Args:
            box: Box dimensions are given as a list of the format [len_x, len_y, len_z].
                 A list of box dimensions can be given for multiple frames.
            folder_path: folder to write box file to
        """
        with open('{}/box.txt'.format(folder_path), 'w') as output_file:
            output_file.write('#iter Lx Ly Lz\n')
            if hasattr(box[0], '__iter__'):
                for i,b in enumerate(box):
                    output_file.write('%d\t' % (i+1))
                    for dimension in b:
                        output_file.write('{}\t'.format(dimension))
                    if i+1 < len(box):
                        output_file.write('\n')
            else:
                output_file.write('1\t')
                for dimension in box:
                    output_file.write('{}\t'.format(dimension))

    @staticmethod
    def _write_clusters_to_analyse(clusters_to_include, folder_path):
        """Write an ini file which specifies which clusters to analyse.

        Args:
            clusters_to_include: A list specifying which cluster types to turn on
            folder_path: folder to write clusters_to_analyse file to
        """
        with open("{}/clusters_to_analyse.ini".format(folder_path), 'w') as output_file:
            output_file.write('[Clusters]\n')
            for cluster in structures.cluster_list:
                if cluster in clusters_to_include:
                    output_file.write("{}\t=\t1\n".format(cluster))
                else:
                    output_file.write("{}\t=\t0\n".format(cluster))

    def _parse_static_clusters(self):
        """Retrieve the static cluster information after running the TCC.
        Returns:
            Pandas dataframe containing the static cluster information
        """
        summary_file = glob('%s/*.static_clust' % self.working_directory)[0]
        table = pandas.read_table(summary_file, index_col='Cluster type', skiprows=1, nrows=len(structures.cluster_list))
        table = table[numpy.isfinite(table['Number of clusters'])]
        return table

    def _parse_cluster_file(self, structure):
        """Retrieve detailed breakdown of clusters after running the TCC.

        Args:
            structure: which structure to analyse
        Returns:
            Numpy array (int) giving the particles in each cluster
        """
        cluster_path = glob('%s/cluster_output/*_%s' % (self.working_directory, structure))[0]

        # Wrap this in a with statement to ignore empty file warnings when no clusters found
        with numpy.warnings.catch_warnings():
            numpy.warnings.simplefilter("ignore")
            clusters = numpy.loadtxt(cluster_path, skiprows=1, dtype=int)

        if len(clusters) == 0: return []
        if len(clusters.shape) == 1: clusters = clusters.reshape(1,-1)

        return clusters

    def _parse_particle_clusters(self, natoms):
        """Determine whether each particle belongs to a certain cluster or not
        Returns:
            Numpy array (bool) saying whether each particle (row) belong to each cluster (column)."""
        nclusters = len(self.active_clusters)
        table = numpy.zeros((natoms, nclusters), dtype=bool)

        for i,structure in enumerate(self.active_clusters):
            found_clusters = self._parse_cluster_file(structure)
            table[found_clusters.reshape(-1), i] = True

        return table

    @property
    def active_clusters(self):
        """Returns: list of clusters active in the analysis."""
        if self.clusters_to_analyse: return self.clusters_to_analyse
        else: return structures.cluster_list

    def get_cluster_dict(self, cluster_names=None):
        """
        Getting the result of particles and the clusters they are in. The
            example output would be like.

        ..code-block::

            {
                'FCC': [frame_1, frame_2, ...]
                '10B': [frame_1, frame_2, ...]
                '8A':  [frame_1, frame_2, ...]
                ...
            }

        and each frame is a numpy array with shape (n_particle, )

        Args:
            cluster_names (list or str): A cluster name, or a list of
                cluster names to be the head of the output table.
                If not provided, all the clusters in the raw_output
                folder will be analysed.

        Return:
            dict: a dict with structure {cluster_name : [frame_1, frame_2, ...]},
                and each frame is a Boolean array with shape (n_particle, 1).

        Example:
            >>> import numpy as np
            >>>
            >>> # 10 frames of 1000 random particles in 3D
            >>> coords = np.random.uniform(0, 10, (10, 1000, 3))
            >>> box = [10] * 3
            >>>
            >>> tcc = TCCWrapper()
            >>> tcc.input_parameters['Output']['Raw'] = True
            >>> _ = tcc.run(box, coords)
            >>>
            >>> cluster_dict = tcc.get_cluster_dict(['13A', '10B'])
            >>>
            >>> cluster_dict.keys()
            dict_keys(['13A', '10B'])
            >>> len(cluster_dict['13A'])
            10
            >>> cluster_dict['13A'][0].shape
            (1000, 1)
            >>> cluster_dict['13A'][0][:3].ravel().tolist()
            [False, False, False]

        """
        raw_out_folder = os.path.join(self.working_directory, 'raw_output')
        if not os.path.isdir(raw_out_folder):
            raise FileNotFoundError(
                "No raw_output data, set [Output][Raw] to True and run again"
            )

        # collect the cluster names if not provided
        if isinstance(cluster_names, type(None)):
            cluster_name_pattern = re.compile(r'sample\.xyz.*raw_(.+)')
            filenames = glob(
                "{folder}/sample.xyz*raw_*".format(folder=raw_out_folder)
            )
            filenames = [os.path.basename(fn) for fn in filenames]
            filenames.sort()
            cluster_names = [
                cluster_name_pattern.match(fn).group(1) for fn in filenames
            ]
        elif isinstance(cluster_names, str):
            cluster_names = [cluster_names]
        elif hasattr(cluster_names, '__len__'):
            if len(cluster_names) == 0:
                raise RuntimeError(
                    "Empty list is illegal for cluster_names"
                )
        else:
            raise TypeError("\n".join((
                "Invalid cluster_names type, the accepted types are",
                "\t1. None",
                "\t2. The string of cluster name",
                "\t3. Common container of cluster names",
                "\t   (list, tuple, numpy array, ...)"
            )))

        # collect the result into a dict, {cluster : [frames_1, frame_2, ...]}
        result = {}
        for cn in cluster_names:
            fn = glob(
                "{folder}/sample.xyz*raw_{cluster_name}".format(
                    folder=raw_out_folder, cluster_name=cn
                )
            )
            if len(fn) == 0:
                raise FileNotFoundError(
                    "Raw output file for cluster {cluster} not found".format(
                        cluster = cn
                    )
                )
            else:
                fn = fn[0]
            is_in_cluster = xyz.get_frames_from_xyz(
                filename=fn, usecols=[0],
                # particle labelled as C and D were considered in the cluster
                convert_func=lambda x : x in ['C', 'D']
            )
            result.update({cn : is_in_cluster})

        return result

    def get_cluster_table(self, cluster_names=None):
        """
        Getting the result of particles and the clusters they are in. The
            output is a list of Boolean tables for many frames. A typical
            table would be like the following.

        ..code-block::

            id, FCC, 13A, 12E, 11F, 10B
            1,    0,   0,   0,   0,   1   # particle 1 is in 10B
            2,    1,   0,   0,   0,   0   # particle 2 is in FCC
            3,    0,   0,   0,   0,   0   # particle 3 is not in any cluster
            ...

        Args:
            cluster_names (list or str): A cluster name, or a list of
                cluster names to be the head of the output table.
                If not provided, all the clusters in the raw_output
                folder will be analysed.

        Return:
            list: a list of pandas table, each table represent one frame


        Example:
            >>> import numpy as np
            >>>
            >>> # 10 frames of 1000 random particles in 3D
            >>> coords = np.random.uniform(0, 10, (10, 1000, 3))
            >>> box = [10] * 3
            >>>
            >>> tcc = TCCWrapper()
            >>> tcc.input_parameters['Output']['Raw'] = True
            >>> _ = tcc.run(box, coords)
            >>>
            >>> cluster_table = tcc.get_cluster_table(["FCC", "13A"])
            >>> len(cluster_table)
            10
            >>> print(cluster_table[0].head())
                 FCC    13A
            0  False  False
            1  False  False
            2  False  False
            3  False  False
            4  False  False
        """
        cluster_dict = self.get_cluster_dict(cluster_names)
        frame_nums = []
        for cn, frames in cluster_dict.items():
            frame_nums.append(len(frames))
        if len(set(frame_nums)) == 1:
            frame_num = frame_nums[0]
        else:
            raise RuntimeError(
                "Inconsistent frame sizes, this is a bug of the TCC wrapper"
            )
        # convert the dictionary to a list of table
        result = []
        for f in range(frame_num):
            data = { key : frame[f].flatten() for key, frame in cluster_dict.items() }
            result.append(
                pandas.DataFrame.from_dict(
                    data=data,
                    orient='columns'
                )
            )
        return result
