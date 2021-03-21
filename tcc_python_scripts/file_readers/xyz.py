""" Module for reading and writing snapshots from and to XYZ (.xyz) file formats."""

import io
import numpy
import pandas
from tcc_python_scripts.file_readers.snapshot import stream_safe_open, NoSnapshotError, SnapshotIncompleteError, Snapshot


class XYZSnapshot(Snapshot):
    """Snapshot of a system of particles in XYZ (.xyz) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def _read(self, path_or_file):
        """ Read a single XYZ snapshot from a file.

        Overwrites any exisiting data in the Snaphsot object.

        Raises:
            NoSnapshotError if file could not be read.
        Args:
            path_or_file: file stream or path to read snapshot from
        """
        with stream_safe_open(path_or_file) as input_file:
            # Read the header - check for EOF.
            number_of_particles = self._process_number_of_particles(input_file)

            # Read and igore the comment line
            input_file.readline()

            # Use pandas to read the main table.
            string_buffer = io.StringIO()
            for line_number in range(number_of_particles):
                line = (input_file.readline())
                if len(line.split()) != 4:
                    raise SnapshotIncompleteError("Error reading XYZ file on line number {}".format(line_number))
                string_buffer.write(line)
            string_buffer.seek(0)
            table = pandas.read_table(string_buffer, sep='\s+', names=('atom', 'x', 'y', 'z'), nrows=number_of_particles)
            if table.shape[0] != number_of_particles:
                raise SnapshotIncompleteError
            self.particle_coordinates = table[['x', 'y', 'z']].values.copy('c').astype(numpy.longdouble)
            self.species = table['atom'].tolist()
            self.time = self.box = None

    @staticmethod
    def _process_number_of_particles(file):
        """ Sanitises the number of particles read from an XYZ file.

        Args:
            an open file handle to the xyz file
        Raises:
            SnapshotIncompleteError if number of particles cannot be interpreted
        Returns:
            an integer number of particles.
        """
        line = file.readline()
        if not line:
            raise NoSnapshotError
        if len(line.split()) > 1:
            raise SnapshotIncompleteError("Can't read number of particles from XYZ file.")
        try:
            number_of_particles = int(line)
        except ValueError:
            raise SnapshotIncompleteError("Can't read number of particles from XYZ file.")
        if not number_of_particles > 0:
            raise NoSnapshotError
        return number_of_particles

    def __str__(self):
        """String representation of the snapshot in XYZ (.xyz) format"""

        f = io.StringIO()
        # Header states number of particles (we have ignored comment line)
        f.write('%d\n' % self.num_particles)

        # Single component system
        if type(self.species) is str:
            for i in range(self.num_particles):
                f.write('\n')
                f.write('%s ' % self.species)
                f.write(' '.join(map(str, self.particle_coordinates[i, :])))
        # Handle particle species separately
        else:
            for i in range(self.num_particles):
                f.write('\n')
                f.write('%s ' % self.species[i])
                f.write(' '.join(map(str, self.particle_coordinates[i, :])))

        return f.getvalue()


def read(file_name):
    """ Returns a generator that reads one or more snapshots from an xyz file.

    Args:
        file_name: Name of XYZ file to read.
    Returns:
         A generator object that generates Snapshots.
    """
    return XYZSnapshot.read_trajectory(file_name)


def write_snapshot(output_filename, snapshot_list):
    """ Write one or more snapshots to the disk.

    Args:
        output_filename: The filename to write the coordinates to.
        snapshot_list: A list of XYZSnapshots.
    """
    for snapshot in snapshot_list:
        snapshot.write(output_filename)


def write(output_filename, particle_coordinates, species=None):
    """ Write a single configuration to the disk.

    Args:
        output_filename: The filename to write the coordinates to.
        particle_coordinates: A list of particle coordinates
        species: A list of particle species. Defaults all particles to 'A' if not provided.
    """
    snapshot = XYZSnapshot(particle_coordinates, species=species)
    snapshot.write(output_filename)


def check_coordinates_type(coordinates):
    """
    Determine if a given particle coordinates are multiple frames or
        just a single frame. Three possbilities were considered:

    1. a numpy array for a single frame, shape (N_paritcle, 3)
    2. a numpy array for multiple frames, shape (N_frame, N_particle, 3)
    3. a list of array for multiple frames, shape (N_frame, N_particle*, 3)
        (*the number of particles in different frames might be different)

    all other possibilities were considered to be invalid
    (TODO: making more data types legal, such as the pandas dataframe)

    Args:
        coordinates (iterable): Particle coordinates with different possible
            shapes / types.

    Return:
        iterable: for #1, return [coordinates];
            for #2 & #3, return the origional coordinates,
    """
    if isinstance(coordinates, numpy.ndarray):
        if (coordinates.ndim == 2) and (coordinates.shape[-1] == 3):  # case 1
            return [coordinates]
        elif (coordinates.ndim == 3) and (coordinates.shape[-1] == 3): # case 2
            return coordinates
        else:
            raise TypeError("Invalid (numpy) coordinates shape of ", coordinates.shape)
    elif isinstance(coordinates, list):
        if numpy.all([isinstance(c, numpy.ndarray) for c in coordinates]):
            if set([c.ndim for c in coordinates]) == set([2]):
                return coordinates
            else:
                raise TypeError("Invalid (list) coordinates containing numpy array with different dimensions")
        else:
            raise TypeError("Invalid (list) coordinates, required: list of numpy arrays")
    else:
        raise TypeError("Invalid coordinates, a list or a numpy array is required")


def write_multiple(output_filename, particle_coordinates, species=None):
    """ Write multiple configurations to the disk.

    Args:
        output_filename: The filename to write the coordinates to.
        particle_coordinates: A list of particle coordinates
        species: A list of particle species. Defaults all particles to 'A' if not provided.
    """
    frames = check_coordinates_type(particle_coordinates)
    for frame in frames:
        snapshot = XYZSnapshot(frame, species=species)
        snapshot.write(output_filename, mode='a')
