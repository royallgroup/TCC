""" Module for reading and writing snapshots from and to XYZ (.xyz) file formats."""

import io
import numpy
import pandas
from tcc_python_scripts.file_readers.snapshot import stream_safe_open, NoSnapshotError, SnapshotIncompleteError, Snapshot


class XYZSnapshot(Snapshot):
    """Snapshot of a system of particles in XYZ (.xyz) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def read(self, path_or_file):
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


def read(file_name, num_frames=1):
    """Read one or more snapshots from an xyz file.

    Args:
        file_name: Name of XYZ file to read.
        num_frames: Number of snapshots to read from the xyz file.
    Returns:
         A list of XYZSnapshot objects.
    """
    return XYZSnapshot.read_trajectory(file_name, num_frames)


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
