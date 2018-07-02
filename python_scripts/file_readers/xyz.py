#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to XYZ (.xyz) file formats.

The module defines:
  - XYZSnapshot: the class the defining the file interface to this file format
  - read: shorthand for XYZSnapshot.read_trajectory
  - write: create a snapshot from coordinates to write to disk
"""

import sys
import io
import numpy
import pandas
from python_scripts.file_readers.snapshot import stream_safe_open, NoSnapshotError, SnapshotIncompleteError, Snapshot


class XYZSnapshot(Snapshot):
    """Snapshot of a system of particles in XYZ (.xyz) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def read(self, path_or_file):
        """Read a snapshot from a file, overwriting any existing data.

        Args:
            path_or_file: file stream or path to read snapshot from
        Raises:
            NoSnapshotError: if could not read from file
        """
        with stream_safe_open(path_or_file) as input_file:
            # Read the header - check for EOF.
            number_of_particles = self.process_number_of_particles(input_file)

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
    def process_number_of_particles(f):
        line = f.readline()
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
    """Read one or more snapshots from the xyz file."""
    return XYZSnapshot.read_trajectory(file_name, num_frames)


def write(x, out=sys.stdout, species=None):
    """ Write a single configuration to the disk."""
    snapshot = XYZSnapshot(x, species=species)
    snapshot.write(out)
