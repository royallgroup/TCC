#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to XYZ (.xyz) file formats. The main class is XYZSnapshot, but some additional functions are defined to provide a simplified interface to this class.

The module defines:
  - XYZSnapshot: the class the defining the file interface to this file format
  - read: shorthand for XYZSnapshot.read_single
  - read_trajectory: shorthand for XYZSnapshot.read_trajectory
  - write: create a snapshot from coordinates to write to disk
"""

import sys, io, numpy, pandas
from .snapshot import stream_safe_open, NoSnapshotError, Snapshot

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
        with stream_safe_open(path_or_file) as f:
            # Read the header - check for EOF.
            line = f.readline()
            if not line: raise NoSnapshotError

            # Process the rest of the header rows.
            number_of_atoms = int(line)
            if not number_of_atoms > 0: raise NoSnapshotError
            comment = f.readline()

            # Use pandas to read the main table.
            c = io.StringIO()
            for i in range(number_of_atoms): c.write(f.readline())
            c.seek(0)
            table = pandas.read_table(c, sep='\s+', names=('atom','x','y','z'), nrows=number_of_atoms)

            self.x = table[['x','y','z']].values.copy('c').astype(numpy.longdouble)
            self.species = table['atom'].tolist()
            self.time = self.box = None

    def __str__(self):
        """String representation of the snapshot in XYZ (.xyz) format"""
        f = io.StringIO()

        # Header states number of particles (we have ignored comment line)
        f.write('%d\n' % self.n)

        # Single component system
        if type(self.species) is str:
            for i in range(self.n):
                f.write('\n')
                f.write('%s ' % self.species)
                f.write(' '.join(map(str, self.x[i,:])))
        # Handle particle species separately
        else:
            for i in range(self.n):
                f.write('\n')
                f.write('%s ' % self.species[i])
                f.write(' '.join(map(str, self.x[i,:])))

        return f.getvalue()

def read(*args, **kwargs):
    """Read a single snapshot from the disk."""
    return XYZSnapshot.read_single(*args, **kwargs)
def read_trajectory(*args, **kwargs):
    """Read a trajectory (i.e. multiple snapshots) from the disk."""
    return XYZSnapshot.read_trajectory(*args, **kwargs)

def write(x, out=sys.stdout, species=None):
    """ Wrtie a single configuration to the disk."""
    snapshot = XYZSnapshot(x, species=species)
    snapshot.write(out)
