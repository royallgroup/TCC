#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to LAMMPS (.atom) file formats. The main class is AtomSnapshot,
but some additional functions are defined to provide a simplified interface to this class.

The module defines:
  - AtomSnapshot: the class the defining the file interface to this file format
  - read: shorthand for AtomSnapshot.read_single
  - read_trajectory: shorthand for AtomSnapshot.read_trajectory
"""

import io
import numpy
import pandas
from snapshot import stream_safe_open, NoSnapshotError, Snapshot


class AtomSnapshot(Snapshot):
    """Snapshot of a system of particles in LAMMPS (.atom) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def read(self, path_or_file):
        """Read a snapshot from a file, overwriting any existing data.

        Args:
            path_or_file: file stream or path to read snapshot from
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format
        """
        with stream_safe_open(path_or_file) as f:
            self.particle_coordinates = numpy.empty((0, 0))
            self.time = self.box = None

            while True:
                item = f.readline().split()
                if not item:
                    raise NoSnapshotError
                assert item[0] == 'ITEM:'

                # Timestep within a trajectory.
                if item[1] == 'TIMESTEP':
                    self.time = int(f.readline())

                # Number of atoms in the header
                elif ' '.join(item[1:4]) == 'NUMBER OF ATOMS':
                    n = int(f.readline())
                    self.particle_coordinates = numpy.empty((n, self.dimensionality))

                # Item containing the bounding box.
                elif ' '.join(item[1:3]) == 'BOX BOUNDS':
                    d = len(item[3:])
                    self.particle_coordinates = numpy.empty((self.num_particles, d))
                    self.box = numpy.zeros((d, 2), dtype=numpy.longdouble)

                    for c in range(d):
                        boundary = f.readline().split()
                        self.box[c][:] = [float(b) for b in boundary]

                # Main table contains the per-atom data. Should come at the end.
                elif item[1] == 'ATOMS':
                    assert self.num_particles > 0
                    assert 1 <= self.dimensionality <= 3
                    assert self.box is not None

                    headings = item[2:]
                    assert 'id' in headings
                    assert 'x' or 'xs' in headings

                    c = io.StringIO()
                    for i in range(n):
                        c.write(f.readline())
                    c.seek(0)
                    table = pandas.read_table(c, index_col=0, sep='\s+', names=headings, nrows=n)
                    #try: table = table.sort_values('id')
                    #except: table = table.sort('id')

                    if 'xs' in headings:
                        cols = ['xs', 'ys', 'zs'][:self.dimensionality]
                        self.particle_coordinates = table[cols].values.copy('c').astype(numpy.longdouble)
                        for c in range(d):
                            self.particle_coordinates[:, c] *= self.box[c]
                    else:
                        cols = ['x', 'y', 'z'][:self.dimensionality]
                        self.particle_coordinates = table[cols].values.copy('c').astype(numpy.longdouble)

                    self.species = numpy.array(table['type'])
                    return

                else:
                    raise RuntimeError('unknown header: %s' % item)

    def __str__(self):
        """String representation of the snapshot in LAMMPS (.atom) format"""
        string_buffer = io.StringIO()
        string_buffer.write('ITEM: TIMESTEP\n%r\n' % self.time)
        string_buffer.write('ITEM: NUMBER OF ATOMS\n%r\n' % self.num_particles)
        string_buffer.write('ITEM: BOX BOUNDS')
        for _ in self.box:
            string_buffer.write(' pp')
        string_buffer.write('\n')
        for dim in self.box:
            if len(dim) == 2:
                string_buffer.write('%.8f %.8f\n' % tuple(dim))
            else:
                string_buffer.write('0 %.8f\n' % dim)
        string_buffer.write('ITEM: ATOMS id type x y z')
        for i, (name, x) in enumerate(zip(self.species, self.particle_coordinates)):
            string_buffer.write('\n')
            string_buffer.write('%r %s' % (i, name))
            for coord in x:
                string_buffer.write(' %.8f' % coord)
        return string_buffer.getvalue()


def read(*args, **kwargs):
    """Read a single snapshot from the disk."""
    return AtomSnapshot.read_single(*args, **kwargs)


def read_trajectory(*args, **kwargs):
    """Read a trajectory (i.e. multiple snapshots) from the disk."""
    return AtomSnapshot.read_trajectory(*args, **kwargs)
