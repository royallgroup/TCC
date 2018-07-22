""" Module for reading and writing snapshots from and to LAMMPS (.atom) file formats."""

import io
import numpy
import pandas
from tcc_python_scripts.file_readers.snapshot import stream_safe_open, NoSnapshotError, Snapshot


class AtomSnapshot(Snapshot):
    """Snapshot of a system of particles in LAMMPS (.atom) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def _read(self, path_or_file):
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
                    num_particles = int(f.readline())
                    self.particle_coordinates = numpy.empty((num_particles, self.dimensionality))

                # Item containing the bounding box.
                elif ' '.join(item[1:3]) == 'BOX BOUNDS':
                    num_spatial_dimensions = len(item[3:])
                    self.particle_coordinates = numpy.empty((self.num_particles, num_spatial_dimensions))
                    self.box = numpy.zeros((num_spatial_dimensions, 2), dtype=numpy.float64)

                    for dimension in range(num_spatial_dimensions):
                        boundary = f.readline().split()
                        self.box[dimension][:] = [float(b) for b in boundary]

                # Main table contains the per-atom data. Should come at the end.
                elif item[1] == 'ATOMS':
                    assert self.num_particles > 0
                    assert 1 <= self.dimensionality <= 3
                    assert self.box is not None

                    headings = item[2:]
                    assert 'id' in headings
                    assert 'x' or 'xs' in headings

                    particle_buffer = io.StringIO()
                    for i in range(self.num_particles):
                        particle_buffer.write(f.readline())
                    particle_buffer.seek(0)
                    table = pandas.read_table(particle_buffer, index_col=0, sep='\s+', names=headings,
                                              nrows=self.num_particles)

                    if 'xs' in headings:
                        cols = ['xs', 'ys', 'zs'][:self.dimensionality]
                        self.particle_coordinates = table[cols].values.copy('c').astype(numpy.float64)
                        for dimension in range(self.dimensionality):
                            side_length = self.box[dimension][1] - self.box[dimension][0]
                            self.particle_coordinates[:, dimension] *= side_length
                            self.particle_coordinates[:, dimension] += self.box[dimension][0]
                    else:
                        cols = ['x', 'y', 'z'][:self.dimensionality]
                        self.particle_coordinates = table[cols].values.copy('c').astype(numpy.float64)

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


def read(file_name):
    """ Returns a generator that reads one or more snapshots from a .atom file.

    Args:
        file_name: Name of .atom file to read.
    Returns:
         A generator object that generates Snapshots.
    """
    return AtomSnapshot.read_trajectory(file_name)
