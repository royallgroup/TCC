#!/usr/bin/env python3
"""
Module defining interface for reading and writing snapshots in various file formats.

The module defines:
  - stream_safe_open: a context manager to facilitate parsing both open file streams and paths with the same interface
  - NoSnapshotError: exception raised when a snapshot could not be read from the file
  - snapshot: the class containing the general snapshot interface
"""

import sys, numpy, contextlib

@contextlib.contextmanager
def stream_safe_open(path_or_file, mode='r'):
    """Context manager for parsers which accept an open file stream or a file path to open.

    Args:
        path_or_file: either an open file stream (in which case the context manager does nothing and returns it) or a path (in which case the context manager will open this, return a stream and then clean up)
        mode: mode to open file in
    """
    if isinstance(path_or_file, str):
        f = file_to_close = open(path_or_file, mode)
    else:
        f = path_or_file
        file_to_close = None

    try:
        yield f
    finally:
        if file_to_close:
            file_to_close.close()

class NoSnapshotError(RuntimeError):
    """Exception raised when not able to read a snapshot from a file."""

class Snapshot:
    """Snapshot of a system of particles.

    Variables:
        n: number of particles
        d: dimensionality of configuration space
        x: particle coordinates (n by d container)
        box: box containing the particles (d by 2 container)
        species: labels of the particle species (string or container of strings)
        time: time or frame of the snapshot within a trajectory
    """

    def __init__(self, x=numpy.empty((0,0)), box=None, species=None, time=0):
        """Create a new snapshot.

        Args:
            x: particle coordinates
            box: box dimensions containing the particles
            species: labels of the particle species
            time: time or frame of the snapshot within a trajectory
        """
        self.x = x
        self.box = box
        if species is None: self.species = ['A']*self.n
        else: self.species = species
        self.time = time

    def copy(self):
        """Return a deep copy of the snapshot."""
        return snapshot(self.x.copy(), self.box.copy(), self.species.copy(), self.time)

    @property
    def n(self):
        """Number of particles in snapshot."""
        return len(self.x)

    @property
    def d(self):
        """Dimensionality of configuration space."""
        return self.x.shape[1]

    @classmethod
    def read_single(cls, path_or_file):
        """Read a single snapshot from the disk.

        Example:
        >>> snapshot.read('snapshot.atom')
        <snapshot n=10976 t=0>

        Args:
            cls: derived class defining specific file format
            path_or_file: file stream or path to read snapshot from
        Returns:
            snapshot: the snapshot read from disk
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format
        """
        with stream_safe_open(path_or_file) as f:
            snap = cls()
            snap.read(f)
            return snap

    @classmethod
    def read_trajectory(cls, path_or_file, max_frames=None):
        """Read a trajectory (i.e. multiple snapshots) from the disk.

        Example where trajectory.atom is an atom file containing two snapshots:
        >>> list(read('trajectory.atom', 2))
        [<snapshot n=10976 t=0>, <snapshot n=10976 t=1>]

        Args:
            cls: derived class defining specific file format
            path_or_file: file stream or path to read trajectory from
        Returns:
            trajectory (generator): generator iterating through the snapshots in the trajectory
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format
        """
        with stream_safe_open(path_or_file) as f:
            frames = 0
            while True:
                try: snap = cls.read_single(f)
                except NoSnapshotError: break

                yield snap
                frames += 1
                if max_frames is not None and frames is max_frames: break

    def write(self, out=sys.stdout):
        """Dump the snapshot to a file in LAMMPS (.atom) format.

        Args:
            out: file or path to write the snapshot to
        """
        with stream_safe_open(out, 'w') as f:
            f.write(str(self))
            f.write('\n')

    def __repr__(self):
        """Internal representation of the object for printing to debugger."""
        return '<snapshot n=%r t=%r>' % (self.n, self.time)

    def __str__(self):
        """String representation of the snapshot in the chosen format not implemented in the base class, and must be written for specific formats."""
        raise NotImplementedError
