#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to DynamO (.xml) file formats. The main class is DynamoSnapshot, but some additional functions are defined to provide a simplified interface to this class.

The module defines:
  - DynamoSnapshot: the class the defining the file interface to this file format
  - read: shorthand for DynamoSnapshot.read_single
  - read_trajectory: shorthand for DynamoSnapshot.read_trajectory
"""

import sys, io, numpy
from .snapshot import stream_safe_open, NoSnapshotError, Snapshot

from lxml import etree as ElementTree
from lxml.etree import Element

class DynamoSnapshot(Snapshot):
    """Snapshot of a system of particles in DynamO (.xml) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def read_trajectory(*args, **kwargs):
        """Standard interface for reading trajectories should throw error because dynamo trajectories are not stored in the usual simple format."""
        raise NotImplementedError

    def read(self, path_or_file):
        """Read a snapshot from a file, overwriting any existing data.

        Args:
            path_or_file: file stream or path to read snapshot from
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format
        """
        with stream_safe_open(path_or_file) as f:
            parser = ElementTree.XMLParser(remove_blank_text=True)
            self.xml = {}
            self.xml['tree'] = ElementTree.parse(f, parser)

            self.xml['root'] = self.xml['tree'].getroot()
            self.xml['particles'] = self.xml['root'].find("ParticleData")
            self.xml['simulation'] = self.xml['root'].find('Simulation')
            self.xml['interactions'] = self.xml['simulation'].find('Interactions')

            # Particle coordinates
            self.x = numpy.array([p.find('P').attrib.values() for p in self.xml['particles']], dtype=numpy.longdouble)
            if self.n != len(self.xml['particles'].getchildren()):
                raise RuntimeError('inconsistent file!')

            # Particle species
            self.xml['genus'] = self.xml['simulation'].find('Genus')
            species = numpy.array([species.attrib['Name'] for species in self.xml['genus']])
            self.species = numpy.empty(self.n, species.dtype)
            for species in self.xml['genus']:
                id_range = species.find('IDRange').attrib
                start = int(id_range['Start'])
                end = int(id_range['End'])
                self.species[start:end+1] = species.attrib['Name']

            # System size information
            self.box = self.xml['simulation'].find('SimulationSize')
            box_lengths = numpy.array([float(self.box.attrib[dim]) for dim in ['x','y','z']])
            self.box = numpy.array([[0.,length] for length in box_lengths])
            self.volume = numpy.product(box_lengths)
            self.density = self.n / self.volume

            # Find the diameters of the particles, assuming additive interactions.
            self.diameters = numpy.empty(self.n)
            for uij in self.xml['interactions']:
                if uij.attrib['Type'] != 'HardSphere':
                    raise RuntimeException('I can only read additive hard spheres in this version!')

                uij_range = uij.find('IDPairRange')
                if uij_range.attrib['Type'] == 'Pair':
                    raise RuntimeException('I can only read additive hard spheres in this version!')
                if uij_range.attrib['Type'] == 'All':
                    self.diameters[:] = float(uij.attrib['Diameter'])
                    break

                uij_range = uij_range.getchildren()[0].attrib
                start = int(uij_range['Start'])
                end = int(uij_range['End'])
                self.diameters[start:end+1] = float(uij.attrib['Diameter'])

            # Compute exact volume fraction.
            intrinsic_volumes = numpy.pi * self.diameters**3 / 6
            self.volume_fraction = numpy.sum(intrinsic_volumes) / self.volume

            del self.__dict__['xml']

    def __str__(self):
        """To be implemented."""
        raise NotImplementedError

def read(*args, **kwargs):
    """Read a single snapshot from the disk."""
    return DynamoSnapshot.read_single(*args, **kwargs)
def read_trajectory(*args, **kwargs):
    """Read a trajectory (i.e. multiple snapshots) from the disk."""
    return DynamoSnapshot.read_trajectory(*args, **kwargs)
