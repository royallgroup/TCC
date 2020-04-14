""" Module for reading and writing snapshots from and to DynamO (.xml) file formats."""

import numpy
import tcc_python.file_readers.snapshot as snapshot
import xml.etree.ElementTree as ElementTree


class NonAdditiveError(RuntimeError):
    """
    Error raised if the given interaction range is a cross-species
    interaction, so potentially the dynamo file describes a non-additive potential.
    """
    pass


class DynamoSnapshot(snapshot.Snapshot):
    """Snapshot of a system of particles in DynamO (.xml) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    @staticmethod
    def _is_hard_sphere(interaction):
        """Determine whether an interaction encoded in the dynamo XML tree is a hard sphere
        interaction or not.

        There are two ways to encode hard spheres: either explicitly as a hard sphere, or as
        a square well with zero well depth (although the latter should be avoided as it results
        in a less efficient simulation). This helper function tests for both eventualities.

        Args:
            interaction: xml entry for the interaction
        Returns:
            boolean stating whether the interaction is a hard sphere
        """
        if interaction.attrib['Type'] == 'HardSphere':
            return True
        return interaction.attrib['Type'] == 'SquareWell' and interaction.attrib['WellDepth'] == '0'

    def _assign_diameters(self, uij_range, diameter):
        """Assign diameters to the particles across the given ID range.

        Diameters will only be assigned if the range describes a single-species interaction.
        Any other range is potentially a nonadditive interaction so an exception is raised as
        a warning.

        Args:
            uij_range: xml entry describing the range of the interaction
            diameter: length scale of these hard spheres
        Raises:
            NonadditiveError: if the given range is a cross-species interaction, so potentially
                              the dynamo file describes a nonadditive potential.
        """
        range_type = uij_range.attrib['Type']

        if range_type == 'All':
            self.diameters[:] = diameter

        elif range_type == 'Single':
            uij_range = uij_range.getchildren()[0].attrib
            start = int(uij_range['Start'])
            end = int(uij_range['End'])
            self.diameters[start:end+1] = diameter

        elif range_type == 'Pair':
            range0, range1 = uij_range.getchildren()
            start0, start1 = [int(r.attrib['Start']) for r in [range0, range1]]
            end0, end1 = [int(r.attrib['End']) for r in [range0, range1]]

            if start0 == start1 and end0 == end1:
                assigning_diameters = self.diameters[start0:end0+1]
                if numpy.all(numpy.isnan(assigning_diameters)):
                    assigning_diameters[:] = diameter
                    return
                else:
                    raise RuntimeError('can only process additive hard spheres!')
            else:
                raise NonAdditiveError

        else:
            raise RuntimeError('unknown range type encountered: %s' % range_type)

    def _verify_cross_species_interactions(self, uij_range, diameter, eps=1e-12):
        """Check whether a cross-species hard sphere interaction is consistent with the known
        diameters of each species, assuming additive interactions.

        Additive cross-species interactions have size \sigma = 0.5*(\sigma_i + \sigma_j) where
        \sigma_{i,j} are the diameters of the individual species. If the diameter does not
        match this expectation then the interaction is not additive.

        Args:
            uij_range: xml entry describing the range of the interaction
            diameter: length scale of this cross-species hard sphere interaction
            eps: machine tolerance for the test of additivity
        Raises:
            NonadditiveError: if the interaction is not additive
        """

        if uij_range.attrib['Type'] != 'Pair':
            return

        range0, range1 = uij_range.getchildren()
        start0, start1 = [int(r.attrib['Start']) for r in [range0, range1]]
        end0, end1 = [int(r.attrib['End']) for r in [range0, range1]]

        if start0 == start1 and end0 == end1:
            return

        diameters0 = self.diameters[start0:end0+1]
        diameters1 = self.diameters[start1:end1+1]
        if numpy.all(diameters0 == diameters0[0]) and numpy.all(diameters1 == diameters1[0]):
            diameters0 = diameters0[0]
            diameters1 = diameters1[0]
            additive_diameter = 0.5*(diameters0 + diameters1)

            if abs(additive_diameter - diameter) < eps:
                return

        raise NonAdditiveError('cross interaction between ranges (%d,%d) and (%d,%d) is non additive!' % (start0, end0, start1, end1))

    def _read(self, path_or_file):
        """Read a snapshot from a file, overwriting any existing data.

        Args:
            path_or_file: file stream or path to read snapshot from
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format or the data does not describe an
                              additive hard sphere system
        """
        with snapshot.stream_safe_open(path_or_file) as f:
            parser = ElementTree.XMLParser()
            self.xml = {}
            try:
                self.xml['tree'] = ElementTree.parse(f, parser)
            except ElementTree.ParseError:
                print("Unable to read dynamo file.")
                raise snapshot.NoSnapshotError

            self.xml['root'] = self.xml['tree'].getroot()
            self.xml['particles'] = self.xml['root'].find("ParticleData")
            self.xml['simulation'] = self.xml['root'].find('Simulation')
            self.xml['interactions'] = self.xml['simulation'].find('Interactions')

            # Particle coordinates
            self.particle_coordinates = numpy.array([list(p.find('P').attrib.values()) for p in self.xml['particles']], dtype=numpy.longdouble)
            if self.num_particles != len(self.xml['particles'].getchildren()):
                raise RuntimeError('inconsistent file!')

            # Particle species
            self.xml['genus'] = self.xml['simulation'].find('Genus')
            species = numpy.array([species.attrib['Name'] for species in self.xml['genus']])
            self.species = numpy.empty(self.num_particles, species.dtype)
            for species in self.xml['genus']:
                id_range = species.find('IDRange').attrib
                start = int(id_range['Start'])
                end = int(id_range['End'])
                self.species[start:end+1] = species.attrib['Name']

            # System size information
            self.box = self.xml['simulation'].find('SimulationSize')
            box_lengths = numpy.array([float(self.box.attrib[dim]) for dim in ['x', 'y', 'z']])
            self.box = numpy.array([[0., length] for length in box_lengths])
            self.volume = numpy.product(box_lengths)
            self.density = self.num_particles / self.volume

            # Find the diameters of the particles, assuming additive interactions.

            # Initialise diameters to NaN: if any are still NaN at the end we know some were
            # uninitialised and the data file is incomplete
            self.diameters = numpy.full(self.num_particles, numpy.nan)

            # Assign diameters
            definitely_additive = True
            for uij in self.xml['interactions']:
                if not self._is_hard_sphere(uij):
                    raise RuntimeError('can only process additive hard spheres!')

                uij_range = uij.find('IDPairRange')
                diameter = float(uij.attrib['Diameter'])
                try:
                    self._assign_diameters(uij_range, diameter)
                except NonAdditiveError:
                    definitely_additive = False

            # Check that all diameters were assigned above
            uninitialised = numpy.isnan(self.diameters)
            if numpy.any(uninitialised):
                raise RuntimeError('some diameters were uninitialised! indices=%r' %
                                   numpy.where(uninitialised))

            if not definitely_additive:
                for uij in self.xml['interactions']:
                    uij_range = uij.find('IDPairRange')
                    diameter = float(uij.attrib['Diameter'])
                    self._verify_cross_species_interactions(uij_range, diameter)

            # Compute exact volume fraction.
            intrinsic_volumes = numpy.pi * self.diameters**3 / 6
            self.volume_fraction = numpy.sum(intrinsic_volumes) / self.volume

            del self.__dict__['xml']

    def __str__(self):
        """To be implemented."""
        raise NotImplementedError


def read(file_name):
    """ Read a snapshot from the dynamo file.
    At the moment only a single frame can be read.

    Args:
        file_name: Name of Dynamo file to read.
    Returns:
         A generator object that generates Snapshots.
    """
    return DynamoSnapshot.read_trajectory(file_name)
