"""Python interface to the TCC executable."""

import os
import tempfile
import shutil
import numpy
import pandas
import subprocess
import platform
from glob import glob

from tcc_python_scripts.tcc import structures, wrapper

class DynamicCluster:
    def __init__(self, particles, time):
        self.particles = numpy.sort(particles)
        self.particles = tuple(particles.tolist())
        self.creation_time = time
        self.last_seen_time = time

    @property
    def table_entry(self):
        return numpy.concatenate(([self.lifetime, self.creation_time, self.last_seen_time], self.particles))

    @property
    def lifetime(self):
        return 1 + self.last_seen_time - self.creation_time

    def __eq__(self, other):
        same = self.particles == other.particles
        if same:
            self.creation_time = min(self.creation_time, other.creation_time)
            self.last_seen_time = max(self.last_seen_time, other.last_seen_time)
            other.creation_time = self.creation_time
            other.last_seen_time = self.last_seen_time
        return same

    def __hash__(self):
        return hash(self.particles)

    def __repr__(self):
        return '<DynamicCluster t0={} t1={}>'.format(self.creation_time, self.last_seen_time)

class DynamicTCC(wrapper.TCCWrapper):
    def save(self, destination):
        """Save results of TCC analysis to a specified path.

        The working directory is deleted when the wrapper falls out of scope, so this must
        be called to retain all the data.

        Args:
            destination: new directory to save data
        """
        os.makedirs(destination)
        shutil.copy2(os.path.join(self.working_directory, 'inputparameters.ini'), destination)

        self.static_populations.to_csv(os.path.join(destination, 'statics.pop.csv'), '\t')
        self.static_populations_err.to_csv(os.path.join(destination, 'statics.err.csv'), '\t')

        for structure, clusters in self.clusters.items():
            out_path = os.path.join(destination, 'clusters.%s' % structure)
            numpy.save(out_path, clusters)

    def run(self, trajectory, decay_threshold=1):
        """Run the TCC on a trajectory, averaging the static data and performing the dynamic TCC
        algorithm to determine cluster lifetimes.

        Args:
            trajectory: container of snapshots to analyse.
                The snapshots must contain the box information and coordinates.
            decay_threshold: threshold number of frames that a structure must disappear for before
                it is recognised as having dissociated.
        """
        self.static_populations = None
        static_populations_squ = None
        self.clusters = {}
        surviving_clusters = {}

        for structure in self.active_clusters:
            self.clusters[structure] = []
            surviving_clusters[structure] = set()

        for frame,snap in enumerate(trajectory):
            super().run(snap.box_dimensions, snap.x, output_clusters=True)

            # Average the static data.
            statics = self._parse_static_clusters()
            if self.static_populations is None:
                self.static_populations = statics
                static_populations_squ = statics**2
            else:
                self.static_populations += statics
                static_populations_squ += statics**2

            # The dynamic TCC algorithm to determine the lifetimes of the clusters.
            for structure in self.active_clusters:
                for cluster in self._parse_cluster_file(structure):
                    cluster = DynamicCluster(cluster, frame)
                    surviving_clusters[structure].add(cluster)

                decayed = []
                for cluster in surviving_clusters[structure]:
                    if frame >= (cluster.last_seen_time + decay_threshold):
                        decayed += [cluster]
                for cluster in decayed:
                    self.clusters[structure] += [cluster.table_entry]
                    surviving_clusters[structure].discard(cluster)

        nframes = frame+1
        self.static_populations /= nframes
        self.static_populations_err = (static_populations_squ - self.static_populations**2)**0.5 / (nframes*(nframes-1))**0.5

        for structure in self.active_clusters:
            self.clusters[structure] = numpy.array(self.clusters[structure], dtype=int)

    def lifetimes(self, structure):
        """Lifetimes of clusters from the dynamic TCC analysis.

        Args:
            structure: which cluster type
        Returns:
            numpy array (int) of the lifetimes of each cluster from the trajectory, in units of frames.
        """
        return self.clusters[structure][:,0]
