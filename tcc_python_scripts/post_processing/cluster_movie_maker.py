"""A script to take the RAW file output from the TCC and produce a combined
cluster XYZ for rendering. Takes a cluster list argument and labels particles
according to the list priority."""

import sys
import io
import string


def add_cluster_to_xyz(xyz_frame, cluster_types):
    return 0


def write_xyz(xyz_frame):
    pass


def main():
    process_arguments()
    raw_file_handles = []
    for cluster_type in (sys.argv[3]):
        raw_file_handles.append(RawFileReader(sys.argv[2], cluster_type))
    xyz_file = XyzFileReader(sys.argv[1])
    for xyz_frame in xyz_file:
        for index, raw_file in enumerate(raw_file_handles):
            cluster_types = raw_file.get_frame()
            xyz_frame = add_cluster_to_xyz(xyz_frame, cluster_types)
            print("Cluster type {} is output with the letter {}"
                  .format(raw_file.cluster_name,
                          string.ascii_uppercase()[index]))
        write_xyz(xyz_frame)


def process_arguments():
    if len(sys.argv) != 4:
        print("Syntax: ./michael11Afinder.py simulationxmolfile.xyz "
              "simulationxmol.raw_ cluster_priority_list")
        sys.exit()


class XyzFileReader:
    def __init__(self, file_name):
        self.file_handle = open(file_name, 'r')

    def __iter__(self):
        line = self.file_handle.readline()
        while line != "" and line != "\n":
            xyz_snapshot = Snapshot()
            xyz_snapshot.process_num_particles(line)
            xyz_snapshot.comment = self.file_handle.readline()

            for particle in range(xyz_snapshot.num_particles):
                line = self.file_handle.readline().split()
                xyz_snapshot.particle_species.append(line[0])
                xyz_snapshot.x_coordinates.append(float(line[1]))
                xyz_snapshot.y_coordinates.append(float(line[2]))
                xyz_snapshot.z_coordinates.append(float(line[3]))
            yield xyz_snapshot
            # Read the first line of the next frame
            line = self.file_handle.readline()


class RawFileReader:
    def __init__(self, file_stub, cluster_type):
        file_name = "{}{}".format(file_stub, cluster_type)
        self.file_handle = open(file_name, "r")
        self.num_particles = 0
        self.cluster_name = cluster_type

    def get_frame(self):
        self.num_particles = int(self.file_handle.readline())
        # Skip the comment line
        _ = self.file_handle.readline()

        type_list = []

        for particle in range(self.num_particles):
            line = self.file_handle.readline().split()
            type_list.append(line[0])
        yield type_list


class Snapshot:

    def __init__(self):
        self.num_particles = 0
        self.comment = ""
        self.particle_species = []
        self.x_coordinates = []
        self.y_coordinates = []
        self.z_coordinates = []

    def __str__(self):
        buffer = io.StringIO()
        for particle in range(self.num_particles):
            buffer.write("{} {} {} {}\n".format(self.particle_species[particle],
                                                self.x_coordinates[particle],
                                                self.y_coordinates[particle],
                                                self.z_coordinates[particle]))
        return buffer.getvalue()

    def process_num_particles(self, line):
        try:
            self.num_particles = int(line)
        except TypeError:
            print("Error reading number of particles from XYZ file.")
            raise TypeError


if __name__ == "__main__":
    main()
