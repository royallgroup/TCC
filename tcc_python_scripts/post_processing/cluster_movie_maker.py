""" A script to take the RAW file output from the TCC and produce a combined
    cluster XYZ for rendering. Takes a cluster list argument and labels
    particles according to the list priority."""

import sys
import io
import string


def add_cluster_to_xyz(xyz_frame, particle_types, cluster_number):
    """
    Given a list of particle types, overwrite the cluster type in xyz_frame if
    a particle is in a cluster.

    Args:
        xyz_frame (cluster_movie_maker.Snapshot): A snapshot of the system.
        particle_types (list of string): Whether particles are in a cluster
            or not. "A" and "B" correspond to particles not in a cluster while
            "C" and "D" correspond to particles in a cluster
        cluster_number (integer): The index of the cluster in the priority

    Returns:
        cluster_movie_maker.Snapshot: an updated XYZ snapshot.
    """
    for index, particle in enumerate(particle_types):
        if particle == "C" or particle == "D":
            xyz_frame.particle_species[index] = \
                string.ascii_uppercase[cluster_number + 1]
    return xyz_frame


def prepare_output_file(output_filename):
    """
    Open the output file once at the beginning of the analysis to delete
    any previous copy.
    """
    open(output_filename, 'w').close()


def main(xyz_name, raw_stub, cluster_list):
    """
    The main loop.

    Returns:
        int: 0 if script ran successfully.
    """
    output_filename = "output.xyz"

    cluster_list = list(reversed(cluster_list.split()))
    raw_file_handles = open_raw_files(cluster_list, raw_stub)
    xyz_reader = XyzFileReader(xyz_name)

    print("Particles not in any cluster are labelled with the letter A")

    prepare_output_file(output_filename)

    for frame_number, xyz_frame in enumerate(xyz_reader):
        for index, raw_file in enumerate(raw_file_handles):
            cluster_types = raw_file.get_frame()
            xyz_frame = add_cluster_to_xyz(xyz_frame, cluster_types, index)
            if frame_number == 0:
                print("Cluster type {} is labelled with the letter {}".format(
                    raw_file.cluster_name, string.ascii_uppercase[index + 1]))
        xyz_frame.write_xyz(output_filename)


def open_raw_files(cluster_list, raw_stub):
    """
    Given a list of clusters, open the corresponding TCC RAW files.

    Args:
        cluster_list (list of string): TCC cluster names in reverse order
            of priority.
        raw_stub (string):  The base of the RAW file name relative to the
            working directory.

    Returns:
        list of RawFileReader: Objects representing the open RAW files.
    """
    raw_file_handles = []
    for cluster_type in cluster_list:
        raw_file_handles.append(RawFileReader(raw_stub, cluster_type))
    return raw_file_handles


def process_arguments():
    """
    Process the command line arguments and return variables with the values.

    Returns:
        strings: Processed command line parameters.
    """
    if len(sys.argv) != 4:
        print("Syntax: ./cluster_movie_maker.py simulation_data.xyz "
              "simulation_data.xyz.raw_ cluster_priority_list")
        sys.exit()
    xyz_name = sys.argv[1]
    raw_stub = sys.argv[2]
    cluster_list = sys.argv[3]

    return xyz_name, raw_stub, cluster_list


class XyzFileReader:
    """
    Generator for xyz files. Opens a file handle when the object is
    instantiated. Iteration then returns the file frame by frame.
    """
    def __init__(self, file_name):
        """
        On instantiation, open a file handle to file_name.

        Args:
            file_name (string): The xyz file to open.
        """
        self.file_handle = open(file_name, 'r')

    def __iter__(self):
        """
        Reads a single frame from the xyz and returns it.

        Returns:
            Snapshot: A snapshot object representing a single frame from the
                xyz.
        """
        line = self.file_handle.readline()
        while line != "" and line != "\n":
            xyz_snapshot = Snapshot()
            xyz_snapshot.num_particles = self.process_num_particles(line)
            xyz_snapshot.comment = self.file_handle.readline()

            for particle in range(xyz_snapshot.num_particles):
                line = self.file_handle.readline().split()
                # All particles are labelled "A" as this represents not being
                # present in a cluster. This value is later overwritten if it
                # is in a cluster.
                xyz_snapshot.particle_species.append("A")
                xyz_snapshot.x_coordinates.append(float(line[1]))
                xyz_snapshot.y_coordinates.append(float(line[2]))
                xyz_snapshot.z_coordinates.append(float(line[3]))
            yield xyz_snapshot
            # Read the first line of the next frame
            line = self.file_handle.readline()

    @staticmethod
    def process_num_particles(line):
        """
        Process the first line of an XYZ frame to make sure that it is a
        valid number.

        Args:
            line (string): The first line of an xyz frame.

        Returns:
            integer: The number of particles from line.
        """
        try:
            num_particles = int(line)
        except TypeError:
            print("Error reading number of particles from XYZ file.")
            raise TypeError
        return num_particles


class RawFileReader:
    """
    Reader for TCC .RAW files. Instantiation opens a file handle and the get
    frame method then allows the file to be read a frame at a time.
    """
    def __init__(self, file_stub, cluster_type):
        """
        Open the RAW file specified by file_name.

        Args:
            file_stub (string): The base of the file name relative to the
         working directory.
            cluster_type (string): The TCC cluster type e.g. "10A"
        """
        file_name = "{}{}".format(file_stub, cluster_type)
        self.file_handle = open(file_name, "r")
        self.num_particles = 0
        self.cluster_name = cluster_type

    def get_frame(self):
        """
        Get a single frame of data from the raw file.

        Returns:
            list of string: Whether the particles are in a cluster or not.
        """
        self.num_particles = int(self.file_handle.readline())
        # Skip the comment line
        _ = self.file_handle.readline()

        type_list = []

        for particle in range(self.num_particles):
            line = self.file_handle.readline().split()
            type_list.append(line[0])
        return type_list


class Snapshot:
    """
    Object representing a single configuration of a system with particle types
    and coordinates.
    """
    def __init__(self):
        """
        Initialise variables.
        """
        self.num_particles = 0
        self.comment = ""
        self.particle_species = []
        self.x_coordinates = []
        self.y_coordinates = []
        self.z_coordinates = []

    def __str__(self):
        """"
        Returns:
            string: String representation of the Snapshot.
        """
        buffer = io.StringIO()
        for particle in range(self.num_particles):
            buffer.write("{} {} {} {}\n".format(self.particle_species[particle],
                                                self.x_coordinates[particle],
                                                self.y_coordinates[particle],
                                                self.z_coordinates[particle]))
        return buffer.getvalue()

    def write_xyz(self, output_name):
        """
        Write the Snapshot to an xyz file.
        """
        with open(output_name, 'a') as output_file:
            output_file.write("{}\n".format(self.num_particles))
            output_file.write("{}".format(self.comment))
            for particle in range(self.num_particles):
                output_file.write("{}\t{}\t{}\t{}\n".
                                  format(self.particle_species[particle],
                                         self.x_coordinates[particle],
                                         self.y_coordinates[particle],
                                         self.z_coordinates[particle]))


if __name__ == "__main__":

    XYZ_NAME, RAW_STUB, CLUSTER_LIST = process_arguments()

    main(XYZ_NAME, RAW_STUB, CLUSTER_LIST)
