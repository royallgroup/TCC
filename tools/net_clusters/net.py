import numpy as np
from sys import argv, exit
from glob import glob


def set_up():
    # Read XYZ file name from argv
    if len(argv) != 2:
        print("Usage python net.py file.xyz")
        exit(0)
    xyz_name = argv[1]
    return xyz_name


def load_cluster_data(num_particles, priority_list, xyz_name):
    # Load data from the raw file into a dictionary
    raw_data = {}

    for species in priority_list:
        print("Reading", species, "...")
        filename = glob("raw_output/" + xyz_name + "*" + species)[0]
        raw_data[species] = np.genfromtxt(filename, skip_header=2, invalid_raise=False, usecols=[0], dtype='U1',
                                          max_rows=num_particles)

    return raw_data


def is_particle_in_cluster(particle_identifier):
    # A cluster is found if the particle identifier is the letter C or D.
    return np.logical_or(particle_identifier == 'C', particle_identifier == 'D')


def write_output_file(gross_percentage, net_percentage, priority_list, xyz_name):
    with open(xyz_name + "_net.txt", 'w') as output_file:
        output_file.write("Species\tGross\tNet\n")
        for species in priority_list:
            output_file.write(species + ":\t%f\t%f\n" % (gross_percentage[species], net_percentage[species]))


def main():
    xyz_name = set_up()

    # Modify this list in order to have your favourite hierarchy
    priority_list = ['FCC', '13A', '12E', '11F', '10B', '9B', '8B', 'sp5c', 'sp4c', 'sp3c']

    # read number of particles
    with open(glob("raw_output/" + xyz_name + "*" + priority_list[0])[0], 'r') as xyz:
        num_particles = int(xyz.readline())

    gross_list = {}
    net_list = {}
    gross_percentage = {}
    net_percentage = {}

    raw_data = load_cluster_data(num_particles, priority_list, xyz_name)

    cluster_tracker = np.full(num_particles, False)

    # Loop through the clusters in priority list and process each
    for species in priority_list:
        gross_list[species] = is_particle_in_cluster(raw_data[species])
        net_list[species] = np.logical_and(gross_list[species], np.logical_not(cluster_tracker))
        cluster_tracker += gross_list[species]
        net_percentage[species] = net_list[species].sum(axis=0) / float(num_particles)
        gross_percentage[species] = gross_list[species].sum(axis=0) / float(num_particles)

    write_output_file(gross_percentage, net_percentage, priority_list, xyz_name)


if __name__ == '__main__':
    main()

main()
