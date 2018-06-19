import numpy as np
from sys import argv, exit
from glob import glob
import os


def set_up():
    # Read XYZ file name from argv
    if len(argv) != 3:
        print("Usage python net.py directory_name priority_list")
        exit(0)
    dir_name = argv[1]
    pritority_list = argv[2]
    return dir_name, pritority_list


def load_cluster_data(num_particles, priority_list, dir_name):
    # Load data from the raw file into a dictionary
    raw_data = {}

    print("Reading data from raw files...")
    for species in priority_list:
        raw_data[species] = []
        filename = glob(dir_name + "/*" + species)[0]
        lines_read = 0
        for frame_particles in num_particles:
            raw_data[species].append(np.genfromtxt(filename, skip_header=lines_read+2, invalid_raise=False,
                                                   usecols=[0], dtype='U1', max_rows=frame_particles))
            lines_read += (frame_particles + 2)

    print("Data read complete...")
    return raw_data


def is_particle_in_cluster(particle_identifier, frame_number):
    # A cluster is found if the particle identifier is the letter C or D.
    return np.logical_or(particle_identifier[frame_number] == 'C', particle_identifier[frame_number] == 'D')


def write_output_file(gross_percentage, net_percentage, priority_list, dir_name):
    with open(dir_name + "/net_clusters.txt", 'w') as output_file:
        output_file.write("Species\tGross\tNet\n")
        for species in priority_list:
            output_file.write(species + ":\t%f\t%f\n" % (gross_percentage[species], net_percentage[species]))
    print("Analysis complete. Output file written.")


def get_particles_per_frame(dir_name, priority_list):
    # Returns a list of particle numbers, one for each time frame
    num_particles = []
    filename = glob(dir_name + "/*" + priority_list[0])[0]
    with open(filename, 'r') as xyz_file:
        line = xyz_file.readline()
        while line != "":
            num_particles.append(int(line))
            # Skip the comment and all the data
            for i in range(num_particles[-1] + 1):
                xyz_file.readline()
            line = xyz_file.readline()

    return num_particles


def net_cluster_calculation(dir_name, priority_list):
    # read number of particles and the data
    priority_list = priority_list.strip('()').split(", ")
    frame_particles_list = get_particles_per_frame(dir_name, priority_list)
    total_particles = sum(frame_particles_list)
    raw_data = load_cluster_data(frame_particles_list, priority_list, dir_name)
    gross_percentage = {}
    net_percentage = {}

    # intitialse totals
    for species in priority_list:
        gross_percentage[species] = 0
        net_percentage[species] = 0

    # Loop through the clusters in priority list and process each
    for frame_number, particles_in_frame in enumerate(frame_particles_list):
        cluster_tracker = np.full(particles_in_frame, False)
        gross_list = {}
        net_list = {}
        for species in priority_list:
            gross_list[species] = is_particle_in_cluster(raw_data[species], frame_number)
            net_list[species] = np.logical_and(gross_list[species], np.logical_not(cluster_tracker))
            cluster_tracker += gross_list[species]
            net_percentage[species] += net_list[species].sum(axis=0) / float(total_particles)
            gross_percentage[species] += gross_list[species].sum(axis=0) / float(total_particles)

    write_output_file(gross_percentage, net_percentage, priority_list, dir_name)


if __name__ == '__main__':
    directory_name, cluster_list = set_up()
    net_cluster_calculation(directory_name, cluster_list)
