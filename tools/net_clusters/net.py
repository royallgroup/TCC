import numpy as np
import sys
import glob


def read_and_split(filename, num_particles):

    data = np.genfromtxt(filename, skip_header=1, invalid_raise=False, usecols=[0], dtype='U1')
    frames = data.reshape((num_particles + 2, len(data) / (num_particles + 2)), order='F')

    return frames[2:, :]


def cluster(info):
    return (info == 'C') + (info == 'D')


def main():
    if len(sys.argv) < 2:
        print("Usage python net.py file.xyz")
        sys.exit(0)
    # input a file name
    xyz_name = sys.argv[1]

    # Modify this list in order to have your favourite hierarchy
    priority_list = ['FCC', '13A', '12E', '11F', '10B', '9B', '8B', '7A', '6A', '5A']

    # read number of particles
    with open(xyz_name, 'r') as xyz:
        num_particles = int(xyz.readline())

    TCCinfo = {}
    netTCCinfo = {}
    gross, net = {}, {}

    for species in priority_list:
        print("Reading", species, "...")
        TCCinfo[species] = read_and_split(glob.glob("raw_output/" + xyz_name + "*" + species)[0], num_particles)

    # Load first cluster in priority list
    test = cluster(TCCinfo[priority_list[0]])
    netTCCinfo[priority_list[0]] = test
    net[priority_list[0]] = netTCCinfo[priority_list[0]].sum(axis=0) / float(num_particles)
    gross[priority_list[0]] = net[priority_list[0]]

    # Load the rest of the clusters in priority list
    for species in priority_list[1:]:
        c = cluster(TCCinfo[species])
        netTCCinfo[species] = np.logical_not(test) * c
        test += c
        net[species] = netTCCinfo[species].sum(axis=0) / float(num_particles)
        gross[species] = c.sum(axis=0) / float(num_particles)

    for species in priority_list:
        if net[species].shape[0] > 1:
            net[species] = np.average(net[species])
            gross[species] = np.average(gross[species])

    with open(xyz_name+"_net", 'w') as fout:
        fout.write("Species\tGross\tNet\n")
        for species in priority_list:
            fout.write(species + ":\t%f\t%f\n" % (gross[species], net[species]))


main()
