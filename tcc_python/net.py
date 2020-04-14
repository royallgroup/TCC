import numpy as np
from sys import argv
import pathlib
from typing import Union, List


def _set_up():
    # Read XYZ file name from argv
    if len(argv) != 3:
        print("Usage tcc_python net.py directory_name priority_list")
        raise IndexError
    dir_name = argv[1]
    priority_list = argv[2]
    return dir_name, priority_list


def _load_cluster_data(num_particles: List[int], priority_list: List[str], dir_name: pathlib.Path):
    # Load data from the raw file into a dictionary
    raw_data = {}

    print("Reading data from raw files...")
    for species in priority_list:
        raw_data[species] = []
        filename = list(dir_name.glob("*" + species))
        if not filename:
            raise FileNotFoundError(f"{species} RAW file not found in directory {dir_name}")
        lines_read = 0
        for frame_particles in num_particles:
            # noinspection PyTypeChecker
            raw_data[species].append(np.genfromtxt(filename[0], skip_header=lines_read+2, invalid_raise=False,
                                                   usecols=[0], dtype='U1', max_rows=frame_particles))
            lines_read += (frame_particles + 2)

    print("Data read complete...")
    return raw_data


def _is_particle_in_cluster(particle_identifier, frame_number):
    # A cluster is found if the particle identifier is the letter C or D.
    return np.logical_or(particle_identifier[frame_number] == 'C', particle_identifier[frame_number] == 'D')


def _write_output_file(gross_percentage, net_percentage, priority_list, dir_name: pathlib.Path):
    output_path = dir_name / pathlib.Path("net_clusters.txt")
    with output_path.open('w') as output_file:
        output_file.write("Species\tGross\tNet\n")
        for species in priority_list:
            output_file.write(species + ":\t%f\t%f\n" % (gross_percentage[species], net_percentage[species]))
    print("Analysis complete. Output file written.")


def _get_particles_per_frame(dir_name: pathlib.Path, priority_list) -> List[int]:
    # Returns a list of particle numbers, one for each time frame
    num_particles = []
    filename = list(dir_name.glob("*" + priority_list[0]))
    if not filename:
        raise FileNotFoundError(f"Cannot find {priority_list[0]} file in RAW directory {dir_name}.")
    with filename[0].open('r') as xyz_file:
        line = xyz_file.readline()
        while line != "":
            num_particles.append(int(line))
            # Skip the comment and all the data
            for i in range(num_particles[-1] + 1):
                xyz_file.readline()
            line = xyz_file.readline()

    return num_particles


def net_cluster_calculation(dir_name: Union[pathlib.Path, str], priority_list: str):
    """
    Take gross TCC cluster population and calculate net cluster population.

    Args:
        dir_name: Directory containing tcc_python RAW output files.
        priority_list: List of cluster names in order of priority
    """
    if isinstance(dir_name, str):
        dir_name = pathlib.Path(dir_name)
    priority_list = priority_list.strip('()').split(", ")
    frame_particles_list = _get_particles_per_frame(dir_name, priority_list)
    total_particles = sum(frame_particles_list)
    raw_data = _load_cluster_data(frame_particles_list, priority_list, dir_name)
    gross_percentage = {}
    net_percentage = {}

    # initialise totals
    for species in priority_list:
        gross_percentage[species] = 0
        net_percentage[species] = 0

    # Loop through the clusters in priority list and process each
    for frame_number, particles_in_frame in enumerate(frame_particles_list):
        cluster_tracker = np.full(particles_in_frame, False)
        gross_list = {}
        net_list = {}
        for species in priority_list:
            gross_list[species] = _is_particle_in_cluster(raw_data[species], frame_number)
            net_list[species] = np.logical_and(gross_list[species], np.logical_not(cluster_tracker))
            cluster_tracker += gross_list[species]
            net_percentage[species] += net_list[species].sum(axis=0) / float(total_particles)
            gross_percentage[species] += gross_list[species].sum(axis=0) / float(total_particles)

    _write_output_file(gross_percentage, net_percentage, priority_list, dir_name)


if __name__ == '__main__':
    directory_name, cluster_list = _set_up()
    net_cluster_calculation(directory_name, cluster_list)
