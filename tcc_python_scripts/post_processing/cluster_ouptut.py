# A debugging script used to map cluster output files to XYZ files for visualisation.
def main():
	import numpy as np
	
	xyz_file = np.loadtxt("../sample_ka.xyz", skiprows=2, usecols=[1,2,3])
	data = np.loadtxt("sample_ka.xyz.rcAA2.rcAB2.rcBB2.Vor1.fc1.PBCs0.clusts_sp3c", skiprows=1, dtype=np.int)
	num_particles = len(xyz_file)
	
	with open("outputfile.xyz", 'w') as outputfile:
		output_char = ""
		for line in data:
			outputfile.write("{}\ncomment\n".format(num_particles))
			for particle in range(num_particles):
				if line[0] == particle or line[1] == particle or line[2] == particle:
					output_char = "A"
				elif line[3] == particle or line[4] == particle:
					output_char = "B"
				else:
					output_char = "C"
				outputfile.write("{}\t{}\t{}\t{}\n".format(output_char, xyz_file[particle][0], xyz_file[particle][1], xyz_file[particle][2]))
				
if __name__ == "__main__":
	main()