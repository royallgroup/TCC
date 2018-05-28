from glob import glob
import os
import subprocess
import platform
import math
import pandas as pd


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, new_path):
        self.newPath = os.path.expanduser(new_path)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class FileOperations:
    @staticmethod
    def build_tcc():
        try:
            if platform.system() == "Windows":
                make = subprocess.run(['cmake', '..', '-G', 'MinGW Makefiles', '-DCMAKE_INSTALL_PREFIX:PATH=../bin'])
                build = subprocess.run(['mingw32-make.exe'])
                install = subprocess.run(['mingw32-make.exe', 'install'])
            elif platform.system() == "Linux" or platform.system() == 'Darwin':
                make = subprocess.run(['cmake', '..', '-DCMAKE_INSTALL_PREFIX:PATH=../bin'])
                build = subprocess.run(['make'])
                install = subprocess.run(['make', 'install'])
            else:
                print("I dont know how to build for your system:%s", platform.system())
                return 1
            if make.returncode == 0 and build.returncode == 0 and install.returncode == 0:
                return 0
            else:
                return 1
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def run_tcc():
        try:
            if platform.system() == "Windows":
                tcc_call_result = subprocess.run(glob("../../../bin/tcc.exe")[0], stdout=subprocess.DEVNULL,
                                                 stderr=subprocess.DEVNULL)
            else:
                tcc_call_result = subprocess.run(glob("../../../bin/tcc")[0], stdout=subprocess.DEVNULL,
                                                 stderr=subprocess.DEVNULL)
            return tcc_call_result.returncode
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def tidy():
        # Remove the files we have created
        for file in (glob("sample.xyz.rc*")):
            os.remove(file)
        return 0


class FileChecks:
    @staticmethod
    def check_bonds():
        measured_results = pd.read_table(glob('sample.xyz*.bonds')[0], skiprows=1, header=None)
        measured_results.fillna(0., inplace=True)
        known_results = pd.read_table('sample.bonds', skiprows=1, header=None)
        known_results.fillna(0., inplace=True)

        for measured_row in measured_results.iterrows():
            if not pd.Series.equals(known_results.iloc[measured_row[0]], measured_row[1]):
                return False

        return True

    @staticmethod
    def check_static_clust():
        measured_results = pd.read_table(glob('sample.xyz*.static_clust')[0], index_col='Cluster type', skiprows=1)
        measured_results.fillna(0., inplace=True)
        known_results = pd.read_table('sample.static_clust', index_col='Cluster type', skiprows=1)
        known_results.fillna(0., inplace=True)

        return FileChecks.compare_cluster_files_by_row(measured_results, known_results)

    @staticmethod
    def check_pop_per_frame():
        measured_results = pd.read_table(glob("sample.xyz*pop_per_frame")[0], index_col='frame')
        measured_results.fillna(0., inplace=True)
        known_results = pd.read_table('sample.pop_per_frame', index_col='frame')
        known_results.fillna(0., inplace=True)

        return FileChecks.compare_cluster_files_by_row(measured_results, known_results)

    @staticmethod
    def compare_cluster_files_by_row(measured_results, known_results):
        for column_name in measured_results:
            for measured_particle_type in measured_results[column_name].items():
                particle_type_found = 0
                for known_particle_type in known_results[column_name].items():
                    if measured_particle_type[0] == known_particle_type[0]:
                        particle_type_found = 1
                        if math.isclose(measured_particle_type[1], known_particle_type[1], abs_tol=0.00001) == False:
                            print("\nColumn: \"{0}\" Particle type: \"{1}\", does not match known quantity."
                                  .format(column_name, measured_particle_type[0]))
                            return False
                        break
                if particle_type_found == 0:
                    print("\nParticle type {0} not found in TCC output".format(measured_particle_type[0]))
                    return False
        return True


def test_build():
    # Build the binary before executing tests
    build_directory = "build"
    if not os.path.exists(build_directory):
        os.makedirs(build_directory)
    with cd(build_directory):
        assert FileOperations.build_tcc() == 0


def test_simple_bonds():
    # Test a relatively large file with simple bonds that finds most clusters
    with cd("./test/integration_tests/simple_bonds"):
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_basic_voronoi():
    # Test a small file with multiple frames
    with cd("./test/integration_tests/basic_voronoi"):
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_cubic_voronoi_with_cell_list():
    # Test a medium file with cubic boundaries and cell list turned on
    with cd("./test/integration_tests/voronoi_cells_cubic"):
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_non_cubic_voronoi_with_cell_list():
    # Test a medium file with cubic boundaries and cell list turned on
    with cd("./test/integration_tests/voronoi_cells_non_cubic"):
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0
