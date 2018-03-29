import subprocess
from glob import glob
import shutil
import filecmp
import os
from subprocess import run
from platform import system


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class FileOperations:
    @staticmethod
    def build_tcc():
        try:
            if system() == "Windows":
                make = run(['cmake', '..', '-G', 'MinGW Makefiles'])
                build = run(['mingw32-make.exe'])
            elif system() == "Linux":
                make = run(['cmake', '..'])
                build = run(['make'])
            else:
                print("I dont know how to build for your system:%s", system())
                return 1
            if make.returncode == 0 and build.returncode == 0:
                return 0
            else:
                return 1
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def copy_tcc():
        # Copy the exectuable to the current directory
        try:
            if system() == "Windows":
                shutil.copy(glob("../../../bin/tcc.exe")[0], os.getcwd())
            else:
                shutil.copy(glob("../../../bin/tcc")[0], os.getcwd())
            return 0
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def run_tcc():
        try:
            return subprocess.call(glob("tcc*")[0], shell=True)
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def tidy():
        # Remove the files we have created
        os.remove(glob("tcc*")[0])
        for file in (glob("sample.xyz.rc*")):
            os.remove(file)
        return 0


class FileChecks:
    @staticmethod
    def check_static_clust():
        return filecmp.cmp("sample.static_clust", glob("sample.xyz*static_clust")[0], shallow=False)

    @staticmethod
    def check_bonds():
        return filecmp.cmp("sample.bonds", glob("sample.xyz*bonds")[0], shallow=False)

    @staticmethod
    def check_pop_per_frame():
        return filecmp.cmp("sample.pop_per_frame", glob("sample.xyz*pop_per_frame")[0], shallow=False)


def test_build():
    # Build the binary before executing tests
    build_directory = "../../build"
    if not os.path.exists(build_directory):
        os.makedirs(build_directory)
    with cd(build_directory):
        assert FileOperations.build_tcc() == 0


def test_simple_bonds():
    # Test a relatively large file with simple bonds that finds most clusters
    with cd("./simple_bonds"):
        assert FileOperations.copy_tcc() == 0
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_basic_voronoi():
    # Test a small file with multiple frames
    with cd("./basic_voronoi"):
        assert FileOperations.copy_tcc() == 0
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_cubic_voronoi_with_cell_list():
    # Test a medium file with cubic boundaries and cell list turned on
    with cd("./voronoi_cells_cubic"):
        assert FileOperations.copy_tcc() == 0
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0


def test_non_cubic_voronoi_with_cell_list():
    # Test a medium file with cubic boundaries and cell list turned on
    with cd("./voronoi_cells_non_cubic"):
        assert FileOperations.copy_tcc() == 0
        assert FileOperations.run_tcc() == 0
        assert FileChecks.check_bonds() is True
        assert FileChecks.check_pop_per_frame() is True
        assert FileChecks.check_static_clust() is True
        assert FileOperations.tidy() == 0
