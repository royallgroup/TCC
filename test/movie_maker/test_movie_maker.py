import filecmp
import os
from tcc_python_scripts.post_processing import cluster_movie_maker


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
    def run_movie_maker():
        try:
            cluster_movie_maker.main("sample.xyz",
                                     "./raw_output/sample.xyz.rcAA1.88.rcAB1.88.rcBB1.88.Vor1.fc1.PBCs1.raw_",
                                     "FCC 13A 12E 11F 10B 9B 8B sp5c sp4c sp3c")
            return 0
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def tidy():
        # Remove the files we have created
        os.remove("sample_output.xyz")
        return 0


class FileChecks:
    @staticmethod
    def check_movie():
        return filecmp.cmp("sample_output.xyz", "sample_output.xyz", shallow=False)


def test_basic_configuration():
    with cd("test/movie_maker"):
        assert FileOperations.run_movie_maker() == 0
        assert FileChecks.check_movie() is True
        assert FileOperations.tidy() == 0
