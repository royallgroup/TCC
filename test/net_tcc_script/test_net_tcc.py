import subprocess
from glob import glob
import filecmp
import os
import sys

sys.path += [os.path.abspath('../../lib/net_clusters')]
dir_path = os.path.dirname(os.path.realpath(__file__))

from net import net_cluster_calculation

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
    def run_net():
        try:
            net_cluster_calculation("./raw_output",  "(FCC, 13A, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c)")
            return 0
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def tidy():
        # Remove the files we have created
        os.remove("./raw_output/net_clusters.txt")
        return 0


class FileChecks:
    @staticmethod
    def check_net():
        return filecmp.cmp("sample_net.txt", "./raw_output/net_clusters.txt", shallow=False)


def test_basic_configuration():
    with cd("./basic_configuration"):
        assert FileOperations.run_net() == 0
        assert FileChecks.check_net() is True
        assert FileOperations.tidy() == 0
