import subprocess
from glob import glob
import shutil
import filecmp
import os


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
    def copy_net():
        # Copy the exectuable to the current directory
        try:
            shutil.copy("../../../tools/net_clusters/net.py", os.getcwd())
            return 0
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def run_net():
        try:
            return subprocess.call(["python", "net.py", "sample.xyz"])
        except Exception as e:
            print(e)
            return 1

    @staticmethod
    def tidy():
        # Remove the files we have created
        os.remove("net.py")
        for file in (glob("sample.xyz*")):
            os.remove(file)
        return 0


class FileChecks:
    @staticmethod
    def check_net():
        return filecmp.cmp("sample_net.txt", glob("sample.xyz*net.txt")[0], shallow=False)


def test_basic_configuration():
    with cd("./basic_configuration"):
        assert FileOperations.copy_net() == 0
        assert FileOperations.run_net() == 0
        assert FileChecks.check_net() is True
        assert FileOperations.tidy() == 0
        