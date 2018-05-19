import sys
sys.path.append('../../python_scripts/file_readers')

import snapshot
import atom
import dynamo
import xyz


class TestXYZ:
    @staticmethod
    def test_empty_file():
        xyz.XYZSnapshot("test_files/empty_file.xyz")


class TestAtom:
    @staticmethod
    def test_empty_file():
        atom.AtomSnapshot("test_files/empty_file.xyz")
    

class TestDynamo:
    @staticmethod
    def test_empty_file():
        dynamo.Snapshot("test_files/empty_file.xyz")


class TestSnapshot:
    @staticmethod
    def test_empty_file():
        snapshot.Snapshot("test_files/empty_file.xyz")
