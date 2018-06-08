from python_scripts.file_readers import snapshot, atom, dynamo, xyz
import pytest


class TestXYZ:
    @staticmethod
    def test_empty_file():
        with pytest.raises(snapshot.NoSnapshotError):
            xyz.read("test/file_readers/test_files/empty_file.xyz")

    @staticmethod
    def test_partial_frame():
        with pytest.raises(snapshot.SnapshotIncompleteError):
            xyz.read("test/file_readers/test_files/partial_frame.xyz")

    @staticmethod
    def test_partial_multiple_frames():
        with pytest.raises(snapshot.SnapshotIncompleteError):
            xyz.read("test/file_readers/test_files/partial_multiple_frames.xyz")


class TestAtom:
    @staticmethod
    def test_empty_file():
        with pytest.raises(snapshot.NoSnapshotError):
            atom.read("test/file_readers/test_files/empty_file.xyz")


class TestDynamo:
    @staticmethod
    def test_empty_file():
        with pytest.raises(snapshot.NoSnapshotError):
            dynamo.read("test/file_readers/test_files/empty_file.xyz")
