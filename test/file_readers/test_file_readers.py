from python_scripts.file_readers import snapshot, atom, dynamo, xyz
import pytest
import numpy as np


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

    @staticmethod
    def test_sample_file():
        data = atom.read("test/file_readers/test_files/sample.atom")
        assert np.array_equal(data.box, np.array([[-1, 1], [-1, 1], [-1, 1]]))
        assert data.dimensionality == 3
        assert data.num_particles == 8
        assert np.array_equal(data.particle_coordinates[0], np.array([-1, -1, -1]))
        assert np.array_equal(data.species, np.array([1, 2, 1, 2, 1, 2, 1, 2]))


class TestDynamo:
    @staticmethod
    def test_empty_file():
        with pytest.raises(snapshot.NoSnapshotError):
            dynamo.read("test/file_readers/test_files/empty_file.xyz")
