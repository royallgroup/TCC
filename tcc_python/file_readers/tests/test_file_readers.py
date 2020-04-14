from tcc_python.file_readers import snapshot, atom, dynamo, xyz
import pytest
import numpy as np
import math

FILE_DIR = "tcc_python/file_readers/tests/test_files"


class TestXYZ:
    @staticmethod
    def test_empty_file():
        import os
        cwd = os.getcwd()
        print(cwd)
        assert list(xyz.read(f"{FILE_DIR}/empty_file.xyz")) == []

    @staticmethod
    def test_partial_frame():
        with pytest.raises(snapshot.SnapshotIncompleteError):
            list(xyz.read(f"{FILE_DIR}/partial_frame.xyz"))

    @staticmethod
    def test_partial_multiple_frames():
        with pytest.raises(snapshot.SnapshotIncompleteError):
            list(xyz.read(f"{FILE_DIR}/partial_multiple_frames.xyz"))


class TestAtom:
    @staticmethod
    def test_empty_file():
        assert list(atom.read(f"{FILE_DIR}/empty_file.xyz")) == []

    @staticmethod
    def test_sample_file():
        data = list(atom.read(f"{FILE_DIR}/sample.atom"))[0]
        assert np.array_equal(data.box, np.array([[-1, 1], [-1, 1], [-1, 1]]))
        assert data.dimensionality == 3
        assert data.num_particles == 8
        assert np.array_equal(data.particle_coordinates[0], np.array([-1, -1, -1]))
        assert np.array_equal(data.species, np.array([1, 2, 1, 2, 1, 2, 1, 2]))


class TestDynamo:
    @staticmethod
    def test_empty_file():
        assert list(dynamo.read(f"{FILE_DIR}/empty_file.xyz")) == []

    @staticmethod
    def test_dynamo_snapshot():
        data = list(dynamo.read(f"{FILE_DIR}/sample_dynamo_snapshot.xml"))[0]
        assert np.allclose(data.box, np.array([[0, 8.71188], [0, 8.71188], [0, 8.71188]]), atol=0.00001)
        assert math.isclose(data.density, 2.075)
        assert data.diameters.shape == (1372,)
        assert data.dimensionality == 3
        assert data.num_particles == 1372
        assert data.particle_coordinates.shape == (1372, 3)
        assert data.species.shape == (1372,)
        assert math.isclose(data.volume, 661.204, rel_tol=0.001)
        assert math.isclose(data.volume_fraction, 0.5746, rel_tol=0.0001)

    @staticmethod
    def test_dynamo_config():
        data = list(dynamo.read(f"{FILE_DIR}/sample_dynamo_config.end.xml"))[0]
        assert np.allclose(data.box, np.array([[0, 9.68024], [0, 9.68024], [0, 9.68024]]), atol=0.00001)
        assert math.isclose(data.density, 1.5125)
        assert data.diameters.shape == (1372,)
        assert data.dimensionality == 3
        assert data.num_particles == 1372
        assert data.particle_coordinates.shape == (1372, 3)
        assert data.species.shape == (1372,)
        assert math.isclose(data.volume, 907.107, rel_tol=0.001)
        assert math.isclose(data.volume_fraction, 0.58617, rel_tol=0.0001)
