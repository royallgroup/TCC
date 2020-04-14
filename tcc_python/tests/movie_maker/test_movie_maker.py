import filecmp
from tcc_python import movie_maker
import pathlib

TEST_DIR = pathlib.Path("tcc_python/tests/movie_maker")
INPUT_FILE = TEST_DIR / pathlib.Path("sample.xyz")
OUTPUT_FILE = TEST_DIR / pathlib.Path("output.xyz")
CORRECT_OUTPUT = TEST_DIR / pathlib.Path("sample_output.xyz")


def test_basic_configuration():
    movie_maker.main(INPUT_FILE, f"{TEST_DIR}/raw_output/sample.xyz.rcAA1.88.rcAB1.88.rcBB1.88.Vor1.fc1.PBCs1.raw_",
                     "FCC 13A 12E 11F 10B 9B 8B sp5c sp4c sp3c", OUTPUT_FILE)
    assert filecmp.cmp(CORRECT_OUTPUT, OUTPUT_FILE, shallow=False) is True
    OUTPUT_FILE.unlink()
