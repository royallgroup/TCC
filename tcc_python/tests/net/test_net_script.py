import filecmp
from tcc_python.net import net_cluster_calculation
import pathlib

TEST_DIR = pathlib.Path("tcc_python/tests/net")
RAW_DIR = (TEST_DIR / pathlib.Path("raw_output")).absolute()
OUTPUT_FILE = TEST_DIR / pathlib.Path("raw_output/net_clusters.txt")
SAMPLE_OUTPUT = TEST_DIR / pathlib.Path("sample_net.txt")


def test_basic_configuration():
    net_cluster_calculation(RAW_DIR,  "(FCC, 13A, 12E, 11F, 10B, 9B, 8B, sp5c, sp4c, sp3c)")
    assert filecmp.cmp(SAMPLE_OUTPUT, OUTPUT_FILE, shallow=False) is True
    OUTPUT_FILE.unlink()
