import subprocess
from os import getcwd, remove
from glob import glob
import numpy
import shutil


def run_tcc():
    try:
        # Copy the exectuable to the current directory and run it
        shutil.copy(glob("../bin/tcc*")[0], getcwd())
        subprocess.call(glob("tcc*")[0], shell=True)
        # Load the sample and measured data
        sample_data = numpy.genfromtxt("sample.static_clust", skip_header=2, usecols=1, max_rows=43, dtype=int)
        result = numpy.genfromtxt(glob("sample.xyz.rc*")[0], skip_header=2, usecols=1, max_rows=43, dtype=int)
        # Check for differences
        delta_check = numpy.absolute(sample_data - result).max()
        # Remove the files we have created
        remove(glob("tcc*")[0])
        remove(glob("sample.xyz.rc*")[0])

        return delta_check

    except Exception as e:
        print(e)
        return 1


def test_tcc():
    assert run_tcc() == 0


test_tcc()
