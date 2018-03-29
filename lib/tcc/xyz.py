#!/usr/bin/env python3
import sys, os, io, pandas
import numpy as np

def write(xyz, out=sys.stdout, comment='', atoms='A'):
    if type(out) is str: out = open(out, 'w')
    if type(xyz) is not list: xyz = [xyz]

    for x in xyz:
        out.write(str(len(x)) + '\n')
        out.write(comment + '\n')
        if type(atoms) is str:
            for i in range(len(x)):
                out.write(atoms + ' ' + ' '.join(map(str, x[i,:])) + '\n')
        else:
            for i in range(len(x)):
                out.write(atoms[i] + ' ' + ' '.join(map(str, x[i,:])) + '\n')

def read(stream, read_atoms=False):
    if type(stream) is str: stream = open(stream)

    # Read the header - check for EOF.
    line = stream.readline()
    if not line: return None
    # Process the rest of the header rows.
    number_of_atoms = int(line)
    assert number_of_atoms is not None and number_of_atoms > 0
    comment = stream.readline()

    # Use pandas to read the main table.
    c = io.StringIO()
    for i in range(number_of_atoms): c.write(stream.readline())
    c.seek(0)
    table = pandas.read_table(c, sep='\s+', names=('atom','x','y','z'), nrows=number_of_atoms)

    coordinates = table[['x','y','z']].values.copy('c').astype(np.longdouble)

    if read_atoms: return table['atom'].tolist(), coordinates
    else: return coordinates

def multi_read(path, max_frames=None, read_atoms=False):
    with open(path) as f:
        frames = 0
        while True:
            configuration = read(f, read_atoms)
            if configuration is None: break

            yield configuration
            frames += 1
            if max_frames is not None and frames is max_frames: break

def trajectory(path, max_frames=None, read_atoms=False):
    return list(multi_read(path, max_frames, read_atoms))
