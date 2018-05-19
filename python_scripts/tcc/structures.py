"""Information on clusters used by the TCC."""

# Commonly used names for structures and their equivalents in the TCC naming conventions.
cluster_index = {
    'sp3':  'sp3a',
    'sp3a': 'sp3a',
    'sp3b': 'sp3b',
    'sp3c': 'sp3c',
    '5A':   'sp3c',
    'sp4':  'sp4a',
    'sp4a': 'sp4a',
    'sp4b': 'sp4b',
    'sp4c': 'sp4c',
    '6A':   'sp4c',
    '6Z':   '6Z',
    'sp5':  'sp5a',
    'sp5a': 'sp5a',
    'sp5b': 'sp5b',
    'sp5c': 'sp5c',
    '7A':   'sp5c',
    '7T_a':  '7T_a',
    '7T_s':  '7T_s',
    '7K':   '7K',
    '8A':   '8A',
    '8B':   '8B',
    '8K':   '8K',
    '9A':   '9A',
    '9B':   '9B',
    '9K':   '9K',
    '10A':  '10A',
    '10B':  '10B',
    '10K':  '10K',
    '10W':  '10W',
    '11A':  '11A',
    '11B':  '11B',
    '11C':  '11C',
    '11D':  '11C',
    '11E':  '11E',
    '11F':  '11F',
    '11W':  '11W',
    '12A':  '12A',
    '12B':  '12B',
    '12C':  '12B',
    '12D':  '12D',
    '12E':  '12E',
    '12K':  '12K',
    '13A':  '13A',
    '13B':  '13B',
    '13K':  '13K',
    'FCC':  'FCC',
    'HCP':  'HCP',
    'BCC_9': 'BCC_9'
}

clusters = sorted(list(set(cluster_index.values())))
symmetry_number = {cluster: 1 for cluster in clusters}
symmetry_number['sp3b'] = 4
symmetry_number['sp4c'] = 3

# Components of each cluster in terms of elementary structures (and itself) for e.g. unit testing.
# Each cluster 'contains' at least one of itself by virtue of equality, though may find more it has symmetry axes.
composition = {cluster: {cluster: symmetry_number[cluster]} for cluster in clusters}
#for cluster in clusters: composition[cluster][cluster] = 1

# A square pyramid: sp4b
composition['sp4b']['sp3a'] = 4
# A Pentagonal pyramid: sp5b
composition['sp5b']['sp3a'] = 5

# The double tetrahedron: sp3c
composition['sp3c']['sp3b'] = 6

# The octahedron: sp4c
composition['sp4c']['sp3a'] = 8

# The pentagonal ring cluster: sp5c (or 7A)
composition['sp5c']['sp3b'] = 10
composition['sp5c']['sp3c'] = 5
composition['sp5c']['6Z'] = 5
composition['sp5c']['7K'] = 5

# The tritetrahedron: 6Z
composition['6Z']['sp3b'] = 8
composition['6Z']['sp3c'] = 2

# More complex polytetrahedra: 7T_a and 7T_s (4 tetrahedra)
composition['7T_a']['sp3b'] = 10
composition['7T_a']['sp3c'] = 3
composition['7T_a']['6Z'] = 2
composition['7T_a']['7T_a'] = 2
composition['7T_s']['sp3b'] = 10
composition['7T_s']['sp3c'] = 3
composition['7T_s']['6Z'] = 3
composition['7T_s']['7T_s'] = 3

# Crystal structures
composition['BCC_9']['sp3a'] = 12
composition['BCC_9']['sp4b'] = 6
