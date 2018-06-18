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
# Each cluster should normally be detected as a single instance of itself however, some clusters are
# found more than once due to symmetry axes and some are not detected using different bond types.
simple_bond_clusters = {cluster: 1 for cluster in clusters}
simple_bond_clusters['sp3b'] = 4
simple_bond_clusters['sp4c'] = 3
simple_bond_clusters['7T_a'] = 2
simple_bond_clusters['7T_s'] = 3

# Whether clusters are detected by the Voronoi bond method depends on the fc parameter which changes the bond detection.
# The detection of clusters is defined for two different values of the fc parameter, 0.82 and 1.
# An fc parameter of 1 allows for a wider range of bond lengths including longer bonds so we call this voronoi_long.
# An fc parameter of 0.82 favours with a narrower range of lengths so we call this voronoi_short.
voronoi_short_clusters = {cluster: 1 for cluster in clusters}
voronoi_short_clusters['sp3b'] = 4
voronoi_short_clusters['sp4c'] = 3
voronoi_short_clusters['7T_a'] = 2
voronoi_short_clusters['7T_s'] = 3
voronoi_short_clusters['7K'] = 0
voronoi_short_clusters['8K'] = 0
voronoi_short_clusters['11A'] = 0

voronoi_long_clusters = {cluster: 1 for cluster in clusters}
voronoi_long_clusters['sp3b'] = 4
voronoi_long_clusters['sp4c'] = 3
voronoi_long_clusters['7T_a'] = 2
voronoi_long_clusters['7T_s'] = 3
voronoi_long_clusters['8A'] = 0
voronoi_long_clusters['9A'] = 0
voronoi_long_clusters['10A'] = 0
voronoi_long_clusters['11B'] = 0
voronoi_long_clusters['12A'] = 0
voronoi_long_clusters['FCC'] = 0
voronoi_long_clusters['HCP'] = 0