"""Information on clusters used by the TCC."""

# Commonly used names for structures and their equivalents in the TCC naming conventions.
cluster_list = ["sp3a", "sp3b", "sp3c", "sp4a", "sp4b", "sp4c", "sp5a", "sp5b", "sp5c",
                "6A", "6Z", "7K", "7T_a", "7T_s", "8A", "8B", "8K", "9A", "9B", "9K", "10A", "10B",
                "10K", "10W", "11A", "11B", "11C", "11E", "11F", "11W", "12A", "12B", "12D",
                "12E", "12K", "13A", "13B", "13K", "FCC", "HCP", "BCC_9"]

# Clusters are defined as the minimum energy structure for a particular potential.
# Each cluster should normally be detected as a single instance of itself however, some clusters are
# found more than once due to symmetry axes and some are not detected using different bond types.

# Whether clusters are detected by the Voronoi bond method depends on the fc parameter which changes the bond detection.
# Detection of clusters is defined for Voronoi bond detection with two different values of the fc parameter, 0.82 and 1.
# An fc parameter of 1 allows for a wider range of bond lengths including longer bonds so we call this voronoi_long.
# An fc parameter of 0.82 favours with a narrower range of lengths so we call this voronoi_short.
# The performance of the TCC is not defined if using a simple bond cutoff.
voronoi_short_clusters = {cluster: 1 for cluster in cluster_list}
voronoi_short_clusters['sp3b'] = 4      # Detected 4 times due to rotational symmetry
voronoi_short_clusters['sp4c'] = 3      # Detected 3 times due to rotational symmetry (the unique cluster is the 6A)
voronoi_short_clusters['7T_a'] = 2      # Detected twice due to symmetry
voronoi_short_clusters['7T_s'] = 3      # Detected 3 times due to symmetry
voronoi_short_clusters['7K'] = 0        # Not detected due to lack of long bonds
voronoi_short_clusters['8K'] = 0        # Not detected due to lack of long bonds

voronoi_long_clusters = {cluster: 1 for cluster in cluster_list}
voronoi_long_clusters['sp3b'] = 4       # Detected 4 times due to rotational symmetry
voronoi_long_clusters['sp4c'] = 0       # Not detected due to preference for sp3 rings over sp4 rings
voronoi_long_clusters['6A'] = 0         # Not detected due to preference for sp3 rings over sp4 rings
voronoi_long_clusters['7K'] = 3         # Detected 3 times due to deformed sp3 clusters detected when cluster is in isolation (not bulk)
voronoi_long_clusters['7T_a'] = 2       # Detected twice due to symmetry
voronoi_long_clusters['7T_s'] = 3       # Detected 3 times due to symmetry
voronoi_long_clusters['8A'] = 0         # Not detected due to interfering long bonds
voronoi_long_clusters['9A'] = 0         # Not detected due to interfering long bonds
voronoi_long_clusters['10A'] = 0        # Not detected due to interfering long bonds
voronoi_long_clusters['11B'] = 0        # Not detected due to interfering long bonds
voronoi_long_clusters['12A'] = 0        # Not detected due to interfering long bonds
voronoi_long_clusters['HCP'] = 0        # Not detected due to overestimation of number of neighbours
