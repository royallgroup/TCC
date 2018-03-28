"""Information on clusters used by the TCC."""

class Translation:
    """Translations between commonly used structure names to the
    internal names used by the TCC and its outputs (which also
    contains redundant point group information.
    """

    structure_to_tcc_id = {
        'sp3': 'sp3', 'sp3a': 'sp3a', 'sp3b': 'sp3b', '5A': '5A_D3h',
        'sp4': 'sp4', 'sp4a': 'sp4a', 'sp4b': 'sp4b', 'sp4c': 'sp4c',
        '6A': '6A_Oh', '6Z': '6Z_C2v', 'tetra': 'tetra',
        'sp5': 'sp5', 'sp5a': 'sp5a', 'sp5b': 'sp5b',
        '7A': '7A_D5h', '7K': '7K',
        '8A': '8A_D2d', '8B': '8B_Cs', '8K': '8K',
        '9A': '9A_D3h', '9B': '9B_C2v', '9K': '9K',
        '10A': '10A_D4d', '10B': '10B_C3v', '10K': '10K', '10W': '10W',
        '11A': '11A_D4d', '11B': '11B_C2v', '11C': '11CD', '11E': '11E_C2',
        '11F': '11F_C2v', '11W': '11W_Cs',
        '12A': '12A_C2v', '12B': '12BC', '12D': '12D_D2d',
        '12E': '12E_D3h', '12K': '12K',
        '13A': '13A_Ih', '13B': '13B_D5h', '13K': '13K',
        'FCC': 'FCC_m13', 'HCP': 'HCP_m13',
        'BCC_9': 'BCC_m9', 'BCC_15': 'BCC_m15'}

    tcc_id_to_structure = {value: key for key,value in structure_to_tcc_id.items()}
