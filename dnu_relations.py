import numpy as np

""" return the coefficient and power based on a keyword

Old code:
    # if Dnu_relation == 'RGB':
    #     dnu_est = 0.3*(numax_est)**0.75
    # elif Dnu_relation == 'RHB':
    #     dnu_est = 0.3*(numax_est)**0.86
    # elif Dnu_relation == 'AGB':
    #     dnu_est = 0.3*(numax_est)**0.77

    This is designed to be edited, if you have you own dnu-numax relation
"""

def Find_Dnu_relations(keyword):

    Dnu_relations = {'GC_RGB': [0.3,0.75],
                     'GC_RHB': [0.3, 0.86],
                     'GC_AGB': [0.3, 0.77],
                     'Yu18': [0.267, 0.764],
                     }

    return Dnu_relations[keyword]
