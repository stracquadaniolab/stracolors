import sys
import matplotlib as mpl
from matplotlib import pyplot as plt


from palettable.cmocean.sequential import *
from palettable.colorbrewer.diverging import *
from palettable.colorbrewer.sequential import *
from palettable.colorbrewer.qualitative import *
from palettable.cartocolors.sequential import *
from palettable.cartocolors.qualitative import Prism_10
from palettable.cmocean.diverging import *

import seaborn as sns
mpl.use("TkAgg")

def reset_rc_params(backend= 'TkAgg',
                    font_family='Arial',
			        spine_left=False, spine_right=False, spine_bottom=False, spine_top=False, 
			        axes_facecolor='white', axes_edgecolor='#264258',
                    grid_alpha=0.5, grid_color='#264258', grid_linestyle='dotted', grid_linewidth=0.5,
                    args_dict={}):
    """
    Function that sets the rc params for matplotlib

    :param backend: backed to use for matplotlib
    :param spine_XXX: set the appearance of each spine
    :param axes_facecolor: background color
    :param axes_edgecolor: borders colors
    :param args_dict: pass a dictionary for 'mpl.rcParams[key]= value'
    """

    # Set all default values.
    mpl.rcdefaults()
    # Force agg backend.
    plt.switch_backend(backend)
    mpl.rcParams['font.family'] = font_family

    #mpl.rcParams['font.family'] = 'serif'
    #mpl.rcParams['font.serif'] = ['Charter']

    mpl.rcParams['axes.spines.left'] = bool(spine_left)
    mpl.rcParams['axes.spines.right'] = bool(spine_right)
    mpl.rcParams['axes.spines.bottom'] = bool(spine_bottom)
    mpl.rcParams['axes.spines.top'] = bool(spine_top)

    mpl.rcParams['axes.facecolor']= str(axes_facecolor)
    mpl.rcParams['axes.edgecolor']= str(axes_edgecolor)


    mpl.rcParams['grid.alpha']= float(grid_alpha)
    mpl.rcParams['grid.color']= str(grid_color)
    mpl.rcParams['grid.linestyle']= str(grid_linestyle)
    mpl.rcParams['grid.linewidth']= float(grid_linewidth)      ## in points

    for k,val in args_dict.items():
        mpl.rcParams[k]= val      ## in points


    

