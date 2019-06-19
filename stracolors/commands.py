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
mlp.use("TkAgg")

def reset_rc_params(backend:'matplotlib backend'= 'TkAgg',
                    font_family='Arial'):
    """
    Function that sets the rc params for matplotlib
    """
    # Set all default values.
    mpl.rcdefaults()
    # Force agg backend.
    plt.switch_backend(backend)
    #matplotlib.use("TkAgg")
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    mpl.rcParams['font.family'] = font_family

    #mpl.rcParams['font.family'] = 'serif'
    #mpl.rcParams['font.serif'] = ['Charter']

    mpl.rcParams['axes.spines.left'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.bottom'] = False
    mpl.rcParams['axes.spines.top'] = False

    mpl.rcParams['axes.facecolor']= 'white'
    mpl.rcParams['axes.edgecolor']= '#264258'


    mpl.rcParams['grid.alpha']=0.5
    mpl.rcParams['grid.color']='#264258' ## grid color
    mpl.rcParams['grid.linestyle']='dotted'        ## solid
    mpl.rcParams['grid.linewidth']= 0.5      ## in points

    


def call_palette(palette_type:'specify one between divergent/sequential/qualitative',
                number_of_colors: 'specify the number of colors for the palette'=9,
                darkness: 'specify one between dark/normal/light, only for sequential'='normal',
                shade: 'if the palette is sequential chooses one between red and blue, only for sequential'='red',
                reverse= False
                ):
    '''
    Function that returns the desired palette
    according to the necessities specifies.
    '''
    matplotlib.use("TkAgg")

    if palette_type == 'divergent':
        if number_of_colors%2:
            nc=int((number_of_colors-1)/2)
            palette=RdBu_11.mpl_colors[4-nc:5+nc]
        else:
            nc=int((number_of_colors)/2)
            palette=RdBu_10.mpl_colors[5-nc:5+nc]

    elif palette_type == 'sequential':

        step=9//number_of_colors
        print(step)
        if shade == 'red':
            if darkness == 'normal':
                palette = OrRd_9.mpl_colors            
            elif darkness == 'dark':
                palette = Amp_9.mpl_colors
            elif darkness == 'light':
                palette = Reds_9.mpl_colors
            else:
                print('ERROR: choose one between red and blue')
                sys.exit()

        elif shade == 'blue':
            
            if darkness == 'normal':
                palette = PuBu_9.mpl_colors
            elif darkness == 'dark':
                palette = Tempo_9.mpl_colors
            elif darkness == 'light':
                palette = Blues_9.mpl_colors
                
            else:
                print('ERROR: choose one between red and blue')
                sys.exit()
        
        else:
            print('ERROR: choose one between red and blue')
            sys.exit()

        palette = palette[0:10:step] #PuBu_9
        print(palette)
    
    elif palette_type == 'qualitative':
        print('qualitative chosen')
        palette=Prism_10.mpl_colors[0:number_of_colors:]
    
    else:
        print('ERROR: choose one between divergent/sequential/qualitative')
        sys.exit()
    
    if reverse:
        palette=palette.reverse()
    
    return(palette)

def call_palette_binary(shade: 'choose one between color and bw' = 'color',
                        reverse: 'choose the order (default is dark to light)' = False):
    '''
    This function returns a binary palette, 
    the default order is dark to light (non-significant/significant).
    '''

    if shade == 'color':
        palette_binary = ['#0571b0','#ca0020']    
    elif shade == 'bw':
        palette_binary = ['#000000', '#FFFFFF' ]
    else:
        print('ERROR: choose one between color/bw')
        sys.exit()
    if reverse:
        palette_binary=palette_binary.reverse()
    return palette_binary



def uoe_colors():
    '''
    Returns blue, red, white colors 
    of The Univerity of Edinburgh
    '''
    blue = '#00325F'
    red = '#C10043'
    white = '#FFFFFF'

    return [blue, red, white]

def show_choices(output_folder: 'output_folder'):

    f,ax=plt.subplots(1)
    uoe=uoe_colors()
    sns.palplot(sns.color_palette(uoe))
    plt.savefig(output_folder+'palette_uoe.png', f='png')


    f,ax=plt.subplots(1)
    pal=call_palette_binary()
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_binary_color.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette_binary(shade='bw')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_binary_bw.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='qualitative')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_qualitative.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='divergent')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_divergent.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='divergent', number_of_colors=4)
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_divergent_4.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', shade='blue')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_blue.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', darkness='dark')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_red_dark.png', f='png')


    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', shade='blue',darkness='dark')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_blue_dark.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', darkness='light')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_red_light.png', f='png')


    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', shade='blue',darkness='light')
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_blue_light.png', f='png')

    f,ax=plt.subplots(1)
    pal=call_palette(palette_type='sequential', shade='blue',darkness='light', number_of_colors=3)
    sns.palplot(sns.color_palette(pal))
    plt.savefig(output_folder+'palette_sequential_blue_light_3.png', f='png')
