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

def call_palette(palette_type,
                number_of_colors=9,
                darkness='normal',
                shade='red',
                reverse= False
                ):

    """
    Function that returns the desired palette
    according to the specifications 

    :param palette_type: specify one between divergent/sequential/qualitative
    :param number_of_colors: specify the number of colors for the palette
    :param darkness: specify one between dark/normal/light, only for sequential
    :param shade: if the palette is sequential chooses one between red and blue, only for sequential
    :param reverse: choose to reverse the palette

    :return: a palette as a list of colors
    """
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
    
    elif palette_type == 'qualitative':
        palette=Prism_10.mpl_colors[0:number_of_colors:]
    
    else:
        print('ERROR: choose one between divergent/sequential/qualitative')
        sys.exit()
    
    if reverse:
        palette.reverse()
    
    return(palette)

def call_palette_binary(shade = 'color',
                        reverse= False):
    '''
    This function returns a binary palette, 
    the default order is dark to light (non-significant/significant).

    :param shade: choose one between color and bw
    :param reverse: choose the order (default is dark to light)

    :return: list of two colors
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

def call_palette_paired( number_of_classes = 5,
                        reverse = False):
    '''
    This function returns a palette of two shades for each class

    :param number_of_classes: number of hues needed, two shades each are returned, max is 5
    :param reverse: choose the order (default is dark to light)

    :return: list of coupled colors
    
    '''

    palette = ['#2166ac','#92c5de', '#b2182b', '#f4a582', '#1b7837','#a6dba0', '#762a83', '#c2a5cf' , '#4d4d4d','#bababa' ]  
    if reverse:
        palette=palette.reverse()
    noc=int(number_of_classes)
    return palette[0:noc*2]

def call_palette_triple( number_of_classes=5,
                        reverse= False):

    '''
    This function returns a palette of three shades for each class

    :param number_of_classes: number of hues needed, 3 shades each are returned, max is 5
    :param reverse: choose the order (default is dark to light)

    :return: list of color triplets
    
    '''

    palette = ['#2166ac','#67a9cf', '#d1e5f0', '#b2182b', '#ef8a62', '#fddbc7', 
                '#1b7837','#7fbf7b','#d9f0d3', '#762a83', '#af8dc3', '#e7d4e8' , '#4d4d4d','#999999', '#e0e0e0' ]  

    if reverse:
        palette=palette.reverse()

    noc=int(number_of_classes)
    return palette[0:noc*3]

def uoe_colors():

    '''
    Returns blue, red, white colors 
    of The Univerity of Edinburgh
    '''

    blue = '#00325F'
    red = '#C10043'
    white = '#FFFFFF'

    return [blue, red, white]
