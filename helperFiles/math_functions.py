from constants import *

'''
colorToProj:  color -> projected
dataToProj:   data -> projected
colorCoord:   data, projected -> color
dataToColor:  data -> color
colorToData:  color -> data
dataCoord:    color, projected -> data
'''


def colorToProj(x, y):
    return ((150 * x) + float(mapCoord['cmap_proj_x0'])), ((-150 * y) + float(mapCoord['cmap_proj_y0']))


def colorCoord(x, y):
    return ((-float(mapCoord['cmap_proj_x0']) + x) / 150.0), (-(-float(mapCoord['cmap_proj_y0']) + y) / 150.0)