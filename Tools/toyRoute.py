import numpy as np
import math
import os
import matplotlib.pyplot as plt

def transition_x(l, L, R):
    # Get x~l on trainsition curve
    w = R*L
    x = l - l**5/(40*(w**2)) # 
    return x

def transition_y(l, L, R):
    # Get y~l on trainsition curve
    w = R*L
    y = l**3/(6*w) - l**7/(336*(w**3))
    return y

def transition_angle(l, L, R):
    # Get angle~l on trainsition curve
    w = R*L
    alpha = l**2/(2*w)
    return alpha

def ToyRouteInfor(L, R):
    ## Node coordinate, curvature, and direction of (C1-B2-A3-B4-C5)

    ToyRouteInfor = {}

    ## Length
    # L = [150, 100, 400, 100, 150]

    ## Junction
    # end node of A3 / start node of B4
    x4, y4 = L[2]/2, 0
    # end node of B4 / start node of C5
    x5, y5 = x4 + transition_x(L[3], L[3], R), y4 + transition_y(L[3], L[3], R)
    # center of the circle
    xo, yo = x5 - R*math.sin(transition_angle(L[3], L[3], R)), y5 + R*math.cos(transition_angle(L[3], L[3], R))
    # end node of C5
    x6, y6 = xo + R*math.sin(transition_angle(L[3], L[3], R) + L[4]/R), yo - R*math.cos(transition_angle(L[3], L[3], R) + L[4]/R)
    # start node of C1
    x1, y1 = -1*x6, -1*y6
    # end node of C1 / start node of B2
    x2, y2 = -1*x5, -1*y5
    # end node of B2 / start node of A3
    x3, y3 = -1*x4, -1*y4
    junctions = [[x1, y1], [x2, y2], [x3, y3], [x4, y4], [x5, y5], [x6, y6]]
    ## Node Coordinates (1 meter/node) and curvature
    # A3 
    section3_x = [i for i in np.arange(-L[2]/2, L[2]/2, 1)]
    section3_y = [0. for i in np.arange(len(section3_x))]
    curvature3 = [0. for i in range(len(section3_x))]
    # B4
    section4_x = [x4 + transition_x(i, L[3], R) for i in np.arange(0, L[3], 1)]
    section4_y = [0 + transition_y(i, L[3], R) for i in np.arange(0, L[3], 1)]
    curvature4 = [0 + i/R/L[3] for i in np.arange(0, L[3], 1.)]
    # C5
    section5_x = [xo + R*math.sin(transition_angle(L[3], L[3], R)+i/R) for i in np.arange(0, L[4]+1, 1)]
    section5_y = [yo - R*math.cos(transition_angle(L[3], L[3], R)+i/R) for i in np.arange(0, L[4]+1, 1)]
    curvature5 = [1/400 for i in range(len(section5_x))]
    # C1
    section1_x, section1_y = [-1*i for i in section5_x[::-1]], [-1*i for i in section5_y[::-1]]
    curvature1 = curvature5[::-1]
    # B2
    section2_x, section2_y = [-1*i for i in section4_x[::-1]], [-1*i for i in section4_y[::-1]]
    curvature2 = curvature4[::-1]

    nodes_x = section1_x + section2_x + section3_x + section4_x + section5_x
    nodes_y = section1_y + section2_y + section3_y + section4_y + section5_y
    curvature = curvature1 + curvature2 + curvature3 + curvature4 + curvature5

    ToyRouteInfor['junction'] = junctions
    ToyRouteInfor['nodex'] = nodes_x
    ToyRouteInfor['nodey'] = nodes_y
    ToyRouteInfor['curvature'] = curvature
    ToyRouteInfor['section'] = {'section1':{'x': section1_x, 'y': section1_y}, 
                                'section2':{'x': section2_x, 'y': section2_y},
                                'section3':{'x': section3_x, 'y': section3_y}, 
                                'section4':{'x': section4_x, 'y': section4_y},
                                'section5':{'x': section5_x, 'y': section5_y}}

    return ToyRouteInfor
