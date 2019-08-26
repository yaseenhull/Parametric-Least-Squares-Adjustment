# -*- coding: utf-8 -*-
"""
Created on Tue May 16 10:24:51 2017

@author: Student
"""


import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

import sympy as sp
import numpy as np
import networkx as nx
import math as mt

def polar(xs,ys,dist,dirt):
    xn = xs + dist*sp.cos(dirt)
    yn = ys + dist*sp.sin(dirt)
    return xn,yn

def polar1(xs,ys,dist,dirt):
    xn = xs + dist*mt.cos(np.radians(dirt))
    yn = ys + dist*mt.sin(np.radians(dirt))
    return xn,yn

def distance(xp,yp,xn,yn):
    return sp.sqrt(((xp-xn)**2)+((yp-yn)**2))
    
def distance1(xp,yp,xn,yn):
    return np.sqrt(((xp-xn)**2)+((yp-yn)**2))

def direction(xp,yp,xn,yn,z):
    return sp.atan((yp-yn)/(xp-xn)) - z
    
def direction1(xp,yp,xn,yn,z):
    return np.arctan(np.radians((yp-yn)/(xp-xn))) - z
  
def edge_S(s,e):
    return G.get_edge_data(s,e)['Distance']
    
def edge_D(s,e):
    return G.get_edge_data(s,e)['Direction']

def edge_sS(s,e):
    return G.get_edge_data(s,e)['Sigma Dist']

def edge_sD(s,e):
    return G.get_edge_data(s,e)['Sigma Dir']

def get_X(dic,Point):
    return dic[Point][0]

def get_Y(dic,Point):
    return dic[Point][1]

def get_Xx(dic,Point):
    return dic[Point]['X']

def get_Yy(dic,Point):
    return dic[Point]['Y']   

    
def joinDir(X2,Y2,X1,Y1):
    dy = Y2 - Y1
    dx = X2 - X1
    
    directn = np.arctan2(dy, dx)
    
    if directn < 0:
        directn = directn + (2*np.pi)
        
    return np.degrees(directn)
    
def joinDist(X2,Y2,X1,Y1):
    dy = Y2 - Y1
    dx = X2 - X1 
    
    dist = np.sqrt(dy**2 + dx**2)
    
    return dist
    
def pdiffDx(Xi,Yi,Xj,Yj):
    return -(Yj - Yi)/((1 + (Yj - Yi)**2/(Xj - Xi)**2)*(Xj - Xi)**2)
    
def pdiffDy(Xi,Yi,Xj,Yj):
    return 1/((1 + (Yj - Yi)**2)/((Xj - Xi)**2)*(Xj - Xi))

def pdiffSx(Xa,Ya,Xb,Yb):
    return (Xb - Xa)/np.sqrt((Xb - Xa)**2 + (Yb - Ya)**2)
    
def pdiffSy(Xi,Yi,Xj,Yj):
    return (Yj - Yi)/np.sqrt((Xj - Xi)**2 + (Yj - Yi)**2)

