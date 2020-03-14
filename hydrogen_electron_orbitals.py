import time

#start time recording
start = time.time()

import datetime
import mpmath
import math
from math import sqrt,sin,cos,tan
import mathutils
from itertools import chain
import bpy, bmesh
import numpy as np
import sympy
from sympy import symbols,I,latex,pi,diff,idiff #"I" is sympy's imaginary number
from sympy.utilities.lambdify import lambdastr
from sympy import factorial as fac
from sympy.functions import Abs,sqrt,exp,cos,sin
from sympy import re, im, simplify
import warnings # in order to suppress divide_by_zero warnings...

#display the latex representation of a symbolic variable by default.
from sympy import init_printing 
init_printing(use_unicode=True)

a_0,z,r,ib=symbols("a_0,z,r,ib")
n,m,l=symbols("n,m,l",integer=True)
int_m=symbols("int_m",integer=True)
theta,phi = symbols("\\theta,\\phi",real=True)


#The variables will used with lambdify...
angle_theta, angle_phi, radius = symbols("angle_theta,angle_phi,radius",real=True)

print("numpy version:   %s"%np.__version__)
print("mpmath version:   %s"%mpmath.__version__)
print("sympy version:   %s"%sympy.__version__)
print("Python version: %s"%bpy.app.version_string)

def P_l(l,theta): #valid for l greater than equal to zero
    """Legendre polynomial"""
    if l>=0:
#        eq=diff((cos(theta)**2-1)**l,cos(theta),l)
        eq=(cos(theta)**2-1)**l
        resultdiff=diff(eq.subs(cos(theta),ib),ib,l)
        eq=resultdiff.subs(ib,cos(theta))
    else:
        print("l must be an integer equal to 0 or greater")
        raise ValueError
    return 1/(2**l*fac(l))*eq

def P_l_m(m,l,theta):
    """Legendre polynomial"""
    #    eq = diff(P_l(l,theta),cos(theta),Abs(m))
    eq = P_l(l,theta)
    resultdiff=diff(eq.subs(cos(theta),ib),ib,Abs(m))
    eq=resultdiff.subs(ib,cos(theta))
    result = sin(theta)**Abs(m)*eq 
    return result

def Y_l_m(l,m,phi,theta):
    """Spherical harmonics"""
    eq = P_l_m(m,l,theta)
    if m>0:
        pe=re(exp(I*m*phi))*sqrt(2)
    elif m<0:
        pe=im(exp(I*m*phi))*sqrt(2)
    elif m==0:
        pe=1
    return abs(sqrt(((2*l+1)*fac(l-Abs(m)))/(4*pi*fac(l+Abs(m))))*pe*eq)

def L(l,n,rho):
    """Laguerre polynomial"""
    _L = 0.
    for i in range((n-l-1)+1): #using a loop to do the summation 
        _L += ((-i)**i*fac(n+l)**2.*rho**i)/(fac(i)*fac(n-l-1.-i)*\
                                          fac(2.*l+1.+i))
    return _L

def R(r,n,l,z=1.,a_0=1.):
    """Radial function"""
    rho = 2.*z*r/(n*a_0)
    _L = L(l,n,rho)
    _R = (2.*z/(n*a_0))**(3./2.)*sqrt(fac(n-l-1.)/\
         (2.*n*fac(n+l)**3.))*exp(-z/(n*a_0)*r)*rho**l*_L
    return _R

def Psi(r,n,l,m,phi,theta,z=1,a_0=1):
    """Wavefunction"""
    _Y = Y_l_m(l,m,phi,theta)
    _R = R(r,n,l)
    return _R*_Y

def P(r,n,l,m,phi,theta):
    """Returns the symbolic equation probability of the location of an electron"""
    return Psi(r,n,l,m,phi,theta)**2*r**2

def main(Prob,r_fun,phi_fun,theta_fun,n_o_c,box_size,res,isostep):
    # define a 3D scalarfield (the function which defines the shape of the isosurface)
    def scalarfield(pos):
        x,y,z=pos[0],pos[1],pos[2]
        w = Prob(r_fun(x,y,z),phi_fun(x,y,z),theta_fun(x,y,z)) * 1e2
        return w
#first point defining the gridbox of the MC-algorithm
    p0 = (-box_size,-box_size,-box_size)
#second point defining the gridbox of the MC-algorithm   
    p1 = (box_size,box_size,box_size)
#resolution in x,y,z direction of the grid (10x10x10 means 1000 cubes)      
    resolution = (res,res,res)           
#create for each isostep an isosurface starting from the outside to inside (low to high probability)      
    isosurface(p0,p1,resolution,isostep,scalarfield,n_o_c)


#
#
#

edgetable=(0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
            0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
            0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
            0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
            0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
            0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
            0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
            0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
            0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
            0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
            0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
            0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
            0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
            0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
            0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
            0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
            0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
            0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
            0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
            0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
            0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
            0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
            0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
            0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
            0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
            0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
            0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
            0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
            0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
            0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
            0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
            0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0)
tritable = [[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
        [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
        [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
        [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
        [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
        [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
        [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
        [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
        [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
        [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
        [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
        [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
        [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
        [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
        [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
        [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
        [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
        [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
        [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
        [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
        [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
        [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
        [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
        [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
        [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
        [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
        [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
        [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
        [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
        [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
        [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
        [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
        [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
        [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
        [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
        [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
        [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
        [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
        [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
        [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
        [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
        [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
        [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
        [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
        [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
        [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
        [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
        [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
        [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
        [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
        [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
        [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
        [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
        [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
        [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
        [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
        [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
        [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
        [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
        [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
        [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
        [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
        [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
        [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
        [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
        [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
        [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
        [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
        [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
        [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
        [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
        [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
        [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
        [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
        [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
        [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
        [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
        [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
        [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
        [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
        [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
        [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
        [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
        [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
        [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
        [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
        [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
        [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
        [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
        [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
        [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
        [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
        [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
        [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
        [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
        [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
        [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
        [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
        [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
        [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
        [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
        [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
        [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
        [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
        [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
        [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
        [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
        [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
        [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
        [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
        [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
        [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
        [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
        [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
        [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
        [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
        [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
        [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
        [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
        [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
        [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
        [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
        [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
        [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
        [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
        [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
        [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
        [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
]


def polygonise(cornervalues, isolevel, x1, y1, z1, x2, y2, z2):


    #   Determine the index into the edge table which
    #   tells us which vertices are inside of the surface
    cubeindex = 0
    if cornervalues[0] < isolevel: cubeindex = cubeindex | 1
    if cornervalues[1] < isolevel: cubeindex = cubeindex | 2
    if cornervalues[2] < isolevel: cubeindex = cubeindex | 4
    if cornervalues[3] < isolevel: cubeindex = cubeindex | 8
    if cornervalues[4] < isolevel: cubeindex = cubeindex | 16
    if cornervalues[5] < isolevel: cubeindex = cubeindex | 32
    if cornervalues[6] < isolevel: cubeindex = cubeindex | 64
    if cornervalues[7] < isolevel: cubeindex = cubeindex | 128

    # Cube is entirely in/out of the surface
    if edgetable[cubeindex] == 0: return []

    vertlist=[[]]*12
    # Find the vertices where the surface intersects the cube
    if (edgetable[cubeindex] & 1):    vertlist[0]  = vertexinterp(isolevel,[x1,y1,z1],[x1,y2,z1],cornervalues[0],cornervalues[1])
    if (edgetable[cubeindex] & 2):    vertlist[1]  = vertexinterp(isolevel,[x1,y2,z1],[x2,y2,z1],cornervalues[1],cornervalues[2])
    if (edgetable[cubeindex] & 4):    vertlist[2]  = vertexinterp(isolevel,[x2,y2,z1],[x2,y1,z1],cornervalues[2],cornervalues[3])
    if (edgetable[cubeindex] & 8):    vertlist[3]  = vertexinterp(isolevel,[x2,y1,z1],[x1,y1,z1],cornervalues[3],cornervalues[0])
    if (edgetable[cubeindex] & 16):   vertlist[4]  = vertexinterp(isolevel,[x1,y1,z2],[x1,y2,z2],cornervalues[4],cornervalues[5])
    if (edgetable[cubeindex] & 32):   vertlist[5]  = vertexinterp(isolevel,[x1,y2,z2],[x2,y2,z2],cornervalues[5],cornervalues[6])
    if (edgetable[cubeindex] & 64):   vertlist[6]  = vertexinterp(isolevel,[x2,y2,z2],[x2,y1,z2],cornervalues[6],cornervalues[7])
    if (edgetable[cubeindex] & 128):  vertlist[7]  = vertexinterp(isolevel,[x2,y1,z2],[x1,y1,z2],cornervalues[7],cornervalues[4])
    if (edgetable[cubeindex] & 256):  vertlist[8]  = vertexinterp(isolevel,[x1,y1,z1],[x1,y1,z2],cornervalues[0],cornervalues[4])
    if (edgetable[cubeindex] & 512):  vertlist[9]  = vertexinterp(isolevel,[x1,y2,z1],[x1,y2,z2],cornervalues[1],cornervalues[5])
    if (edgetable[cubeindex] & 1024): vertlist[10] = vertexinterp(isolevel,[x2,y2,z1],[x2,y2,z2],cornervalues[2],cornervalues[6])
    if (edgetable[cubeindex] & 2048): vertlist[11] = vertexinterp(isolevel,[x2,y1,z1],[x2,y1,z2],cornervalues[3],cornervalues[7])

    #Create the triangle
    triangles = []
    i=0
    while tritable[cubeindex][i] != -1:
        triangles.append([vertlist[tritable[cubeindex][i  ]],
                                   vertlist[tritable[cubeindex][i+1]],
                                   vertlist[tritable[cubeindex][i+2]]])
        i+=3

    return triangles

def vertexinterp(isolevel,p1,p2,valp1,valp2):
   if (ABS(isolevel-valp1) < 0.00001):
      return p1
   if (ABS(isolevel-valp2) < 0.00001):
      return p2
   if (ABS(valp1-valp2) < 0.00001):
      return p1
   mu = (isolevel - valp1) / (valp2 - valp1);
   x = p1[0] + mu * (p2[0] - p1[0]);
   y = p1[1] + mu * (p2[1] - p1[1]);
   z = p1[2] + mu * (p2[2] - p1[2]);

   return x,y,z

def create_mesh_for(objname,verts,faces):
    me = bpy.data.meshes.new(objname)  # create a new mesh
    me.from_pydata(verts,[],faces)
    me.update()      # update the mesh with the new data


    bm = bmesh.new()
    bm.from_mesh(me)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.01)
    bm.to_mesh(me)

    ob = bpy.data.objects.new(objname,me) # create a new object
    ob.data = me          # link the mesh data to the object
    return ob

def creategeometry(verts):
    faces=[]
    faceoffset=0
    for ver in verts:
        if len(ver)==4:
            faces.append((faceoffset+0,faceoffset+1,faceoffset+2,faceoffset+3))
            faceoffset+=4
        elif len(ver)==3:
            faces.append((faceoffset+0,faceoffset+1,faceoffset+2))
            faceoffset+=3
    return list(chain.from_iterable(verts)),faces

def make_object_in_scene(verts, scene,contour):
    verts,faces=creategeometry(verts)
    object_name = "orb_n=" + str(n) + "_" + " " + "l=" + str(l) + " " + "_" + "m=" + str(m)
    block=create_mesh_for(object_name,verts,faces)

    bpy.context.collection.objects.link(block)
    selectobj(block)

    return block

def selectobj(obj):
    for o2 in bpy.context.scene.objects:
        if o2==obj: 
            o2.select_set(state=True)

def arange(start, stop, step):
     r = start
     while r < stop:
        yield r
        r += step

def cellloop(p0,p1,r):
    for z in arange(p0[2],p1[2],r[2]):
     for y in arange(p0[1],p1[1],r[1]):
      for x in arange(p0[0],p1[0],r[0]):
        yield x,y,z

def cornerloop(x,y,z):
    for cz in (0,z):
        for cy,cx in zip((0,y,y,0),(0,0,x,x)):
             yield cx,cy,cz

def isosurface(p0,p1,resolution,isostep,isofunc,n_o_c):
    r=[(x1-x0)/sw for x0,x1,sw in zip(p0,p1,resolution)]

    triangles=[]
    z_a = p0[2]
    z_plane_a = [ [ isofunc([x,y,z_a]) for y in arange(p0[1], p1[1], r[1]) ] for x in arange(p0[0], p1[0], r[0])]
    
#    print("z_plane_a = ")
#    print(z_plane_a)
#    print(len(z_plane_a[0]))
#    print(" ")
    
    c_loop_1 = list( cornerloop(1,1,1) )

    cornervalues = [0]*8
# for each z plane value do:    
    for z in arange(p0[2], p1[2], r[2]):
        z2 = z + r[2]
        z_plane_b = [ [ isofunc([x,y, z2]) for y in arange(p0[1], p1[1], r[1])] for x in arange(p0[0], p1[0], r[0])]
#for each y plane do:
        for yi in range(len(z_plane_a[0])-1):
            y = p0[1]+yi*r[1]
            y2 = y + r[1]
#for each x plane do:
            for xi in range(len(z_plane_a)-1):
                x = p0[0]+xi*r[0]
                x2 = x + r[0]
                if True:
                    cornervalues = [
                        z_plane_a[xi][yi],
                        z_plane_a[xi][yi+1],
                        z_plane_a[xi+1][yi+1],
                        z_plane_a[xi+1][yi],
                        z_plane_b[xi][yi],
                        z_plane_b[xi][yi+1],
                        z_plane_b[xi+1][yi+1],
                        z_plane_b[xi+1][yi],
                    ]
                else:
                    cornervalues = [ (z_plane_a if cz==0 else z_plane_b)[xi+cx][yi+cy] for cx,cy,cz in c_loop_1]
                for contour in range(1, n_o_c + 1, 1):    
                    isolevel = (contour) * isostep
                    triangles.extend(polygonise(cornervalues, isolevel, x,y,z, x2, y2, z2))

        z_plane_a = z_plane_b

    return make_object_in_scene(triangles, bpy.context.scene, n_o_c)

def find_3dview_space():
    # Find 3D_View window and its scren space
    area = None
    for a in bpy.data.window_managers[0].windows[0].screen.areas:
        if a.type == 'VIEW_3D':
            area = a
            break

    if area:
        space = area.spaces[0]
    else:
        space = bpy.context.space_data

    return space


def display_orbital(n,l,m,n_o_c,box_size,resolution,isostep):
    """Diplays a 3D view of electron orbitals"""

    P_tex = "" #initialize the LaTex string of the probabilities
    
    #Validate the quantum numbers
    assert(n>=1), "n must be greater or equal to 1"   #validate the value of n
    assert(0<=l<=n-1), "l must be between 0 and n-1"  #validate the value of l
    assert(-l<=m<=l), "p must be between -l and l"    #validate the value of p
        
    #Determine the probability equation symbolically and convert
    #it to a string
    prob = lambdastr((radius,angle_phi,angle_theta), P(radius,n,l,m,angle_phi,angle_theta))
#    print(prob)
    #record the probability equation as a LaTex string
    P_eq = simplify(P(r,n,l,m,phi,theta))
    P_tex+="$$P ="+latex(P_eq)+"$$ \n\n " 
    
    if '(nan)' in prob: #Check for errors in the equation
        print("There is a problem with the probability function.")
        raise ValueError
    
    #Convert the finctions in the probability equation from the sympy  
    #library to the numpy library to allow for the use of matrix 
    #calculations
    prob = prob.replace('math.sin','np.sin') #convert to numpy
    prob = prob.replace('math.cos','np.cos') #convert to numpy
    prob = prob.replace('math.Abs','np.abs') #convert to numpy
    prob = prob.replace('math.pi','np.pi')   #convert to numpy
    prob = prob.replace('math.exp','np.exp') #convert to numpy
    
#    print("Sybolic Prob function: ")
#    print(prob)

    #convert the converted string to a callable functio
    Prob = eval(prob)
    
    #go and let the marching boxes do their thing and create the isosurface mesh
    main(Prob,r_fun,phi_fun,theta_fun,n_o_c,box_size,resolution,isostep)      
    
    return

def create_blender_objects(n_o_c, isostep,n,l,m):
            #box-size is based on quantum number n as the size of the final generated               
            #object changes with n. If you get divide by zero error ensure that list entry sizes are uneven
            box_size_list = [5,13,35,70,110,180,250]
            box_size = box_size_list[n-1]

            #mesh resolution in x,y & z direction of the grid (eg. 25 means 25x25x25 = 15,625 cubes)
            #If a resolution less than 150 is used the marching cubes algiritm has difficulty creating smooth                #meshes in the higher n values
            resolution = 150

            P_tex = "" #initialize the LaTex string of the probabilities

            #Create isosurface meshes for each isostep
            display_orbital(n,l,m,n_o_c,box_size,resolution,isostep)


            #add material to the generated isosurface object(s)
            bpy.ops.object.select_all(action='DESELECT')
            objectname = "orb_n=" + str(n) + "_" + " " + "l=" + str(l) + " " + "_" + "m=" + str(m)
            ob = bpy.data.objects[objectname]

            # Get material
            mat = bpy.data.materials.get("Iso 01")
            if mat is None:
                # create material
                mat = bpy.data.materials.new(name="Iso 01")

            # Assign it to object
            if ob.data.materials:
                # assign to 1st material slot
                ob.data.materials[0] = mat
            else:
                # no slots
                ob.data.materials.append(mat)

            #recalculate normals to outside
            ob.select_set(state=True)
            bpy.context.view_layer.objects.active = ob
            # go edit mode
            bpy.ops.object.mode_set(mode='EDIT')
            # select al faces
            bpy.ops.mesh.select_all(action='SELECT')
            # recalculate outside normals 
            bpy.ops.mesh.normals_make_consistent(inside=False)
            # go object mode again
            bpy.ops.object.editmode_toggle()

            #move object to new location based on n,l & m multiplied by offset
            offset=440
            bpy.context.object.location[0] = n*offset #x
            bpy.context.object.location[1] = -m*offset #y
            bpy.context.object.location[2] = l*offset #z
            
            bpy.ops.object.shade_smooth()

            print("orb_n=" + str(n) + "_" + " " + "l=" + str(l) + " " + "_" + "m=" + str(m) + " created")
            
#Recursivly transverse layer_collection for a particular name
def recurLayerCollection(layerColl, collName):
    found = None
    if (layerColl.name == collName):
        return layerColl
    for layer in layerColl.children:
        found = recurLayerCollection(layer, collName)
        if found:
            return found
        
##################################################################################################
# Start

np.seterr(divide='ignore', invalid='ignore')
#Needed in display_orbital
r_fun = lambda _x,_y,_z: (np.sqrt(_x**2+_y**2+_z**2))
theta_fun = lambda _x,_y,_z: (np.arccos(_z/r_fun(_x,_y,_z)))
phi_fun = lambda _x,_y,_z: (np.arctan(_y/_x)*(1+_z-_z))

vec=mathutils.Vector
ABS=abs

#Change the Active LayerCollection to 'Orbitals' (the n.l.m Blender onjects will live here)
layer_collection = bpy.context.view_layer.layer_collection
layerColl = recurLayerCollection(layer_collection, 'Orbitals')
bpy.context.view_layer.active_layer_collection = layerColl

#Delete previous generated Blender MESH objects
bpy.ops.object.select_all(action='DESELECT')
for ob in bpy.context.scene.objects:
    if ob.type == 'MESH' and ob.name.startswith("orb_"):
        #Select the object
        ob.select_set(state=True)
#Delete all objects selected above 
bpy.ops.object.delete()

#number of isosurfaces to build within the gridbox 
n_o_c = 1
#probability space between isosurfaces            
isostep = 0.1


#if what is 'single' only one n,l & m blender object will be created, if it is anything else multiple blender objects for n,l,m will be created
what='single'

if what == 'single':
    #Single n,l,m blender object will be created
    
    #n is  the principle quantum number and relates to the period the element is in, or the shell.
    n = 1
    #l is the angular momentum quantum number which defines the sub shell s, p, d, f, of which there are 
    #n subshells whose values are 0 <= l <= (n-1)
    l = 0
    #m is the magnetic quantum number which further subdivides the subshell into orbitals, of which                  #there are 2l + 1 orbitals whose values are -l <= m <= +l
    m = 0
    create_blender_objects(n_o_c, isostep,n,l,m)
else: 
    #multiple n,l,m blender objects will be created
    
    #n is  the principle quantum number and relates to the period the element is in, or the shell.
    for n in range(1,8):
     
        #l is the angular momentum quantum number which defines the sub shell s, p, d, f, of which there are 
        #n subshells whose values are 0 <= l <= (n-1)
        for l in range(0,n):
     
            #m is the magnetic quantum number which further subdivides the subshell into orbitals, of which                  #there are 2l + 1 orbitals whose values are -l <= m <= +l
            for m in range(-l,l+1):            
                create_blender_objects(n_o_c, isostep,n,l,m)
      

bpy.ops.object.select_all(action='DESELECT')

elapsed = time.time()-start
elapsed =round(elapsed)
conversion = datetime.timedelta(seconds=elapsed)
converted_time = str(conversion)
print("Elapsed Time %r"%converted_time)