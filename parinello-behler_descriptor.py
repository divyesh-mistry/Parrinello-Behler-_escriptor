import numpy as np
from mpmath import *
import sympy as sp
from numpy import linalg as LA

# Load Data file which contain position of each atom within cutoff radious

pos = np.loadtxt("position.txt")
x=[]
y=[]
z=[]

for line in pos:
    x.append(line[0])
    y.append(line[1])
    z.append(line[2])
# position of selected atom
xc = 60.8857
yc = 77.0608
zc = 2.59301


#N =len(x)   # total number of atom inside the cutoff radious r=10
N=3
l=2


def position(x,y,z,xc,yc,zc,N):
    """Position function takes all atom position of all atoms about reference
    atom (xc,yc,zc) in terms of (x,y,z) co-ordinates and return position
    matrix.here theta(θ) and phi(φ) which represent colatitude (polar angle)
    and the longitude φ, or azimuth.

    Args:
        x ([float]): [position x]
        y ([float]): [position y]
        z ([float]): [position z]
        xc ([float]): [x position of reference atom xc ]
        yc ([float]): [y position of reference atom yc ]
        zc ([float]): [z position of reference atom zc ]
        N ([int]): [Number atoms]

    Returns:
        [list]: [position of each atom: r_ij,theta,phi]
    """
    r_ij = np.empty(shape=[N,N],dtype='object')
    r=np.empty(shape=[N,N],dtype='object')
    phi_ij= np.empty(shape=[N,N],dtype='object')
    theta_ij = np.empty(shape=[N,N],dtype='object')

    for i in range(N):
        for j in range(N):
            r_ij[i][j]=(x[j]-x[i]-xc),(y[j]-y[i]-yc),(z[j]-z[i]-zc)

            r [i][j]= LA.norm(r_ij[i][j])
            theta_ij[i][j] = float(asin(sqrt((x[j]-x[i]-xc)*(x[j]-x[i]-xc)+(y[j]-y[i]-yc)*(y[j]-y[i]-yc))/r[i][j]))
            phi_ij[i][j] = float(acot((y[j]-y[i]-yc)/(x[j]-x[i]-xc)))
    return (r_ij,r,theta_ij,phi_ij)
r_ij,r,theta_ij,phi_ij=position(x,y,z,xc,yc,zc, N)
