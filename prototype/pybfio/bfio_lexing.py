#!/usr/bin/env python
"""Butterfly algorithm for evaluating Fourier Integral Operators.

"""

import math
import itertools
import numpy
import interpolation


__author__ = "Nicholas Maxwell, Laurent Demanet, Lexing Ying"
__copyright__ = "not sure"
__credits__ = ["Nicholas Maxwell", "Laurent Demanet", "Lexing Ying"]
__license__ = "MIT"
__version__ = "0.0.0"
__maintainer__ = "Nick Maxwell"
__email__ = "nicholas.maxwell@gmail.com    "
__status__ = "Prototype"


CONST_SQRT2_2 = 0.70710678118654757
CONST_2PI = 6.2831853071795865

CONST_2PI_J = 6.2831853071795865j
CONST_PI_SQRT2_J = 4.4428829381583661J

def log2(x):
    return (int)(math.floor(math.log(x,2)))

def range2(n,m=None):
    if m is None:
        m=n
    return itertools.product(range(n),range(m))

def polar(N, k):
    p1 = math.atan2(k[1],k[0])/CONST_2PI+0.5
    p0 = math.sqrt(k[0]**2 + k[1]**2)/(CONST_SQRT2_2*N)
    return numpy.array((p0, p1))










def bfio_eval(N, input_sources, input_phase, n_cheby, tree_length):

    # some pre-processing

    phase = lambda x,p: CONST_PI_SQRT2_J*N*p[0]*input_phase(x,(numpy.cos(CONST_2PI*p[1]),numpy.sin(CONST_2PI*p[1])))
    sources = [ (polar(N, (float(j1)-N/2, float(j2)-N/2)), input_sources[j1,j2]) for j1,j2 in range2(N) ]

    boxed_sources = {}
    for j1,j2 in range2(2**tree_length):
        a1,a2,b1,b2 = k_spacing*j1, k_spacing*j2, k_spacing*(j1+1), k_spacing*(j2+1)
        leaf=[]
        for (source_point, source) in sources:
            if a1 <= source_point[0] and source_point[0] < b1 and a2 <= source_point[1] and source_point[1] < b2:
                leaf.append( (source_point, source) )
        boxed_sources[j1,j2] = leaf

    sources = boxed_sources


    # preliminaries
    
    unit_cheby_grid_1d = interpolation.chebyshev_grid(n_cheby )*0.5+0.5
    unit_cheby_grid_2d = numpy.zeros((n_cheby,n_cheby))
    for i,j in range2(n_cheby):
        unit_cheby_grid_2d[i,j] = unit_cheby_grid1[i]*unit_cheby_grid2[j]

    lagrange_tensor

    # initialization
    print "initialization"
    
    k_level = tree_length
    x_level = log2(N)-tree_length
    
    next_coefficients = {}
    for i1,i2 in range2(2**x_level):
        next_coefficients[i1,i2] = {}
    
    for j1,j2 in range2(2**k_level):
        
        for i1,i2 in range2(2**x_level):
            next_coefficients[i1,i2][j1,j2] = numpy.zeros((n_cheby,n_cheby))
        
        interp_grid1 = (unit_cheby_grid_1d+j1)*2**-k_level
        interp_grid2 = (unit_cheby_grid_1d+j2)*2**-k_level
        
        if len(sources[j1,j2])>0:
            for i1,i2 in range2(2**x_level):
                x_center = ((i1+0.5)*x_spacing, (i2+0.5)*2**-x_level)
                
                for t1,t2 in range2(n_cheby):
                    interp_point = (interp_grid1[t1], interp_grid2[t2])
                    value = 0.0j
                    for (source_point, source) in sources:
                        prod =  interpolation.lagrange_prod(source_point[0], interp_grid1, t1)
                        prod *= interpolation.lagrange_prod(source_point[1], interp_grid2, t2)
                        prod *= numpy.exp( phase(x_center, interp_point)-phase_p(x_center, source_point) )
                        value += prod*source

                next_coefficients[i1,i2][j1,j2][t1,t2] = value
    
    coefficients = next_coefficients

    k_level -= 1
    x_level += 1
    
    #recurrence, first half
    print "recurrence first half"
    
    while k_level >= log2(N)/2:
        print k_level

        next_coefficients = {}
        for i1,i2 in range2(2**x_level):
            next_coefficients[i1,i2] = {}

        for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
            print '\t', j1,j2
            
            cheby_grid1 = (unit_cheby_grid+j1)*k_spacing
            cheby_grid2 = (unit_cheby_grid+j2)*k_spacing
            
            for i1,i2 in range2(2**x_level): # looping over index of space-boxes
                x_center = ( (0.5+i1)*x_spacing, (0.5+i2)*x_spacing)
                
                points = []
                values = []
                
                for c1, c2 in [ (j1*2,j2*2), (j1*2+1,j2*2), (j1*2,j2*2+1), (j1*2+1,j2*2+1) ]: # looping over index of child frequency-boxes
                    temp = expansion_coefficients_firsthalf(N, phase, x_center, tree[k_level+1][i1/2,i2/2][c1,c2], cheby_grid1, cheby_grid2)
                    
                    for point, value in temp:
                        if point in points:
                            values[points.index(point)] += value
                        else:
                            points.append(point)
                            values.append(value)
                
                tree[k_level][i1,i2][j1,j2] = [ (points[k], values[k]) for k in range(len(points)) ]
        k_level -= 1
        x_level += 1

    print "at middle; switching representation"
    
    switch_representation(N, phase, tree, k_level, n_cheby )

    print "recursion second half"

    while k_level >= 0:
        print k_level
        x_level = log2(N)-k_level
        x_spacing = 2**-x_level
        k_spacing = 2**-k_level

        for i1,i2 in range2(2**x_level):
            tree[k_level][i1,i2] = {}

        for i1,i2 in range2(2**x_level): # looping over index of space-boxes
            print '\t', i1,i2
            
            cheby_grid1 = (unit_cheby_grid+i1)*x_spacing
            cheby_grid2 = (unit_cheby_grid+i2)*x_spacing
            
            for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
                
                points = []
                values = []
                
                for c1, c2 in [ (j1*2,j2*2), (j1*2+1,j2*2), (j1*2,j2*2+1), (j1*2+1,j2*2+1) ]: # looping over index of child frequency-boxes
                    k_center = ( (0.5+j1)*k_spacing, (0.5+j2)*k_spacing)
                    
                    temp = expansion_coefficients_secondhalf(N, phase, k_center, tree[k_level+1][i1/2,i2/2][c1,c2], cheby_grid1, cheby_grid2)
                    
                    for point, value in temp:
                        if point in points:
                            values[points.index(point)] += value
                        else:
                            points.append(point)
                            values.append(value)
                
                tree[k_level][i1,i2][j1,j2] = [ (points[k], values[k]) for k in range(len(points)) ]
        k_level -= 1

    print "termination"

    

    for i1,i2 in range2(N):
        


if __name__=="__main__":
    
    import random
    
    N = 2**3
    n_cheby = 3
    
    input = numpy.zeros((N,N))
    for i in itertools.product(range(N),range(N)):
        input[i] = random.uniform(-1,1)
    
    def phase(x, k):
        return numpy.dot(x,k)/N
    
    bfio_eval(N, phase, input, n_cheby, log2(N)-0)





"""

import random
    
N = 2**4
n_cheby = 3
    
input = numpy.zeros((N,N))
for i in itertools.product(range(N),range(N)):
    input[i] = random.uniform(-1,1)
    
def phase(x, k):
    return numpy.dot(x,k)/N
    
bfio_eval(N, phase, input, n_cheby, log2(N)-0)


"""


