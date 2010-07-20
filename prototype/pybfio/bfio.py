#!/usr/bin/env python
"""Butterfly algorithm for evaluating Fourier Integral Operators.

"""

import math
import itertools
import numpy
import sys
import unittest
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
CONST_PI = math.pi
CONST_2PI = 6.2831853071795865

CONST_2PI_J = 6.2831853071795865j
CONST_PI_SQRT2_J = 4.4428829381583661j

def log2(x):
    return (int)(math.floor(math.log(x,2)))

def range2(n,m=None):
    if m is None:
        m=n
    return itertools.product(range(n),range(m))

def polar(N, k):
    k = numpy.array(k)
    p1 = math.atan2(k[1], k[0])/CONST_2PI+0.5
    p0 = numpy.linalg.norm(k)/(CONST_SQRT2_2*N)
    return numpy.array((p0, p1))




"""

    boxed_output_targets = {}
    for i1,i2 in range2(N):
        a1,a2,b1,b2 = float(j1)/N, float(j2)/N, float(j1+1)/N, float(j2+1)/N
        box = []
        for target_point in output_targets:
            if a1 <= target_point[0] and target_point[0] < b1 and a2 <= target_point[1] and target_point[1] < b2:
                box.append(target_point)
        boxed_output_targets[i1,i2] = box
"""


<<<<<<< HEAD


def bfio_box_sources(level, sources):
    
    boxed_sources = {}
    for j1,j2 in range2(2**level):
        a1,a2,b1,b2 = float(j1)*2**-level, float(j2)*2**-level, float(j1+1)*2**-level, float(j2+1)*2**-level
        box=[]
        next_sources=[]
        while len(sources)>0:
            popped = sources.pop()
            source_point, source = popped
            if source_point[0] <= b1 and source_point[1] <= b2 and ( a1 < source_point[0] or a1 <= source_point[0] and j1==0) and ( a2 < source_point[1] or a2 <= source_point[1] and j2==0):
                box.append(popped )
            else:
                next_sources.append(popped)
        sources = next_sources
        boxed_sources[j1,j2] = box
    
    return boxed_sources



def bfio_polar_prep(N, input_sources, level, input_phase):
    
    phase = lambda x,p: CONST_PI_SQRT2_J*N*p[0]*input_phase(x,(numpy.cos(CONST_2PI*(p[1]-0.5)),numpy.sin(CONST_2PI*(p[1]-0.5))))
    
    halfN = 0.5*N
    sources = [ (polar(N, (float(j1)-halfN, float(j2)-halfN)), input_sources[j1,j2]) for j1,j2 in range2(N) ]
    # recovered_point = numpy.array( (numpy.cos(CONST_2PI*(polar_point[1]-0.5)), numpy.sin(CONST_2PI*(polar_point[1]-0.5))) ) *N*polar_point[0]*CONST_SQRT2_2
    
    boxed_sources = bfio_box_sources(level, sources)
    
    return boxed_sources, phase

=======



def bfio_polar_boxed_sources(N, input_sources, tree_length):

    halfN = 0.5*N
    sources = [ (polar(N, (float(j1)-halfN, float(j2)-halfN)), input_sources[j1,j2]) for j1,j2 in range2(N) ]
    # recovered_point = numpy.array( (numpy.cos(CONST_2PI*(polar_point[1]-0.5)), numpy.sin(CONST_2PI*(polar_point[1]-0.5))) ) *N*polar_point[0]*CONST_SQRT2_2
    
    boxed_sources = {}
    for j1,j2 in range2(2**tree_length):
        a1,a2,b1,b2 = float(j1)*2**-tree_length, float(j2)*2**-tree_length, float(j1+1)*2**-tree_length, float(j2+1)*2**-tree_length
        box=[]
        next_sources=[]
        while len(sources)>0:
            popped = sources.pop()
            source_point, source = popped
            if source_point[0] <= b1 and source_point[1] <= b2 and ( a1 < source_point[0] or a1 <= source_point[0] and j1==0) and ( a2 < source_point[1] or a2 <= source_point[1] and j2==0):
                box.append(popped )
            else:
                next_sources.append(popped)
        sources = next_sources
        boxed_sources[j1,j2] = box

    return boxed_sources



def bfio_prototype(N, input_sources, input_phase, n_cheby, tree_length):
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

    # preliminaries

<<<<<<< HEAD



def bfio_prototype(N, sources, phase, n_cheby, start_level, end_level):

    # preliminaries

    unit_cheby_grid = interpolation.chebyshev_grid(n_cheby, 0.0, 1.0)

    child_tensor = {}
    for c1,c2 in range2(2,2):
        child_tensor[c1,c2] = interpolation.lagrange_tensor(unit_cheby_grid, (unit_cheby_grid+c1)*0.5, unit_cheby_grid, (unit_cheby_grid+c2)*0.5)
    
    def new_interp_grid_zeros():
        return numpy.zeros((n_cheby, n_cheby), dtype=numpy.complex128)
    
    def expansion_coefficients(level):
        coefficients = {}
        for i1,i2 in range2(int(N*2**-level)):
            coefficients[i1,i2] = {}
            for j1,j2 in range2(2**level):
                coefficients[i1,i2][j1,j2] = new_interp_grid_zeros()
        return coefficients
    
    middle_level = (start_level+end_level)/2
    
    # initialization
    print "initialization"
    
    k_level = start_level
    x_level = log2(N)-start_level
    next_coefficients = expansion_coefficients(k_level)
    
    for j1,j2 in range2(2**k_level):
        if len(sources[j1,j2]) > 0:
            for i1,i2 in range2(2**x_level):
                print j1, i1
                lagrange1=interpolation.lagrange_interpolator((unit_cheby_grid+i1)*2**-k_level)
                lagrange2=interpolation.lagrange_interpolator((unit_cheby_grid+i2)*2**-k_level)
                x_center = ((i1+0.5)*2**-x_level, (i2+0.5)*2**-x_level)
                for t1,t2 in range2(n_cheby):
                    target_point = (lagrange1[t1], lagrange2[t2])
                    value = 0
                    for (source_point, source) in sources[j1,j2]:
                        prod =  lagrange1.basis(source_point[0], t1)
                        prod *= lagrange2.basis(source_point[1], t2)
                        prod *= numpy.exp( +phase(x_center, source_point))
                        value += prod*source
                    next_coefficients[i1,i2][j1,j2][t1,t2] = value*numpy.exp( -phase(x_center, source_point))
    
    return 0
    coefficients = next_coefficients
    k_level -= 1
    x_level += 1
    
    #recurrence, first half
    print "recurrence first half"

    while k_level >= middle_level:
        print k_level

        next_coefficients = expansion_coefficients(k_level)
        kernel_child = new_interp_grid_zeros()
        source_points = new_interp_grid_zeros()

        for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
            for i1,i2 in range2(2**x_level): # looping over index of space-boxes
                x_center = ((0.5+i1)*2**-x_level, (0.5+i2)*2**x_level)

                for c1,c2 in [ (0,0), (0,1), (1,0), (1,1) ]: # looping over children
                    jc1, jc2 = j1*2+c1, j2*2+c2 # index of child frequency-boxes, at level k_level+1

                    for t1,t2 in range2(n_cheby):
                        source_point = ((jc1+unit_cheby_grid[t1])*2**-(k_level+1), (jc2+unit_cheby_grid[t2])*2**-(k_level+1))
                        kernel_child[t1,t2] = numpy.exp( +phase(x_center, source_point))

                    next_coefficients[i1,i2][j1,j2] += child_tensor[c1,c2]( kernel_child*coefficients[i1/2,i2/2][jc1,jc2] ) # element-wise pruduct

                for t1,t2 in range2(n_cheby):
                    target_point = ( (j1+unit_cheby_grid[t1])*2**-(k_level), (j2+unit_cheby_grid[t2])*2**-(k_level) )
                    next_coefficients[i1,i2][j1,j2] *= numpy.exp( -phase(x_center, target_point))
        
        coefficients = next_coefficients
        k_level -= 1
        x_level += 1
    
    #at middle; switching representation
    print "at middle; switching representation"

    k_level += 1
    x_level -= 1
    next_coefficients = expansion_coefficients(k_level)

    for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
        for i1,i2 in range2(2**x_level): # looping over index of space-boxes

            for t1,t2 in range2(n_cheby):
                target_point = ((i1+unit_cheby_grid[t1])*2**-(x_level), (i2+unit_cheby_grid[t2])*2**-(x_level))
                value = 0.0j
                for s1,s2 in range2(n_cheby):
                    source_point = ((j1+unit_cheby_grid[s1])*2**-(k_level), (j2+unit_cheby_grid[s2])*2**-(k_level))
                    value += coefficients[i1,i2][j1,j2][s1,s2]*numpy.exp( +phase(target_point, source_point))

                next_coefficients[i1,i2][j1,j2][t1,t2] = value

    coefficients = next_coefficients
    k_level -= 1
    x_level += 1
=======
    phase = lambda x,p: CONST_PI_SQRT2_J*N*p[0]*input_phase(x,(numpy.cos(CONST_2PI*(p[1]-0.5)),numpy.sin(CONST_2PI*(p[1]-0.5))))

    sources = bfio_polar_boxed_sources(N, input_sources, tree_length)

    unit_cheby_grid = interpolation.chebyshev_grid(n_cheby, 0, 1)

    child_tensor = {}
    for c1,c2 in range2(2,2):
        child_tensor[c1,c2] = interpolation.lagrange_tensor((unit_cheby_grid+c1)*0.5, unit_cheby_grid, (unit_cheby_grid+c2)*0.5, unit_cheby_grid)

    def expansion_coefficients(level):
        coefficients = {}
        for i1,i2 in range2(int(N*2**-level)):
            coefficients[i1,i2] = {}
            for j1,j2 in range2(2**level):
                coefficients[i1,i2][j1,j2] = numpy.zeros((n_cheby,n_cheby))
        return coefficients

    # initialization
    print "initialization"

    k_level = tree_length
    x_level = log2(N)-tree_length
    next_coefficients = expansion_coefficients(k_level)
    
    for j1,j2 in range2(2**k_level):
        for i1,i2 in range2(2**x_level):
            interp_grid1 = (unit_cheby_grid+j1)*2**-k_level
            interp_grid2 = (unit_cheby_grid+j2)*2**-k_level
            if len(sources[j1,j2]) > 0:
                for i1,i2 in range2(2**x_level):
                    x_center = ((i1+0.5)*2**-x_level, (i2+0.5)*2**-x_level)

                    for t1,t2 in range2(n_cheby):
                        interp_point = (interp_grid1[t1], interp_grid2[t2])
                        value = 0
                        for (source_point, source) in sources[j1,j2]:
                            prod =  interpolation.lagrange_prod(source_point[0], interp_grid1, t1)
                            prod *= interpolation.lagrange_prod(source_point[1], interp_grid2, t2)
                            prod *= numpy.exp( +phase(x_center, source_point))
                            value += prod*source
                        next_coefficients[i1,i2][j1,j2][t1,t2] = value*numpy.exp( -phase(x_center, interp_point))
    coefficients = next_coefficients
    k_level -= 1
    x_level += 1

    #recurrence, first half
    print "recurrence first half"

    while k_level >= log2(N)/2:
        print k_level

        next_coefficients = expansion_coefficients(k_level)
        kernel_child = numpy.zeros((n_cheby, n_cheby))
        source_points = numpy.zeros((n_cheby, n_cheby))

        for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
            for i1,i2 in range2(2**x_level): # looping over index of space-boxes
                x_center = ((0.5+i1)*2**-x_level, (0.5+i2)*2**x_level)

                for c1,c2 in [ (0,0), (0,1), (1,0), (1,1) ]: # looping over children
                    jc1, jc2 = j1*2+c1, j2*2+c2 # index of child frequency-boxes, at level k_level+1

                    for t1,t2 in range2(n_cheby):
                        source_point = ((jc1+unit_cheby_grid[t1])*2**-(k_level+1), (jc2+unit_cheby_grid[t2])*2**-(k_level+1))
                        kernel_child[t1,t2] = numpy.exp( +phase(x_center, source_point))

                    next_coefficients[i1,i2][j1,j2] += child_tensor[c1,c2]( kernel_child*coefficients[i1/2,i2/2][jc1,jc2] ) # element-wise pruduct

                for t1,t2 in range2(n_cheby):
                    target_point = ( (j1+unit_cheby_grid[t1])*2**-(k_level), (j2+unit_cheby_grid[t2])*2**-(k_level) )
                    next_coefficients[i1,i2][j1,j2] *= numpy.exp( -phase(x_center, target_point))
        
        coefficients = next_coefficients
        k_level -= 1
        x_level += 1
    
    #at middle; switching representation
    print "at middle; switching representation"

    k_level += 1
    x_level -= 1
    next_coefficients = expansion_coefficients(k_level)

    for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
        for i1,i2 in range2(2**x_level): # looping over index of space-boxes

            for t1,t2 in range2(n_cheby):
                target_point = ((i1+unit_cheby_grid[t1])*2**-(x_level), (i2+unit_cheby_grid[t2])*2**-(x_level))
                value = 0.0j
                for s1,s2 in range2(n_cheby):
                    source_point = ((j1+unit_cheby_grid[s1])*2**-(k_level), (j2+unit_cheby_grid[s2])*2**-(k_level))
                    value += coefficients[i1,i2][j1,j2][s1,s2]*numpy.exp( +phase(target_point, source_point))

                next_coefficients[i1,i2][j1,j2][t1,t2] = value

    coefficients = next_coefficients
    k_level -= 1
    x_level += 1

    parent_tensor = {}
    for c1,c2 in range2(2,2):
        parent_tensor[c1,c2] = interpolation.lagrange_tensor( unit_cheby_grid, (unit_cheby_grid+c1)*0.5, unit_cheby_grid, (unit_cheby_grid+c2)*0.5)
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

    # recursion, second half
    print "recursion, second half"

<<<<<<< HEAD
    while k_level >= end_level:
        print k_level
        
        next_coefficients = expansion_coefficients(k_level)
        kernel_parent = new_interp_grid_zeros()
        source_points = new_interp_grid_zeros()

        for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
            for i1,i2 in range2(2**x_level): # looping over index of space-boxes
                
                for c1,c2 in [ (0,0), (0,1), (1,0), (1,1) ]: # looping over children
                    jc1,jc2 = j1*2+c1, j2*2+c2 # index of child frequency-boxes, at level k_level+1
                    child_center = ((0.5+jc1)*2**-(k_level+1), (0.5+jc2)*2**-(k_level+1))
=======
    while k_level >= 0:
        print k_level

        next_coefficients = expansion_coefficients(k_level)
        kernel_parent = numpy.zeros((n_cheby, n_cheby))
        source_points = numpy.zeros((n_cheby, n_cheby))

        for j1,j2 in range2(2**k_level): # looping over index of frequency-boxes
            for i1,i2 in range2(2**x_level): # looping over index of space-boxes

                for c1,c2 in [ (0,0), (0,1), (1,0), (1,1) ]: # looping over children
                    jc1,jc2 = j1*2+c1, j2*2+c2 # index of child frequency-boxes, at level k_level+1
                    child_center = ((0.5+jc1)**2**-(k_level+1), (0.5+jc2)*2**-(k_level+1))
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

                    for t1,t2 in range2(n_cheby):
                        source_point = ((i1/2+unit_cheby_grid[t1])*2**-(x_level-1), (i2/2+unit_cheby_grid[t2])*2**-(x_level-1))
                        kernel_parent[t1,t2] = numpy.exp( -phase(source_point, child_center))
                        
<<<<<<< HEAD
                    temp = child_tensor[c1,c2]( kernel_parent*coefficients[i1/2,i2/2][jc1,jc2] )
=======
                    temp = parent_tensor[c1,c2]( kernel_parent*coefficients[i1/2,i2/2][jc1,jc2] )
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

                    for t1,t2 in range2(n_cheby):
                        target_point = ((i1+unit_cheby_grid[t1])*2**-x_level, (i2+unit_cheby_grid[t2])*2**-x_level)
                        temp[t1,t2] *= numpy.exp( +phase(target_point, child_center))
#                    print j1,j2,i1,i2,c1,c2, numpy.linalg.norm( temp )
                    next_coefficients[i1,i2][j1,j2] += temp

        coefficients = next_coefficients
        k_level -= 1
        x_level += 1
    
    k_level += 1
    x_level -= 1

    result = {}
    for i1,i2 in range2(N):
        result[i1,i2] = coefficients[i1,i2][0,0]

    return result, N, n_cheby, unit_cheby_grid, phase

<<<<<<< HEAD
    # termination
    print "termination"

=======

    # termination
    print "termination"

>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

def bfio_interpolate_results(target_points, bfio_output):
    coefficients, N, n_cheby, unit_cheby_grid, phase = bfio_output[0:5]
    k_center = (0.5, 0.5)
    
    result=[]
    for target in target_points:
        i1 = int(numpy.floor(target[0]*N))
        i2 = int(numpy.floor(target[1]*N))
        if 0 <= i1 and i1 < N and 0 <= i2 and i2 < N:
            interp_grid1 = (unit_cheby_grid+i1)*2**-log2(N)
            interp_grid2 = (unit_cheby_grid+i2)*2**-log2(N)
            
            value = 0
            for t1,t2 in range2(n_cheby):
                source_point = (interp_grid1[t1], interp_grid2[t2])
<<<<<<< HEAD
                prod =  interpolation.lagrange_basis(target[0], interp_grid1, t1)
                prod *= interpolation.lagrange_basis(target[1], interp_grid2, t2)
=======
                prod =  interpolation.lagrange_prod(target[0], interp_grid1, t1)
                prod *= interpolation.lagrange_prod(target[1], interp_grid2, t2)
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
                prod *= numpy.exp( -phase(source_point, k_center))
                prod *= coefficients[i1,i2][t1,t2]
                value += prod
            result.append( value*numpy.exp( +phase(target, k_center)) )
        else:
            result.append(0)
    return result
<<<<<<< HEAD





def test_case():
    n = 6
    N = 2**n
    
    import scipy
    import scipy.io
    
    fname = "/home/nick/Desktop/radar/f_%d.mat"%N
    print fname
    loaded = scipy.io.loadmat(fname )
    loaded_sources = loaded['f']
    sources = []
    for i,j in range2(N):
        source_point = numpy.array((i,j), dtype=numpy.float64)/N
        source = loaded_sources[i,j]
        sources.append( (source_point, source) )
    
    n_cheby = 5
    start_level = n-2
    end_level = 2
    phase = lambda x,k: numpy.dot(x,k)
    
    boxed_sources = bfio_box_sources(start_level, sources)
    bfio_prototype(N, boxed_sources, phase, n_cheby, start_level, end_level)






class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_polar(self):
        norm = numpy.linalg.norm
        """
        from math import *
        X = numpy.arange(-1,1,0.1)
        Y = itertools.product(X,X)
        for x in Y:
            x = numpy.array((x[0], x[1]))
            x /= norm(x)
            print atan2(x[1], x[0]), norm(x - numpy.array((cos(atan2(x[1], x[0])), sin(atan2(x[1], x[0])))))
        """

        for n in range(4,7):
            N = 2**n
            for point in [ numpy.array((float(k1-N/2),float(k2-N/2))) for k1,k2 in range2(N) ]:
                polar_point = polar(N, point)
                self.assertTrue( 0.0 <= polar_point[0] and polar_point[0] <= 1.0 )
                self.assertTrue( 0.0 <= polar_point[1] and polar_point[1] <= 1.0 )
                recovered_point = numpy.array( (numpy.cos(CONST_2PI*(polar_point[1]-0.5)), numpy.sin(CONST_2PI*(polar_point[1]-0.5))) ) *N*polar_point[0]*CONST_SQRT2_2
                if norm(point-recovered_point) >= 2E-13:
                    print norm(point-recovered_point)
                self.assertTrue( norm(point - recovered_point) < 2E-13 )

    def test_interpolation(self):
        norm = numpy.linalg.norm
        import random
        

        grid = interpolation.chebyshev_grid(13, 0.0, 1.0)
        for x in grid:
            self.assertTrue( 0.0 <= x and x <= 1.0 )
        self.assertTrue( abs(grid[0]-0.0) < 1E-15 and abs(grid[-1]-1) < 1E-15 )

        n_in = 12
        rand = lambda a: random.uniform(-a,a)
        input_X = interpolation.chebyshev_grid(n_in, rand(0.3)-1, rand(0.3)+1)
        input_Y = interpolation.chebyshev_grid(n_in, rand(0.3)-1, rand(0.3)+1)
        output_X = numpy.arange(rand(0.25)-0.25, rand(0.25)+0.25, 0.01)
        output_Y = numpy.arange(rand(0.25)-0.25, rand(0.25)+0.25, 0.01)
        T = interpolation.lagrange_tensor(input_X, output_X, input_Y, output_Y)

        #print T.left_matrix
        #print T.right_matrix
        test_function = lambda x,y: numpy.cos(x*y*6)
        
        input = numpy.zeros((n_in, n_in))
        output = numpy.zeros((len(output_X), len(output_Y)))
        for i,j in range2(n_in):
            input[i,j] = test_function(input_X[i], input_Y[j])
        for i,j in range2(len(output_X), len(output_Y)):
            output[i,j] = test_function(output_X[i], output_Y[j])
        self.assertTrue( norm( T(input)-output )/norm(output) < 1E-4 )

    def test_child_tensor(self ):
        norm = numpy.linalg.norm
        n_cheby = 20
        unit_cheby_grid = interpolation.chebyshev_grid(n_cheby, 0.0, 1.0)

        child_tensor = {}
        for c1,c2 in range2(2,2):
            child_tensor[c1,c2] = interpolation.lagrange_tensor(unit_cheby_grid, (unit_cheby_grid+c1)*0.5, unit_cheby_grid, (unit_cheby_grid+c2)*0.5)
        
        test_function = lambda x,y: numpy.cos(6.3*x)*numpy.sin(6.3*y)
        
        for c1,c2 in range2(2):
            source = numpy.zeros((n_cheby, n_cheby))
            for t1,t2 in range2(n_cheby):
                x = unit_cheby_grid[t1]
                y = unit_cheby_grid[t2]
                source[t1, t2] = test_function(x,y)
            output = child_tensor[c1,c2](source)
            exact = numpy.zeros((n_cheby, n_cheby))
            for t1,t2 in range2(n_cheby):
                x=(unit_cheby_grid[t1]+c1)*0.5
                y=(unit_cheby_grid[t2]+c2)*0.5
                exact[t1, t2] = test_function(x,y)
            self.assertTrue( norm(exact-output)/norm(exact) <= 1E-14 )

    def test_bfio_functionality(self ):
        test_case()

        
        
        
        
        
        
        
"""
        
=======






class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_polar(self):
        norm = numpy.linalg.norm
        """
        from math import *
        X = numpy.arange(-1,1,0.1)
        Y = itertools.product(X,X)
        for x in Y:
            x = numpy.array((x[0], x[1]))
            x /= norm(x)
            print atan2(x[1], x[0]), norm(x - numpy.array((cos(atan2(x[1], x[0])), sin(atan2(x[1], x[0])))))
        """

        for n in range(4,7):
            N = 2**n
            for point in [ numpy.array((float(k1-N/2),float(k2-N/2))) for k1,k2 in range2(N) ]:
                polar_point = polar(N, point)
                self.assertTrue( 0.0 <= polar_point[0] and polar_point[0] <= 1.0 )
                self.assertTrue( 0.0 <= polar_point[1] and polar_point[1] <= 1.0 )
                recovered_point = numpy.array( (numpy.cos(CONST_2PI*(polar_point[1]-0.5)), numpy.sin(CONST_2PI*(polar_point[1]-0.5))) ) *N*polar_point[0]*CONST_SQRT2_2
                if norm(point-recovered_point) >= 2E-13:
                    print norm(point-recovered_point)
                self.assertTrue( norm(point - recovered_point) < 2E-13 )

    def test_interpolation(self):
        norm = numpy.linalg.norm
        import random
        

        grid = interpolation.chebyshev_grid(13, 0.0, 1.0)
        for x in grid:
            self.assertTrue( 0.0 <= x and x <= 1.0 )
        self.assertTrue( abs(grid[0]-0.0) < 1E-15 and abs(grid[-1]-1) < 1E-15 )

        n_in = 12
        rand = lambda a: random.uniform(-a,a)
        input_X = interpolation.chebyshev_grid(n_in, rand(0.3)-1, rand(0.3)+1)
        input_Y = interpolation.chebyshev_grid(n_in, rand(0.3)-1, rand(0.3)+1)
        output_X = numpy.arange(rand(0.25)-0.25, rand(0.25)+0.25, 0.01)
        output_Y = numpy.arange(rand(0.25)-0.25, rand(0.25)+0.25, 0.01)
        T = interpolation.lagrange_tensor(input_X, output_X, input_Y, output_Y)
        test_function = lambda x,y: numpy.cos(x*y*6)
        
        input = numpy.zeros((n_in, n_in))
        output = numpy.zeros((len(output_X), len(output_Y)))
        for i,j in range2(n_in):
            input[i,j] = test_function(input_X[i], input_Y[j])
        for i,j in range2(len(output_X), len(output_Y)):
            output[i,j] = test_function(output_X[i], output_Y[j])
        self.assertTrue( norm( T(input)-output )/norm(output) < 1E-4 )



    def test_bfio_polar_boxed_sources(self ):
        n = 4
        N = 2**n
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
        sources = numpy.zeros((N,N))
        result = bfio_polar_boxed_sources(N, sources, log2(N))
        
        self.assertTrue( N*N == sum([len(result[i]) for i in result]))
<<<<<<< HEAD


    def test_test(self ):
        norm = numpy.linalg.norm
        import random
        
        N = 16
        
        halfN = 0.5*N
        X0 = numpy.array((7,7), dtype=numpy.float64)
        
        sources = numpy.zeros((N,N), dtype=numpy.complex128)
        for i,j in range2(N):
            K = numpy.array((i,j), dtype=numpy.float64) - halfN
            sources[i,j] = float( norm(X0-K) < 1E-3 )
            #numpy.exp(CONST_2PI_J*numpy.dot(X,K))/N**2
        
        targets = []
        for i,j in range2(N):
            targets.append((float(i)/N, float(j)/N))
        
        c_phase = lambda x1,x2,k1,k2: -x1*(k1+halfN) - x2*(k2+halfN)
        phase = lambda x,k: c_phase(x[0],x[1],k[0],k[1])

        n_cheby = 4
        tree_length = log2(N)-0
        bfio_output = bfio_prototype(N, sources,  phase, n_cheby, tree_length)
        flat_out = bfio_interpolate_results(targets, bfio_output)
        bfio_out = numpy.zeros((N,N), dtype=numpy.complex128)
        for index,target in enumerate(targets):
            bfio_out[int(target[0]*N), int(target[1]*N)] = flat_out[index]
        
        print "computing exact fio"
        import exact_fio
        fio_out = exact_fio.exact_fio_2d(sources, N, c_phase)
        for i in range2(N):
            if abs(fio_out[i].real)< 1E-10:
                fio_out[i] = fio_out[i].imag*1.0j
            if abs(fio_out[i].imag)< 1E-10:
                fio_out[i] = fio_out[i].real
            
        print "done."


        print "bfio_out: \n"
        print bfio_out[0:N,0:N]
        print "\n\n\nexact fio:\n"
        print fio_out[0:N,0:N]
=======


    def test_bfio_eval(self ):
        norm = numpy.linalg.norm
        import random
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160

        N = 2**4
        n_cheby = 4

        sources = numpy.zeros((N,N), dtype=numpy.complex128)
        for i in itertools.product(range(N),range(N)):
            sources[i] = random.uniform(-1,1)+random.uniform(-1,1)*1.0j

<<<<<<< HEAD
"""

    def _test_bfio_eval(self ):
        norm = numpy.linalg.norm
        import random

        N = 16
        n_cheby = 6

        sources = numpy.zeros((N,N), dtype=numpy.complex128)
        for i in itertools.product(range(N),range(N)):
            sources[i] = random.uniform(-1,1)+random.uniform(-1,1)*1.0j

        targets = []
        for i,j in range2(N):
            targets.append((float(i)/N, float(j)/N))
        
        halfN = 0.5*N
        c_phase = lambda x1,x2,k1,k2: -x1*(k1+halfN) - x2*(k2+halfN)
        phase = lambda x,k: c_phase(x[0],x[1],k[0],k[1])

=======
        targets = []
        for i,j in range2(N):
            targets.append((float(i)/N, float(j)/N))
        
        halfN = 0.5*N
        def phase(x, k):
            return numpy.dot(x,numpy.array(k)+halfN)/N

>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
        tree_length = log2(N)-0

        bfio_output = bfio_prototype(N, sources,  phase, n_cheby, tree_length)
        flat_out = bfio_interpolate_results(targets, bfio_output)
        bfio_out = numpy.zeros((N,N), dtype=numpy.complex128)
        for index,target in enumerate(targets):
            bfio_out[int(target[0]*N), int(target[1]*N)] = flat_out[index]
        
        print "computing exact fio"
        import exact_fio
<<<<<<< HEAD
=======
        c_phase = lambda x1,x2,k1,k2: phase((x1,x2), (k1,k2))
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
        fio_out = exact_fio.exact_fio_2d(sources, N, c_phase)
        print "done."


        print "bfio_out: \n"
<<<<<<< HEAD
        print bfio_out[0:4,0:4]
        print "\n\n\nexact fio:\n"
        print fio_out[0:4,0:4]
=======
        print bfio_out
        print "\n\n\nexact fio:\n"
        print fio_out
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160



<<<<<<< HEAD
if __name__=="__main__":

    args = sys.argv[1:]

    if len(args)>0 and args[0] == "profile":
        import cProfile
        cProfile.run('test_case ()', sort=1)
    else:
        unittest.main()
=======

if __name__=="__main__":

    unittest.main()
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
