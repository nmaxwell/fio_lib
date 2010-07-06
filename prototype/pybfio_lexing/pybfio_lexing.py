
import numpy
import itertools

import interpolation

CONST_SQRT2_2 = 0.70710678118654757
CONST_PI = math.pi
CONST_2PI = 6.2831853071795865

CONST_2PI_J = 6.2831853071795865j
CONST_PI_SQRT2_J = 4.4428829381583661j

def log2(x):
    return (int)(math.floor(math.log(x,2)))

def range2(first, second=None, third=None):
    try:
        if second is None:
            return itertools.product(range(first[0]), range(first[1]))
        if third is None:
            return itertools.product(range(first[0], second[0]), range(first[1], second[1]))
        else:
            return itertools.product(range(first[0], second[0], third[0]), range(first[1], second[1], third[1]))
    except:
        if second is None:
            return itertools.product(range(first), range(first))
        if third is None:
            return itertools.product(range(first, second), range(first, second))
        else:
            return itertools.product(range(first, second, third), range(first, second, third))

def chebyshev_grid(n_points, a=-1, b=1):
    """ returns chebyshev-spaced points on [a,b]
    """
    return numpy.array([ (numpy.cos(float(k)*CONST_PI/(n_points-1))*0.5+0.5)*(b-a)+a for k in range(n_points) ])[::-1]

def ndgrid(grid ):
    N = len(grid)
    grid1 = numpy.zeros((N,N), dtype=numpy.float64)
    grid2 = numpy.zeros((N,N), dtype=numpy.float64)
    grid2[:] = grid
    grid1 = numpy.transpose(grid2)
    return grid1, grid2




def expansion_coefficients(level):
    coefficients = {}
    for i1,i2 in range2(int(N*2**-level)):
        coefficients[i1,i2] = {}
        for j1,j2 in range2((2**level)):
            coefficients[i1,i2][j1,j2] = numpy.zeros((n_cheby,n_cheby), dtype=numpy.complex128)
    return coefficients


def eval_kernel(phase, x_list, k_list):
    kernel = numpy.zeros((len(x_list), len(k_list)), dtype=numpy.complex128)
    for index in range2(len(x_list), len(k_list))
        kernel[index] = numpy.exp(CONST_2PI_J*phase(x_list[index[0]], k_list[index[1]]))
    return kernel






def bfio_eval(N, input_data, phase, n_cheby, start_level, stop_level):
    
    output = numpy.zeros((N,N), dtype=numpy.complex128)
    
    unit_cheby_grid = interpolation.chebyshev_grid(n_cheby, 0, 1)
    mats = {}
    mats[0] = lagrange_matrix(unit_cheby_grid, unit_cheby_grid*0.5)
    mats[1] = lagrange_matrix(unit_cheby_grid, unit_cheby_grid*0.5+0.5)
    
    n_zones = 2**end_level
    zone_size = N/n_zones
    middle_level = (start_level+stop_level)/2
    halfN = 0.5*N
    
    for zone1, zone2 in range2(n_zones):
        print "in zone", zone
        
        k_start = numpy.array((zone1*zone_size, zone2*zone_size))
        k_end = numpy.array(((zone1+)*zone_size, (zone2+1)*zone_size))
        
        level = start_level
        while level >= middle_level:
            n_kboxes = 2**level
            n_xboxes = N/n_kboxes
            kbox_size = float(N)/n_kboxes
            xbox_size = 1.0/n_xboxes
            
            next_coefficients = expansion_coefficients(start_level)
            
            unit_kbox = ndgrid(numpy.arange(0,1,1.0/k_box_size))
            unit_cheby_grid_2d = ndgrid(unit_cheby_grid)
            
            for k1,k2 in range2(k_start/int(kbox_size), k_end/int(kbox_size))
                k = numpy.array((k1,k2))
                k_center = (k+0.5)*kbox_size-halfN
                
                for x1,x2 in range2(n_xboxes):
                    x = numpy.array((x1,x2))
                    x_center = (x+0.5)*xbox_size
                    
                    if level == start_level:
                        section = input_data[k1*kbox_size:(k1+1)*kbox_size, k2*kbox_size:(k2+1)*kbox_size]
                        source_points = [ (k+float(j)/kbox_size)*kbox_size - halfN for j in range(kbox_size)]
                        inside_kernel = eval_kernel(phase, [x_center, ], source_points)
                        outside_kernel = eval_kernel(phase, [x_center, ], source_points)
                        
                        interpolator = lagrange_matrix(numpy.arange(0,1,1.0/k_box_size), unit_cheby_grid)
                        next_coefficients[k1,k2][x1,x2] = numpy.dot(numpy.dot(interpolator, (section*inside_kernel)), numpy.transpose(interpolator))/outside_kernel
            
            level -= 1

