
import numpy
numpy_norm = numpy.linalg.norm
import unittest
import ctypes, ctypes.util

c_interpolation = ctypes.cdll.LoadLibrary(ctypes.util.find_library('bfio_interpolation'))

c_interpolation.lagrange_basis.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_void_p, ctypes.c_int]
c_interpolation.lagrange_basis.restype = ctypes.c_double

c_interpolation.lagrange_interpolate.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p]
c_interpolation.lagrange_interpolate.restype = ctypes.c_double


CONST_PI = 3.1415926535897931
CONST_2PI = 6.2831853071795862


class lagrange_interpolator:
    def __init__(self, support_points):
        self.support = numpy.array(support_points, dtype=numpy.float64)
        self.len = ctypes.c_int(len(self.support))
        self.support_pointer = self.support.ctypes.data_as(ctypes.c_void_p)
    
    def set(self, function):
        self.function_values = numpy.array([function(x) for x in self.support], dtype=numpy.float64)
        self.function_values_pointer = self.function_values.ctypes.data_as(ctypes.c_void_p)
    
    def basis(self, point, index):
        return c_interpolation.lagrange_basis(ctypes.c_double(point), self.len, self.support_pointer, ctypes.c_int(index))
    
    def __call__(self, point):
        return c_interpolation.lagrange_interpolate(ctypes.c_double(point), self.len, self.support_pointer, self.function_values_pointer)
    
    def __getitem__(self, index):
        return self.support[index]



def lagrange_basis(point, support, index):
    
    prod=1.0
    for k in range(len(support)):
        if k != index:
            prod *= (point - support[k])/(support[index]-support[k])
    return prod



def lagrange_matrix(support_points, interpolation_points):
    
    n=len(support_points)
    m=len(interpolation_points)
    L = numpy.zeros((m,n), dtype=numpy.float64)
    
    for i in range(m):
        for j in range(n):
            L[i,j] = lagrange_basis(interpolation_points[i], support_points, j)
    
    return L

class lagrange_tensor:
    def __init__(self, support_points1, interpolation_points1, support_points2, interpolation_points2 ):
        self.left_matrix = lagrange_matrix(support_points1, interpolation_points1)
        self.right_matrix = numpy.transpose( lagrange_matrix(support_points2, interpolation_points2))

    def __call__(self, rhs):
        return numpy.dot( self.left_matrix, numpy.dot(rhs, self.right_matrix ))

def chebyshev_grid(n_points, a=-1, b=1):
    """ returns chebyshev-spaced points on [a,b]
    """
    return numpy.array([ (numpy.cos(float(k)*CONST_PI/(n_points-1))*0.5+0.5)*(b-a)+a for k in range(n_points) ])[::-1]



class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_lagrange_cheby_product(self):
        
        import pylab
        import itertools
        import scipy
        
        test_function = lambda x,y: numpy.cos(CONST_2PI*x*2)*numpy.sin(CONST_2PI*x)
        #numpy.cos(CONST_2PI*x*4)*(float(y <= 0.5)*2-1)*abs(y-0.2)
        
        sample_points=range(10,30)
        errors=[]
        
        n_test_points = 30
        for n_sample_points in sample_points:
            Q = chebyshev_grid(n_sample_points)*0.5+0.5
            U = itertools.product(Q,Q)
            X = [ float(i)/n_test_points for i in range(n_test_points) ]
            W = itertools.product(X,X)
            L = lagrange_matrix(Q, X)
            F = numpy.array( [ [test_function(x,y) for x in Q ] for y in Q ] )
            H = numpy.dot( numpy.dot(L, F), numpy.transpose(L) )
            G = numpy.array( [ [test_function(x,y) for x in X ] for y in X ] )
            errors.append( numpy.log10(numpy_norm(G-H)/numpy_norm(G)+1E-17) )
        
        (slope,offset)=scipy.polyfit(sample_points,errors,1)
        mean = scipy.mean(errors)
        #pylab.plot(sample_points,errors)
        #pylab.show()
        #print slope,mean
        self.assertTrue( slope <= -0.4 and mean <= -4 )
    
    def test_lagrange_interpolator(self ):
        import random
        norm = numpy.linalg.norm
        n_cheby = 30
        g = chebyshev_grid(n_cheby, 0, 1)
        T = lagrange_interpolator(g)
        f = lambda x: x*x*numpy.cos(x*6.3*3)
        T.set(f)
        X = numpy.arange(0,1,0.001)
        Y1 = numpy.array([T(x) for x in X])
        Y2 = numpy.array([f(x) for x in X])
        self.assertTrue(norm(Y1-Y2)/norm(Y2) < 1E-10)
        for k in range(n_cheby):
            x = random.uniform(0,1)
            self.assertTrue(abs(lagrange_basis(x, g, k)-T.basis(x,k)) < 1E-15)






if __name__=="__main__":    

    unittest.main()
    








