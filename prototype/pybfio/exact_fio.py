



CONST_2PI_J = 6.2831853071795865j


import ctypes, ctypes.util
import unittest
import itertools
import time as time_module

import numpy

c_exact_fio = ctypes.cdll.LoadLibrary(ctypes.util.find_library('exact_fio'))

<<<<<<< HEAD

C_PHASE_2D = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)

c_exact_fio.exact_fio_dft.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p ]



c_exact_fio.exact_fio_2d.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, C_PHASE_2D ]
=======
C_PHASE_2D = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double)

c_exact_fio.exact_fio_2d.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, C_PHASE_2D ]
c_exact_fio.exact_fio_dft.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p ]
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160


def exact_fio_2d(input_data, N, phase_function):
    output_data = numpy.zeros((N, N), dtype=numpy.complex128)
    if input_data.dtype != numpy.complex128:
        print 'bad'
    
    c_phase = C_PHASE_2D( phase_function)

    try:
        err = c_exact_fio.exact_fio_2d( output_data.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(N), input_data.ctypes.data_as(ctypes.c_void_p), c_phase  )
    except:
        print "exact_fio.py error: exact_fio_2d, exceptioin"

    return output_data


def exact_fio_dft(input_data, N):
    output_data = numpy.zeros((N, N), dtype=numpy.complex128)
    if input_data.dtype != numpy.complex128:
        print 'bad'
    
    try:
        err = c_exact_fio.exact_fio_dft( output_data.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(N), input_data.ctypes.data_as(ctypes.c_void_p)  )
    except:
        print "exact_fio.py error: exact_fio_dft, exception"

    return output_data






def range2(n,m=None):
    if m is None:
        m=n
    return itertools.product(range(n),range(m))



def python_exact_fio_2d(input_data, N, phase):
    halfN = N*0.5
    output_data = numpy.zeros((N,N), dtype=numpy.complex128)
    for i1,i2 in range2(N):
        value = 0
        for j1,j2 in range2(N):
            value += input_data[j1,j2]*numpy.exp(CONST_2PI_J*phase(float(i1)/N, float(i2)/N, j1-halfN, j2-halfN))
        output_data[i1,i2] = value
    return output_data

def python_exact_fio_dft(input_data, N):
    output_data = numpy.zeros((N,N), dtype=numpy.complex128)
    for i1,i2 in range2(N):
        value = 0
        for j1,j2 in range2(N):
            value += input_data[j1,j2]*numpy.exp(-CONST_2PI_J*float(i1*j1+i2*j2 )/N)
        output_data[i1,i2] = value
    return output_data


class TestSequenceFunctions(unittest.TestCase):
    
    def setUp(self ):
        pass

    def test_exact_fio_dft(self ):
<<<<<<< HEAD
        print "test_exact_fio_dft"
=======
        print "text_exact_fio_dft"
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
        import random
        norm = numpy.linalg.norm
        rand = lambda :random.uniform(-1,1)
        N = 32
        input_data = numpy.zeros((N, N), dtype=numpy.complex128)
        for j1,j2 in range2(N):
            input_data[j1,j2] = rand()+rand()*1.0j
        t1 = time_module.clock()
        c_output_data = exact_fio_dft(input_data, N)
        t2 = time_module.clock()
        c_time_elapsed = t2-t1
        t1 = time_module.clock()
        py_output_data = python_exact_fio_dft(input_data, N)
        t2 = time_module.clock()
        py_time_elapsed = t2-1
        print "error: ", norm(c_output_data - py_output_data)/norm(py_output_data)
        print "c time:\t", c_time_elapsed
        print "py time:\t", py_time_elapsed



    def test_exact_fio_2d(self  ):
        print "test_exact_fio_2d"
        import random
        norm = numpy.linalg.norm
        rand = lambda :random.uniform(-1,1)
        N = 32
        input_data = numpy.zeros((N, N), dtype=numpy.complex128)
        for j1,j2 in range2(N):
            input_data[j1,j2] = rand()+rand()*1.0j
        halfN = 0.5*N
        phase = lambda x1,x2,k1,k2: -x1*(k1+halfN) - x2*(k2+halfN)
        t1 = time_module.clock()
        c_output_data = exact_fio_dft(input_data, N)
        t2 = time_module.clock()
        c_time_elapsed = t2-t1
        t1 = time_module.clock()
<<<<<<< HEAD
        py_output_data = exact_fio_2d(input_data, N, phase)
=======
        py_output_data = python_exact_fio_2d(input_data, N, phase)
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
        t2 = time_module.clock()
        py_time_elapsed = t2-1
        print "error: ", norm(c_output_data - py_output_data)/norm(py_output_data)
        print "c time:\t", c_time_elapsed
        print "py time:\t", py_time_elapsed

if __name__ == "__main__":
    
    unittest.main()



<<<<<<< HEAD
=======




























>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
