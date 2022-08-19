# distutils: language = c++
'''
    // choose branching variable based on given input
    switch(branch_var){

        case 1:
            choosebranchingvar_max_var_(Sigma, &n, item.u, &i);
            break;

        case 2:
            choosebranchingvar_random_(Sigma, &n, item.u, &i);
            break;

        default:
            choosebranchingvar_greedy_(Sigma, &n, item.u, &i);
'''

from cpp_version cimport cython_wrapper
cimport numpy as np
import numpy as np
# from libcpp
cdef int run_BnB_flattened(double[::1] Sigma, int n, int k, double tol, int depth, int rounds, str filename, int branch_var):
    '''Syntax of this function is
    run_BnB(Sigma, n, k, tol, depth, rounds, filename, branch_var)
    Inputs:
    - Sigma is the flattened covariance matrix. Type: n^2 ndarray of type float64

    - n is the number of candidate sensor locations. Type: integer

    - k is the number of sensors to place. Type: integer

    - tol is the numerical tolerance to solve the SDP relaxations to. Type: float

    - depth is the number of assignments to make before resorting to brute force. Type: integer

    - rounds is the number of random roundings to perform from the relaxed SDP solution. Type: integer

    - filename gives the filename to save the results to. Type: string
    
    - branch_var gives the branching strategy as below. Type: integer
        case 0:
            Greedy Good
        case 1:
            Greedy Bad
        case 2:
            Max Variance
        case 3:
            Random
    '''
    out = cython_wrapper(&Sigma[0], n, k, tol, depth, rounds, filename.encode('UTF-8'), branch_var)
    return out
        
cpdef int run_BnB(np.ndarray[np.float64_t, ndim=2] Sigma, int k, double tol, int depth,
                 int rounds, str filename, int branch_var):
    '''Syntax of this function is:
    run_BnB(Sigma, k, tol, depth, rounds, filename, branch_var)

    Inputs:
    - Sigma is the covariance matrix. Type: 2D ndarray of type float64

    - k is the number of sensors to place. Type: integer

    - tol is the numerical tolerance to solve the SDP relaxations to. Type: float

    - depth is the number of assignments to make before resorting to brute force. Type: integer

    - rounds is the number of random roundings to perform from the relaxed SDP solution. Type: integer

    - filename gives the filename to save the results to. Type: string
    
    - branch_var gives the branching strategy as below. Type: integer
        case 0:
            Greedy Good
        case 1:
            Greedy Bad
        case 2:
            Max Variance
        case 3:
            Random
    '''
    out = run_BnB_flattened(np.ascontiguousarray(Sigma).flatten(), Sigma.shape[0], k, tol, depth,
                            rounds, filename, branch_var)
    return out

   
