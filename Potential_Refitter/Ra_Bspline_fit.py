import Bspline_fits_2eRM as bf
import numpy as np
import scipy.optimize as OPT

CMperAU = 2.194746313632e5
IRa = 81_842.5/CMperAU

es12 = np.array([    0.  , 43_405.05, 59_165.3,  66_838.72,  71_172.88])

nist_e_s = np.column_stack((es12/CMperAU-IRa,[n-1 for n in range(7,12)]))

ep12 = np.array([ 21_351.3259, 50_606.29 ])


ep32 = np.array([26_208.8538, 52_392.02 ])

nist_e_pav = np.column_stack((((2*ep12+4*ep32)/6)[:]/CMperAU-IRa,[n-2 for n in range(7,10)]))

ed32 = np.array([ 12_084.2721,  48_744.04, 61_735.0, 68_264.07, 72_043.33, 74_434.30 ])

ed52 = np.array([  13_742.99, 49_240.45,  61_973.8,  68_394.87, 72_123.74, 74_489  ])

nist_e_dav = np.column_stack((((4*ed32+6*ed52)/10)/CMperAU-IRa,[n-3 for n in range(6,13)]))


ef52 = np.array([48_987.86, 59_517, 66_521.85])
ef72 = np.array([49_272.14, 59_814, 66_691.18])

nist_e_fav = np.column_stack((((6*ef52+8*ef72)/14)/CMperAU-IRa,[n-4 for n in range(5,8)]))

import time
ts = time.time()
R0 = 160
Alpha = 18.0
Zion = 88.0
kk = 7
def s_opt():
    with open("Optimized_params_s.txt", mode="w") as file:
        x0 = [3.7702,4.9928,1.5179,1.3691]
        print("Starting with s")
        print("Initial value", bf.opt(x0))
        res_s = OPT.minimize(bf.optD, args=(0,nist_e_s,R0, kk, Alpha, Zion),x0=x0,method='Nelder-Mead')
        s_params = res_s.x
        file.write("s params  ")
        file.write(("%.8f   "*len(res_s.x))%tuple(res_s.x))
        file.write("\n")
        file.flush()
        print("S done, check file")
        print("Paramaters", res_s.x)
        print("Success on the minimization: ",res_s.success)
        print('Minimum value: ',res_s.fun)
        print('Optimizer message: ',res_s.message)
    return 0


def p_opt():
    with open("Optimized_params_p.txt", mode="w") as file:
        x0 = [3.9430,5.0552,3.6770,1.0924]
        print("Starting with p")
        print("Initial value", bf.opt(x0))
        res_s = OPT.minimize(bf.opt, args=(1,nist_e_pav,R0, kk, Alpha, Zion),x0=s_params)
        file.write("p params  ")
        file.write(("%.8f   "*len(res_s.x))%tuple(res_s.x))
        file.write("\n")
        file.flush()
        print("P done, check file")
        print("Paramaters", res_s.x)
        print("Success on the minimization: ",res_s.success)
        print('Minimum value: ',res_s.fun)
        print('Optimizer message: ',res_s.message)
    return 0

def d_opt():
    with open("Optimized_params_d.txt", mode="w") as file:
        x0 = [3.7008,4.7748,1.4956,2.2784]
        print("Starting with d")
        print("Initial value", bf.opt(x0))
        res_s = OPT.minimize(bf.opt, args=(2,nist_e_dav,R0, kk, Alpha, Zion),x0=s_params)
        file.write("d params  ")
        file.write(("%.8f   "*len(res_s.x))%tuple(res_s.x))
        file.write("\n")
        file.flush()
        print("D done, check file")
        print("Paramaters", res_s.x)
        print("Success on the minimization: ",res_s.success)
        print('Minimum value: ',res_s.fun)
        print('Optimizer message: ',res_s.message)
    return 0

def f_opt():
    with open("Optimized_params_f.txt", mode="w") as file:
        x0 = [3.8125,5.0332,2.1016,1.2707]
        print("Starting with f")
        print("Initial value", bf.opt(x0))
        res_s = OPT.minimize(bf.opt, args=(3,nist_e_fav,R0, kk, Alpha, Zion),x0=s_params)
        file.write("f params  ")
        file.write(("%.8f   "*len(res_s.x))%tuple(res_s.x))
        file.write("\n")
        file.flush()
        print("F done, check file")
        print("Paramaters", res_s.x)
        print("Success on the minimization: ",res_s.success)
        print('Minimum value: ',res_s.fun)
        print('Optimizer message: ',res_s.message)
