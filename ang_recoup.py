import numpy as np
import sympy.physics.wigner as wg
# Miscellaneous functions

def kij(i,j):
    if(i!=j):
        return 0
    else:
        return 1

def bJcs_kLS(n, lc, sc, jc, le, Jcs, se, L, S, J):
        UFT = np.zeros((n,n))
        for ii in range(n):
            for jj in range(n):
                deltas = kij(lc[ii],lc[jj])*kij(le[ii],le[jj])
                if(deltas==0):
                    UFT[ii,jj] = 0
                else:
                    #print("Values of ang momentum %.1f %.1f %.1f %.1f %.1f %.1f"%(lc(ii), le(ii), L(jj), Jcs(ii), sc(ii), jc(ii))
                    sj1 = wg.wigner_6j(sc[ii], se[ii], S[jj], Jcs[ii], lc[ii], jc[ii])
                    sj2 = wg.wigner_6j(lc[ii], le[ii], L[jj], J, S[jj], Jcs[ii])
                    UFT[ii,jj] = (-1)**(se[ii] + jc[ii] + le[ii] -J)*\
                    np.sqrt((1 + 2*jc[ii])*(1 + 2*Jcs[ii])*(1 + 2*L[jj])*(1 + 2*S[jj]))*sj1*sj2
        
        return UFT

    

def bjj_kJk(n, jc, le, se, je, K, J):
    UFT = np.zeros((n,n))
    for ii in range(n):
        for jj in range(n):
            deltas = kij(jc[ii],jc[jj])*kij(le[ii],le[jj])
            if(deltas==0):
                UFT[ii,jj] = 0
            else:
                #print("Values of ang momentum %.1f %.1f %.1f %.1f %.1f %.1f"%(lc(ii), le(ii), L(jj), Jcs(ii), sc(ii), jc(ii))
                sj1 = wg.wigner_6j(le[ii], se[ii], je[ii], J, jc[ii], K[jj])
                UFT[ii,jj] = (-1)**(je[ii] + jc[ii]-J)*\
                np.sqrt((1 + 2*je[ii])*(1 + 2*K[jj]))*sj1
    return UFT