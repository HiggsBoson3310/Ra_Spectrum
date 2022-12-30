import numpy as np
from scipy.interpolate import BSpline
import scipy.integrate as integ
import scipy.linalg as LA
import matplotlib.pyplot as plt

#Spectrographic notation
angs = {"s":0, "p":1, "d":2, "f":3}
#Fine structure constant
fsc = 1/137.035999084

def v(r, al1, al2, al3, rc, alpha, zion, l):
    f1 = np.exp(-al1*r)*(zion-2)
    f2 = np.exp(-al3*r)*al2
    z = 2 + f1 + r*f2
    expf = np.exp(-(r/rc)**6)
    upol = (-0.5*alpha/r**4)*(1-expf)
    v = upol-z/r +0.5 * l*(l+1)/(r*r)
    return v

def vso(r, al1, al2, al3, rc, alpha, zion, l, j):
    f1 = np.exp(-al1*r)*(zion-2)
    f2 = np.exp(-al3*r)*al2
    z = 2 + f1 + r*f2
    dz = -al1*f1 - al3*f2*r + f2
    expf = np.exp(-(r/rc)**6)
    upol = (-0.5*alpha/r**4)*(1-expf)
    v = upol-z/r 
    dv = -1/r * (dz-z/r) + 2*alpha/r**5 * (1-expf*(1+3.0/2.0*(r/rc)**6))
    vs = fsc*fsc/2.0 * (j*(j+1)-l*(l+1)-3.0/4.0)/2.0
    vs = vs*1/r*dv
    vs = vs/(1-fsc*fsc/2.0*v)/(1-fsc*fsc/2.0*v)
    vt = v + vs + 0.5 * l*(l+1)/(r*r)
    return vt


def exp_mesh(r0,rm, g, N):
    rl = np.zeros(N)
    for i in range(N):
        rl[i] = r0 + (rm-r0)*(np.exp(g*((i)/(N-1)))-1)/(np.exp(g)-1)
    return rl

def bspline_h(r0, N, k, l):
    #Generate the Bspline basis and store it in a list
    r = exp_mesh(0, r0, 20, N+k-1)
    B = []
    Bm  = []
    BM = []
    #print(len(r))
    fig, ax = plt.subplots(1,1, figsize=(20,20))
    for i in range(0,len(r)-k):
        B.append(BSpline.basis_element(r[i:i+k+1],extrapolate=False))
        Bm.append(r[i])
        if((i+k)<len(r)):
            BM.append(r[i+k])
        else:
            BM.append(r[-1])

        if(i>0): 
            rfine = np.linspace(Bm[-1],BM[-1],200)
            ax.plot(rfine,B[-1](rfine))
    plt.savefig("Splines_Hydrogen.png",dpi=150)
    #Compute the kinetic energy term
    H = np.zeros((len(B),len(B)))
    S = np.zeros((len(B),len(B)))
    #print(len(B))
    for i in range(len(B)):
        for j in range(i+1):
            a =  max(Bm[i],Bm[j])
            b =  min(BM[i],BM[j])
            bi = B[i]
            bj = B[j]
            if(a<b):
                fh = lambda x: -0.5*bi(x)*bj.derivative(nu=2)(x) + 0.5*l*(l+1) *  bi(x)*bj(x)/x**2 - bi(x)*bj(x)/x
                H[i,j] = integ.quad(fh ,a,b)[0]
                H[j,i] = H[i,j]
                S[i,j] = integ.quad(lambda x: bi(x)*bj(x), a, b)[0]
                S[j,i] = S[i,j]
    eig1 = LA.eigvals(a = H, b = S, check_finite=False)
    #eig = LA.eigh(a = H, b = S, lower=True, eigvals_only=True, subset_by_value=[-np.inf,0.0], overwrite_a=True, overwrite_b=True, check_finite=False)
    return eig1

def Bspline_spec(al1, al2, al3, rc, l, e0, Nn, R0, kk, Alpha, Zion,plot_op=False):
    #Generate the Bspline basis and store it in a list
    r = exp_mesh(1e-12, R0, 20, Nn+kk)
    #r = np.linspace(1e-12, R0, Nn+kk+1)
    B = []
    Bm  = []
    BM = []
    if(plot_op): fig, ax = plt.subplots(1,1, figsize=(20,20))
    for i in range(0,len(r)-kk):
        B.append(BSpline.basis_element(r[i:i+kk+1],extrapolate=False))
        Bm.append(r[i])
        BM.append(r[i+kk])
        if(plot_op):
            if(i>0): 
                rfine = np.linspace(Bm[-1],BM[-1],200)
                ax.plot(rfine,B[-1](rfine))
            plt.savefig("Splines_YbIon.png",dpi=150)
    if(plot_op): plt.close()
    #Compute the kinetic energy term
    H = np.zeros((len(B),len(B)))
    S = np.zeros((len(B),len(B)))
    #print(len(r), len(B)+kk+1)
    #print(len(B))
    lxs, ws = np.polynomial.legendre.leggauss(2*kk-2)
    lxh, wh = np.polynomial.legendre.leggauss(2*kk-4)
    for i in range(len(B)):
        for j in range(i+1):
            a =  max(Bm[i],Bm[j])
            b =  min(BM[i],BM[j])
            bi = B[i]
            bj = B[j]
            if(a<b):
                #il = np.argwhere(r==a)[0]
                #ih = np.argwhere(r==b)[0]
                #while(il+1<=ih):
                #    lxsi = 0.5*(lxs+1)*(r[il+1]-r[il])+r[il]
                #   S[i,j] += sum(ws*bi(lxsi)*bj(lxsi)) * 0.5*(r[il+1]-r[il])
                #   lxhi = 0.5*(lxh+1)*(r[il+1]-r[il])+r[il]
                #   H[i,j] += -0.5*sum(wh*bi(lxhi)*bj.derivative(nu=2)(lxhi)) * 0.5*(r[il+1]-r[il])
                #   il+=1
                xx = np.linspace(a,b,200)
                fh = v(xx,al1,al2,al3,rc,Alpha,Zion,l)*bi(xx)*bj(xx)-0.5*bi(xx)*bj.derivative(nu=2)(xx)
                H[i,j] += integ.simpson(x = xx, y=fh)
                H[j,i] = H[i,j]
                S[i,j] = integ.simpson(x =xx, y=bi(xx)*bj(xx))
                S[j,i] = S[i,j]
    eig1, ev1 = LA.eig(a = H, b = S, check_finite=True)
    print(np.sort(eig1[np.where(eig1<=0)]))
    #eig = LA.eigh(a = 5000*H, b = S, lower=True, eigvals_only=True, subset_by_value=[-np.inf,0.0], overwrite_a=True, overwrite_b=True, check_finite=False)
    eig1, ev1 = eig1[np.intersect1d(np.where(eig1<0),np.where(eig1>e0))], ev1[:,np.intersect1d(np.where(eig1<0),np.where(eig1>e0))]
    #print(np.shape(ev1))
    xx = exp_mesh(1e-12, R0, 20, 1500)
    nodes = np.zeros(len(eig1))
    psi0 = np.zeros((len(eig1), len(xx)))
    for psi,vec in enumerate(np.transpose(ev1)):
        psiv = []
        po = 0
        for i,x in enumerate(xx):
            pn=0
            for b,xl,xM in zip(enumerate(B),Bm,BM):
                if(x<xM and x>xl):
                    pn+= vec[b[0]]*b[1](x)
            if((pn*po)<0 and eig1[psi]>v(x,al1,al2,al3,rc,Alpha,Zion,0)): 
                #print(x,po,pn)
                nodes[psi]+=1
            po=pn
            psiv.append(pn)
        psi0[psi]=np.array(psiv)
    if(plot_op):
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        for i,psi in enumerate(psi0):
            ax.plot(xx,psi,label="%i"%(i))
        ax.legend(loc=0)
        ax.set_xlim(0,45)
        #ax.set_ylim(-60e-12,60e-12)
        plt.savefig("First_energy_l%i.png"%(l))
        np.savetxt("First_state_l%i.dat"%(l), np.column_stack((xx,psi0[0])))
    return np.column_stack((np.real(eig1),nodes.astype(float)))

def Bspline_spec_so(al1, al2, al3, rc, l, j, e0, Nn, R0, kk, Alpha, Zion, plot_op=False):
    #Generate the Bspline basis and store it in a list
    r = exp_mesh(1e-12, R0, 20, Nn+kk)
    #r = np.linspace(1e-12, R0, Nn+kk+1)
    B = []
    Bm  = []
    BM = []
    if(plot_op): 
        fig, ax = plt.subplots(1,1, figsize=(20,20))
    for i in range(0,len(r)-kk):
        B.append(BSpline.basis_element(r[i:i+kk+1],extrapolate=False))
        Bm.append(r[i])
        BM.append(r[i+kk])
        if(plot_op):
            if(i>0): 
                rfine = np.linspace(Bm[-1],BM[-1],200)
                ax.plot(rfine,B[-1](rfine))
            plt.savefig("Splines_YbIon_so.png",dpi=100)
    if(plot_op):
        plt.close()
    #Compute the kinetic energy term
    H = np.zeros((len(B),len(B)))
    S = np.zeros((len(B),len(B)))
    #print(len(r), len(B)+kk+1)
    #print(len(B))
    lxs, ws = np.polynomial.legendre.leggauss(2*kk-2)
    lxh, wh = np.polynomial.legendre.leggauss(2*kk-4)
    for i in range(len(B)):
        for q in range(i+1):
            a =  max(Bm[i],Bm[q])
            b =  min(BM[i],BM[q])
            bi = B[i]
            bj = B[q]
            if(a<b):
                #il = np.argwhere(r==a)[0]
                #ih = np.argwhere(r==b)[0]
                #while(il+1<=ih):
                #    lxsi = 0.5*(lxs+1)*(r[il+1]-r[il])+r[il]
                #   S[i,j] += sum(ws*bi(lxsi)*bj(lxsi)) * 0.5*(r[il+1]-r[il])
                #   lxhi = 0.5*(lxh+1)*(r[il+1]-r[il])+r[il]
                #   H[i,j] += -0.5*sum(wh*bi(lxhi)*bj.derivative(nu=2)(lxhi)) * 0.5*(r[il+1]-r[il])
                #   il+=1
                xx = np.linspace(a,b,200)
                fh = vso(xx,al1,al2,al3,rc,Alpha,Zion,l,j)*bi(xx)*bj(xx)-0.5*bi(xx)*bj.derivative(nu=2)(xx)
                H[i,q] += integ.simpson(x = xx, y=fh)
                H[q,i] = H[i,q]
                S[i,q] = integ.simpson(x =xx, y=bi(xx)*bj(xx))
                S[q,i] = S[i,q]
    eig1, ev1 = LA.eig(a = H, b = S, check_finite=True)
    #eig = LA.eigh(a = 5000*H, b = S, lower=True, eigvals_only=True, subset_by_value=[-np.inf,0.0], overwrite_a=True, overwrite_b=True, check_finite=False)
    eig1, ev1 = eig1[np.intersect1d(np.where(eig1<0),np.where(eig1>e0))], ev1[:,np.intersect1d(np.where(eig1<0),np.where(eig1>e0))]
    #print(np.shape(ev1))
    xx = exp_mesh(1e-12, R0, 20, 1500)
    nodes = np.zeros(len(eig1))
    psi0 = np.zeros((len(eig1), len(xx)))
    for psi,vec in enumerate(np.transpose(ev1)):
        psiv = []
        po = 0
        for i,x in enumerate(xx):
            pn=0
            for b,xl,xM in zip(enumerate(B),Bm,BM):
                if(x<xM and x>xl):
                    pn+= vec[b[0]]*b[1](x)
            if((pn*po)<0 and eig1[psi]>v(x,al1,al2,al3,rc,Alpha,Zion,0)): 
                #print(x,po,pn)
                nodes[psi]+=1
            po=pn
            psiv.append(pn)
        psi0[psi]=np.array(psiv)
    if(plot_op):
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        for i,psi in enumerate(psi0):
            ax.plot(xx,psi,label="%i"%(i))
        ax.legend(loc=0)
        ax.set_xlim(0,45)
        #ax.set_ylim(-60e-12,60e-12)
        plt.savefig("First_energy_l%i_so.png"%(l))
        plt.close()
        np.savetxt("First_state_l%i_so.dat"%(l), np.column_stack((xx,psi0[0])))
    return np.column_stack((np.real(eig1),nodes.astype(float)))

def opt(p, l, nist_e,R0, kk, Alpha, Zion):
    al1, al2, al3, rc = p
    e0 = -0.5
    Nn = 650

    spline = Bspline_spec_so(al1, al2, al3, rc, l, l+1./2., e0, Nn,R0, kk, Alpha, Zion)
    for n in nist_e[:,1]:
        if(not(n in spline[:,1].astype(int))):
            #print("There is no sol with %i nodes"%(n))
            spline = np.append(spline,[[2,n]],axis=0)
    tot = 0
    #print(spline)
    for i,n in enumerate(nist_e[:,1]):
        tot += (nist_e[i,0]-spline[np.where(spline[:,1].astype(float)==n)[0],0])**2
    return tot

def opt_fs(p, l, nist_e_fs,R0, kk, Alpha, Zion):
    al1, al2, al3, rc = p
    e0 = -0.5
    Nn = 650
    if(l==0):
        print("No spin orbit splitting in s orbital use opt function.")
        return None
    spline_m = Bspline_spec_so(al1, al2, al3, rc, l, l-1./2., e0, Nn,R0, kk, Alpha, Zion)
    spline_p = Bspline_spec_so(al1, al2, al3, rc, l, l+1./2., e0, Nn,R0, kk, Alpha, Zion)
    for n in nist_e_fs[:,1]:
        if(not(n in spline_m[:,1].astype(int))):
            spline_m = np.append(spline_m,[[2,n]],axis=0)
            spline_p = np.append(spline_p,[[2,n]],axis=0)
    tot = 0
    for i,n in enumerate(nist_e_fs[:,2]):
        ipp = np.where(spline_p[:,1].astype(float)==n)[0]
        imm = np.where(spline_m[:,1].astype(float)==n)[0]
        tot += (nist_e_fs[i,0]-spline_p[ipp,0])**2 + (nist_e_fs[i,1]-spline_m[imm, 0])**2 + 10*((nist_e_fs[i,0]-nist_e_fs[i,1])**2 - (spline_p[ipp,0]-spline_m[imm,0])**2)**2
    return tot




'''
res = OPT.minimize(opt_fs, x0=[4.04272, 16.3226, 2.38868, 1.26463], args=(1, nist_e_pfs), method='Nelder-Mead', options={'disp':True, 'maxfev':200, 'return_all':True, 'adaptive': True})
print(res.x)


#print(nist_e_s)
#ts = time.time()
#res = OPT.minimize(opt, x0 = [4.04272, 16.3226, 2.38868, 1.26463],args = (0, nist_e_s), method='Nelder-Mead')
#tf = time.time()
print(nist_e_s
)
opt([4.04272, 16.3226, 2.38868, 1.26463],0, nist_e_s)
ts = time.time()
print(Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 1, -0.5))
te = time.time()
print(te-ts)

print(Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 1, -0.5, 320, plot_op=False))
print(pav)

print(Bspline_spec_so(4.04272, 16.3226, 2.38868, 1.26463, 1, 1./2., -0.5, 320, plot_op=False))
print(p12)

print(Bspline_spec_so(4.04272, 16.3226, 2.38868, 1.26463, 1, 3./2., -0.5, 320, plot_op=False))
print(p32)
te = time.time()


ts = time.time()         #al1,     al2,     al3,      rc, l,     j,   e0,  Nn, plot_op=False

print(640)
Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 0, -0.5, 650, plot_op=False)
print(320)
Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 0, -0.5, 1200, plot_op=False)
print(2500)
Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 0, -0.5, 2500, plot_op=False)
te= time.time()
print(te-ts)


ts= time.time()
print(nist_e_dav)
print(opt([4.04272000,    71.32260000,    2.38868000,    1.26463000], 2, nist_e_dav))
#res_d = OPT.minimize(opt, args=(2,nist_e_dav),x0=[4.04272000,    71.32260000,    2.38868000,    1.26463000])
#print(res_d.x)
#print(opt(, 2, nist_e_dav))
#print(Bspline_spec(4.04272, 16.3226, 2.38868, 1.26463, 0, -0.5,320))
te = time.time()
print(te-ts)



#print(tf-ts,res.x)

'''
