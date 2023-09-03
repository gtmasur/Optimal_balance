# GTM: Updated version of Marcel Oliver's lsgtools.py and 
# randomdata.py codes, and here, the code is nondimensional.

"""
GTM: This file produces all initializations.

For L1-balanced model, we solve eqns (3.15) and (3.18) 
or (3.20) in MO's work: 'Comparison of variational 
balanced models for the rotating shallow water equations'.
(3:15)  : kinematic equation of u
(3:18)  : pv conservation to obtain pv from random h
(3:20)  : inversion of (3.18) to obtain h from random q
Therefore, it is possible to set up L1-balanced model
starting from random pv or h fields. Since \lambda = 1/2
in our case, then physical and balanced-model coordinates
are the same, no need any transformation.
In the code,
N   : Number of grid points in each direction
h0  : mean height
h   : free surface field
eps : Rossby number from simulation
d   : order of spectrum decay
"""

from pylab import *
from numpy.fft import *
import os, sys

#from Linear_BC import *

rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


tol = 1e-10  # to use in Helmholtz solvers

#================================================
# creates random field

def k (ix,iy,sim):
    return int(sqrt(-sim.sLap[ix,iy]))

def spectrum(u,sim):
    s = zeros(N)
    uft = fft2(u)
    for i in range(N):
        for j in range(N): 
            s[k(i,j,sim)] += abs(uft[i,j])**2 
    return s*(2*pi/N)**2  
    
def randfield (sim, s=4258463, ham=0.2):
    #Produce random field which an asymptotic spectral decay of
    #   k^-(d) and a spectral maximum at ak0
    global N
    d = sim.d
    ak0 = sim.ak0
    N = sim.Nx
    
    seed(s) # Seed random number generator
    ur = rand(N,N)
    sr = spectrum(ur,sim)
    
    uft = fft2(ur)
    for i in range(N):
        for j in range(N):
            uft[i,j] /= (sqrt(sr[k(i,j,sim)]) + 1e-12) 
    
    b = (7+d)/4.0     # Controls the scaling in the decaying part of the spectrum
    a = (2*b-3)/3.0   # Ensures that spectrum peaks at k=ak0
    ss = (-sim.sLap)**(3.0/2.0)/(-sim.sLap+a*ak0**2)**b
    ff = ifft2(ss*uft*sim.sFilt).real
    ff *= ham/max( abs(amax(ff)),abs(amin(ff)))  # this gives |ff|<0.2
    #plotfield(ff)
    return  ff
    
def plotfield(f):
    imshow(f, vmin = f.min(), vmax = f.max(), cmap = 'jet')
    #imshow(f, vmin = -0.2, vmax = 0.2, cmap = 'seismic')
    colorbar()
    title('The field in physical space')
    show()    
    
#================================================
# L1-balanced initialization settings...
"""
main functions: 
L1_initialdata : creates L1-balanced init. from random h
"""
       
def set_L1u(sim):
    normref = 1
    tol = 1e-10
    #N = sim.Nx  # assuming Nx=Ny
    #d = sim.d
    h0 = sim.h0
    eps = sim.eps
    
    for ii in range(sim.Nz):  
    
        # take desire h field...
        if sim.system == "Real":
            h = sim.soln.h[:,:,ii]
        elif sim.system == "Ramped":
            h = sim.ramp_soln.h[:,:,ii]
            
        # Be careful! In (3.15), everything is in h:total height.
        (hx,hy) = ( sim.Dx(h), sim.Dy(h))
        g = (eps/2)*(2*(h0+h)*sim.Lap(h)+ hx**2 + hy**2)
        ug = sim.Dy(-h + g)
        vg = sim.Dx(h - g)
        (u,v) = (-hy,hx)
        sim.sIHelm = sim.sFilt/(1-eps*h0*sim.sLap)
        while normref > tol:
            unew = sim.IHelm(ug + eps*( h*sim.Lap(u) + 2*hx*sim.Dx(u) + 2*hy*sim.Dy(u)))
            vnew = sim.IHelm(vg + eps*( h*sim.Lap(v) + 2*hx*sim.Dx(v) + 2*hy*sim.Dy(v)))
            normu = norm(u-unew)/norm(unew)
            normv = norm(v-vnew)/norm(vnew)
            normref = max(normu, normv)
            (u,v) = (unew,vnew)
        print("L1-u field is constructed with precision:", normref)
        
        # project on predefined field
        if sim.system == "Real":
            sim.soln.u[:,:,ii] = u[:,:] 
            sim.soln.v[:,:,ii] = v[:,:]   
        elif sim.system == "Ramped":
            sim.ramp_soln.u[:,:,ii] = u[:,:] 
            sim.ramp_soln.v[:,:,ii] = v[:,:]   
        
       
def h_from_randpv(sim):
    # solve (q-eps*Lap)h=1-qh0
    # to get h from random q field
    
    # as h is free surface height, do not add h0 
    normref = 1
    for ii in range(sim.Nz):
        q = sim.solnp.q[:,:,ii]
        qbar = mean(q)
        sim.sIHelm = sim.sFilt/(qbar-sim.eps*sim.sLap )
        h = zeros(q.shape) 
        while normref > tol:
            hnew = sim.IHelm(1.- q*sim.h0 - (q-qbar)*h)
            normref = norm(h-hnew)/norm(hnew)    
            h = hnew  
        print ("h field is constructed with precision:", normref)    
        sim.soln.h[:,:,ii] = h
    
def L1_initialdata(sim):
    # To set L1-balanced initialization 
    # 1. set random h field
    # 2. set L1-balanced u field
    for ii in range(sim.Nz):
        sim.soln.h[:,:,ii] = randfield(sim)  
    set_L1u(sim)
    

#================================================
"""
geostrophic condition for LBC or for initialdata
"""

def set_geostrophy(sim):

    for ii in range(sim.Nz):  
        if sim.system == 'Real':
            # compute v = dx_h
            # compute u = - dy_h
            h = sim.soln.h[:,:,ii].copy()
            sim.soln.v[:,:,ii] = sim.Dx(h)
            sim.soln.u[:,:,ii] = - sim.Dy(h)
            
        elif sim.system == 'Ramped':
            # compute v = dx_h
            # compute u = -dy_h
            h = sim.ramp_soln.h[:,:,ii].copy()
            sim.ramp_soln.v[:,:,ii] = sim.Dx(h)
            sim.ramp_soln.u[:,:,ii] = - sim.Dy(h)
    

def geost_initialdata(sim):
    for ii in range(sim.Nz):
        sim.soln.h[:,:,ii] = randfield(sim) 
    set_geostrophy(sim)
    

#================================================
"""
zero velocity depending on base point
"""

def zero_initialdata(sim):
    for ii in range(sim.Nz):
        sim.soln.h[:,:,ii] = randfield(sim) 
    # velocity is set as zero as default

#================================================
"jet is chosen as initial condition"

def jet_initialdata(sim):
    # insert Carsten's jet initial data in sim.soln
    path = '../InitialData/balanced_jet.npz'
    data = np.load(path)
    
    for ii in range(sim.Nz):
        sim.soln.u[:,:,ii] = data.f.u
        sim.soln.v[:,:,ii] = data.f.v
        sim.soln.h[:,:,ii] = data.f.h
        
    """
    # the code below is used to balanced Carsten's jet...
    
    sim.linear_bc = "obliqpr"
    build_proj_mats(sim)
    for ii in range(sim.Nz):
        u = data.f.u
        v = data.f.v
        h = data.f.h
        
        # project (u,v,h) to phase space
        zz = 1j*zeros((sim.Nx,sim.Ny,3,1))
        zz[:,:,0,0] = fft2(u)
        zz[:,:,1,0] = fft2(v)
        zz[:,:,2,0] = fft2(h)
            
        # separate into components
        zzr = 1j*zeros((sim.Nx,sim.Ny,3,1))
        zzg = 1j*zeros((sim.Nx,sim.Ny,3,1))
        for i in range(sim.Nx):
            for j in range(sim.Ny):
                zzr[i,j,:,:] = sim.Pr[i,j,:,:].dot(zz[i,j,:,:])   # zzr = Pr*zz
                zzg[i,j,:,:] = sim.Pg[i,j,:,:].dot(zz[i,j,:,:])   # zzg = Pg*zz
        
        if not allclose(zz, zzr + zzg):
            print("Error: Components do not build the phase flow!")
            sys.exit()
            
        # project on physical space
        ur = ifft2(zzr[:,:,0,0]).real    
        vr = ifft2(zzr[:,:,1,0]).real
        hr = ifft2(zzr[:,:,2,0]).real
        ug = ifft2(zzg[:,:,0,0]).real
        vg = ifft2(zzg[:,:,1,0]).real
        hg = ifft2(zzg[:,:,2,0]).real
        # check the relative norms of wave components
        #print absrel_norm(sim,ur+ug,u), absrel_norm(sim,vr+vg,v), absrel_norm(sim,hr+hg,h)
        
        # project on simulation object
        print "Rossby-wave component is projected on sim..."
        sim.soln.u[:,:,ii] = ur[:,:] 
        sim.soln.v[:,:,ii] = vr[:,:] 
        sim.soln.h[:,:,ii] = hr[:,:] 
        
        
        plotfield(sim.soln.u[:,:,ii])
        plotfield(sim.soln.v[:,:,ii])
        plotfield(sim.soln.h[:,:,ii])
        
        # Save data
        fname = './balanced_jet'
        np.savez_compressed(fname, u = ur[:,:], v = vr[:,:], h = hr[:,:], )"""
        
        
    print ("The simulation is initialised with Carsten's jet...")   
    

#def jet_initialdata(sim):
#    for ii in range(sim.Nz):
#        # set initial data for u
#        
#        a = int(sim.Nx/4)
#        b = int(3*sim.Nx/4)
#        u0 = 2.25
#        sim.soln.u[:,:,ii] = u0*( np.exp(-(sim.Y-sim.y[a])**2/(0.1)**2) 
#                               - np.exp(-(sim.Y-sim.y[b])**2/(0.1)**2) )
#        #sim.soln.v[:,:,ii] = 0*sim.soln.u[:,:,ii]
#         
#        #plotfield(sim.soln.u[:,:,ii])        
#                
#        # set geostrophic h for u
#        #def EE_old(y,k): return 1. + special.erf( 10.*sqrt(2)*(y-k) ) 
#        #sim.soln.h[:,:,ii] = cc*( EE(2*pi,a)-EE(0,a)+ EE(2*pi,b)-EE(0,b) 
#        
#        sim.soln.h[:,:,ii] = u0*sim.dx[0]*2*1e-3*sin(sim.X/sim.Lx*10*pi)
#        plotfield(sim.soln.h[:,:,ii])
#        
#        exit()

#================================================
"""
Main function: initialize_fields
depending on initialsetup, initialize fields.
"""

# balanced_wave: initialize the present run with 
# balanced state of the previous test run
def initialize_fields(sim):

    naming = True
    print ("initial setup:....", sim.initialsetup)
    # set random h field and then depending on initialsetup, set u
    if sim.initialsetup == 'geost':
        print ("Setting geostrophic velocity for random h...")
        geost_initialdata(sim)
    elif sim.initialsetup == 'L1':
        print ("Setting L1-balanced velocity for random h...")
        L1_initialdata(sim)
    elif sim.initialsetup == 'zero':
        print ("Setting zero velocity for random h...")
        zero_initialdata(sim)
        #elif sim.initialsetup == 'balanced_wave':
    elif sim.initialsetup == 'geost_balwave':
        naming = False
        print ("Setting previous balanced state as initial...")
        for ii in range(sim.Nz):
            sim.soln.u[:,:,ii] = sim.balu
            sim.soln.v[:,:,ii] = sim.balv
            sim.soln.h[:,:,ii] = sim.balh
    elif sim.initialsetup == 'jet':
        # here, we create jet flow
        jet_initialdata(sim)
            
    # After first normal initialization, switch to balanced_wave	
    if naming:
        if sim.optimal_balance:
            if sim.save_diagimb and sim.balanced_wave:
                sim.initialsetup = sim.initialsetup+'_balwave'
    
def set_conditions(sim):
    # define tools to use for boundary conditions
    sim.set_geostrophy = set_geostrophy
    sim.set_L1u = set_L1u
    
    

