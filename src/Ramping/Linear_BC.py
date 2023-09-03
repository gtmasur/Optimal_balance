"""
GTM: This file sets BC for ramping procedure at
     the linear end (t_ramp=0)
"""

from pylab import *
import sys
from scipy.fftpack import fftfreq, fft2, ifft2
from .Nonlinear_BC import *


def absrel_norm(sim,f,g):
    absnorm = norm(f-g)
    if sim.norm_type == 'absolute':
        return absnorm
    elif sim.norm_type == 'relative':
        # symetrized relative norm
        return absnorm/(0.5*norm(f) + 0.5*norm(g)) 
    
    
#############################################
# Functions to test projection matrices    

def Pr_obliq(sim,k,l):
    # build oblique projection matrix with B
    h0 = sim.h0
    eps = sim.eps
    def B_mat(k,l):
        return array([[ -1j*l,   0.,    -1j/eps],
                         [  1j*k,   1j/eps,   0.  ],
                         [   1.,    -h0*k,   -h0*l]])
                         
    def Binv_mat(k,l):
        #det = linalg.det(B_mat(k,l))
        c = -h0*(k**2+l**2)/eps - 1./eps**2
        return 1./c*array([[ -1j*h0*l/eps,       1j*h0*k/eps,     -1./eps**2],
                           [   1j*h0*k*l,      1j*h0*l**2 + 1j/eps,   k/eps ],
                           [-1j*h0*k**2-1j/eps,  -1j*h0*k*l,          l/eps ]])   

    M = zeros((3,3)); M[0,0] = 1.
    B = B_mat(k,l)
    B_inv = Binv_mat(k,l)
    #I = eye(3); R = B.dot(B_inv)
    #if not allclose(R,I): print "B matrix is wrong"
    return B.dot(dot(M,B_inv))   

def idemp_proj(sim, Pr, Pg, I):
    # check whether Pr is idempotent or not
    PP = zeros((sim.Nx,sim.Ny,3,3), dtype='complex128')
    for i in range(sim.Nx):
        for j in range(sim.Ny):
            PP[i,j,:,:] = dot(Pr[i,j,:,:], Pg[i,j,:,:])       
    print( allclose(PP, 0*I), norm(PP) )
    
def Pr_orth(k,l):
    # build orthogonal projection matrix
    n1 = array([-1j*l, 1j*k, 1])[:,newaxis]
    c = norm(n1)**2
    return 1./c * n1.dot(conjugate(n1).T)
    
#############################################    
   
   
def obliq_proj_mats(sim):
    # main idea:
    # for every physical point index (i,j), take their
    # corresponding wave number and build 3x3 matrix.
    #print( "....building oblique projection matrices....")
    h0 = sim.h0
    eps = sim.eps
    kk = (-1j*sim.sDx).real
    ll = (-1j*sim.sDy).real
    # build Pr: rossby wave proj. matrix
    #Ph = 1j*zeros((sim.Nx,sim.Ny,3,3))
    Pr = 1j*zeros((sim.Nx,sim.Ny,3,3))
    I = zeros((sim.Nx,sim.Ny,3,3))
    for i in range(sim.Nx):
        for j in range(sim.Ny):
            k = kk[i,j]; l = ll[i,j]
            I[i,j,:,:] = eye(3)
            #Ph[i,j,:,:] = Pr_obliq(sim,k,l)
            c = 1./ ( -h0*(k**2+l**2) -1./eps)  # =1/det*1/eps
            Pr[i,j,:,:] = c * array([[-h0*l**2,   h0*k*l,    1j*l/eps],
                                     [ h0*k*l,   -h0*k**2,  -1j*k/eps],
                                     [-1j*h0*l,   1j*h0*k,   -1./eps]])
                                      
    # build Pg: gravity wave proj. matrix        
    Pg = I - Pr  
    # obliq projection is idempotent: P^2=P
    #idemp_proj(sim, Pr, Pg, I) 
    #print( allclose(Pr,Ph) )
    return (Pr, Pg)

    
def orth_proj_mats(sim):
    #print( "....building orthogonal projection matrices....")
    kk = (-1j*sim.sDx).real   
    ll = (-1j*sim.sDy).real
    #Ph = 1j*zeros((sim.Nx,sim.Ny,3,3))
    Pr = 1j*zeros((sim.Nx,sim.Ny,3,3))
    I = zeros((sim.Nx,sim.Ny,3,3))
    for i in range(sim.Nx):
        for j in range(sim.Ny):
            k = kk[i,j];  l = ll[i,j]
            #PP = Pr_orth(k,l) 
            I[i,j,:,:] = eye(3)
            c =  1./(l**2 + k**2 + 1.)  # =1/||n1||
            Pr[i,j,:,:] = c * array([[ l**2,   -k*l,  -1j*l ],
                                     [ -k*l,   k**2,   1j*k ],  
                                     [ 1j*l,   -1j*k,   1. ]])
            #Ph[i,j,:,:] = transpose(conjugate(Pr[i,j,:,:]))
    Pg = I - Pr
    # orthogonal projections are idempotent: P^2=P
    #idemp_proj(sim, Pr, Pg, I)
    # orthogonal projections are Hermitian: P = conj(P^T)
    #print( allclose(Pr,Ph))
    return (Pr, Pg)
    
    
def build_proj_mats(sim):  
    # build and keep proj matrices over one test case
    if sim.linear_bc == "obliqpr":
        (sim.Pr, sim.Pg) = obliq_proj_mats(sim)
    elif sim.linear_bc == "orthpr":
        (sim.Pr, sim.Pg) = orth_proj_mats(sim)  
       
    
def wave_separation(sim):

    for ii in range(sim.Nz):
    
        u = sim.ramp_soln.u[:,:,ii]
        v = sim.ramp_soln.v[:,:,ii]
        h = sim.ramp_soln.h[:,:,ii]
        
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
            print( "Error: Components do not build the phase flow!")
            sys.exit()
            
        # project on physical space
        ur = ifft2(zzr[:,:,0,0]).real    
        vr = ifft2(zzr[:,:,1,0]).real
        hr = ifft2(zzr[:,:,2,0]).real
        ug = ifft2(zzg[:,:,0,0]).real
        vg = ifft2(zzg[:,:,1,0]).real
        hg = ifft2(zzg[:,:,2,0]).real
        # check the relative norms of wave components
        #print( absrel_norm(sim,ur+ug,u), absrel_norm(sim,vr+vg,v), absrel_norm(sim,hr+hg,h))
        
        # project on simulation object
        if sim.wave == 'Rossby':
        
            if sim.optbal_algorithm in ['algorithm1','algorithm3']:
                print( "Rossby-wave component is projected on sim...")
                sim.ramp_soln.u[:,:,ii] = ur[:,:] 
                sim.ramp_soln.v[:,:,ii] = vr[:,:] 
                sim.ramp_soln.h[:,:,ii] = hr[:,:]
                
            elif sim.optbal_algorithm == 'algorithm2':
                if sim.end_time > sim.time: 
                    print( "Rossby-wave component is projected on sim...")
                    sim.ramp_soln.u[:,:,ii] = ur[:,:] 
                    sim.ramp_soln.v[:,:,ii] = vr[:,:] 
                    sim.ramp_soln.h[:,:,ii] = hr[:,:] 
                elif sim.end_time == sim.time:
                    print( "norms are computed at linear end...")
                    compute_norm_linearend(sim,u,v,h,ur,vr,hr)
                    
        elif sim.wave == 'Gravity':
            # This is available just for test case....        
            print( "gravity-wave component is projected on sim...")
            sim.ramp_soln.u[:,:,ii] = ug[:,:] 
            sim.ramp_soln.v[:,:,ii] = vg[:,:] 
            sim.ramp_soln.h[:,:,ii] = hg[:,:] 


def preserve_zeta(sim):
    # compute h as stream function
    for ii in range(sim.Nz):
        u = sim.ramp_soln.u[:,:,ii]
        v = sim.ramp_soln.v[:,:,ii]
        # zeta: vorticity
        zeta = sim.Div_perp(u,v)
        # h as a stream function
        sim.ramp_soln.h[:,:,ii] = sim.ILap(zeta)
    # compute u via geostrophic condition
    sim.set_geostrophy(sim) 
    return
 

 
def compute_norm_linearend(sim,u,v,h,ur,vr,hr):

    def prognostic_alg2(sim,u,v,h):
        # short version of classic_to_prognostic
        zeta = sim.Div_perp(u,v)
        q = sim.Filt((sim.eps*zeta + 1.)/(sim.h0+h))
        delta = sim.Div(u,v)
        gamma = zeta - sim.Lap(h)
        return (q,delta,gamma)

    # added (delta,gamma) norms in algorithm2....
    for ii in range(sim.Nz): 
        if sim.prognostic:
            q   = sim.ramp_solnp.q[:,:,ii]       
            delta = sim.ramp_solnp.delta[:,:,ii] 
            gamma = sim.ramp_solnp.gamma[:,:,ii]  
        else:            
            (q, delta, gamma) = prognostic_alg2(sim,u,v,h)
        
    (ug, vg, hg) = (u-ur, v-vr, h-hr) 
    (qg, deltag, gammag) = prognostic_alg2(sim,ug,vg,hg)
    (qr, deltar, gammar) = prognostic_alg2(sim,ur,vr,hr)

    sim.norm_u = absrel_norm(sim,ur,u)
    sim.norm_v = absrel_norm(sim,vr,v)
    sim.norm_h = absrel_norm(sim,hr,h)
    sim.norm_ug = norm(ug)
    sim.norm_vg = norm(vg)
    sim.norm_hg = norm(hg)
    sim.norm_ur = norm(ur)
    sim.norm_vr = norm(vr)
    sim.norm_hr = norm(hr)
    
    sim.norm_q = absrel_norm(sim,qr,q)
    sim.norm_delta = absrel_norm(sim,deltar,delta)
    sim.norm_gamma = absrel_norm(sim,gammar,gamma)
    sim.norm_qg = norm(qg)
    sim.norm_deltag = norm(deltag)
    sim.norm_gammag = norm(gammag)
    sim.norm_qr = norm(qr)
    sim.norm_deltar = norm(deltar)
    sim.norm_gammar = norm(gammar)
    
    
def linearbc_and_norm(sim):
	# a small trick to use available functions
	
	# project linear wave on physical space
	sim.project_on_physicalspace(sim)
	# find the rossby wave
	# (Not sure if "preserved zeta" will work
	# but we cannot put geostrophy since hg becomes 0)
	preserve_zeta(sim)
	# compute norms 
	for ii in range(sim.Nz):
	    (h,hr) = (sim.soln.h[:,:,ii],sim.ramp_soln.h[:,:,ii])
	    (u,ur) = (sim.soln.u[:,:,ii],sim.ramp_soln.u[:,:,ii])
	    (v,vr) = (sim.soln.v[:,:,ii],sim.ramp_soln.v[:,:,ii])
	    print( "norms are computed at linear end...")
	    compute_norm_linearend(sim,u,v,h,ur,vr,hr)
    

def set_linear_bc(sim):
    # For obliqpr and orthpr, algorithm2 is included inside.
    
    if sim.linear_bc == "obliqpr":
    
        print( "")
        print( "setting linear_bc for oblique projection...")
        #sim.wave = 'Rossby' 
        wave_separation(sim)        
#        check_pv_conser(sim)
            
    elif sim.linear_bc == "orthpr":
  
        print( "")
        print( "setting linear_bc for orthogonal projection...")
        #sim.wave = 'Rossby'
        wave_separation(sim)
        
    # TODO: change the name of altgeost...    
    elif sim.linear_bc == "altgeost":
    
        print( "")
        print( "setting linear_bc as preserved zeta...")
          
        if sim.optbal_algorithm in ['algorithm1','algorithm3']:
            preserve_zeta(sim)  
        elif sim.optbal_algorithm == 'algorithm2':
            if sim.end_time > sim.time: 
                preserve_zeta(sim)  
            elif sim.end_time == sim.time:
                linearbc_and_norm(sim)

    elif sim.linear_bc == "geost":
        print( "")
        print( "setting linear_bc as geost (preserved h)....")
        sim.set_geostrophy(sim) 
        
    return 

#############################################
# Work on this part while checking the paper....
def check_pv_conser(sim):
    # to observe Pv conservation data... 
    # it is assumed that prognostic is on
    qbef = sim.ramp_solnp.q.copy()
    sim.classic_to_prognostic(sim)
    qaf = sim.ramp_solnp.q.copy()
    print( "")
    print( "***** check pv conservation of LBC:")
    print( norm(qbef - qaf)/norm(qbef), allclose(qbef, qaf))
   
def check_pv_conser2(sim):
    # to observe Pv conservation data... 
    # it is assumed that prognostic is on
    def qq(sim):
        h = sim.soln.h[:,:,0].copy()
        u = sim.soln.u[:,:,0].copy()
        v = sim.soln.v[:,:,0].copy()
        zeta = sim.Div_perp(u,v)
        # linear PV 
        print( h)
        q = (sim.eps*zeta + h/sim.h0)
        return q
    """qbef = sim.ramp_solnp.q.copy()
    sim.classic_to_prognostic(sim)
    qaf = sim.ramp_solnp.q.copy()"""
    qbef = qq(sim)
    sim.classic_to_prognostic(sim)
    qaf = qq(sim)
    print( "")
    print( "***** check pv conservation of LBC:")
    print( norm(qbef - qaf), allclose(qbef, qaf))
    
#############################################



