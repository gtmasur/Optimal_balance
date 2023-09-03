"""
GTM : This file sets BC for ramping procedure at
      the nonlinear end (t_ramp = T).
"""

from pylab import *
import numpy as np
import sys
from scipy.fftpack import fftfreq, fft2, ifft2


def absrel_norm(sim,f,g):
    absnorm = norm(f-g)
    if sim.norm_type == 'absolute':
        return absnorm
    elif sim.norm_type == 'relative':
        # symetrized relative norm
        return absnorm/(0.5*norm(f) + 0.5*norm(g))     
   

def IHelm_solver(sim, q, gamma):  
    norm_h = 1.
    tol = 1e-10
    # iterate until ||hh - h|| < thrs
    qbar = mean(q)
    qhat = q-qbar
    sim.sIHelm = sim.sFilt / (- qbar + sim.eps*sim.sLap )  # used in sim.IHelm
    h = np.zeros((sim.Nx,sim.Ny)) 
    while norm_h > tol:
        h_new = sim.IHelm(-sim.eps*gamma + sim.h0*q + h*qhat - 1.)
        norm_h = norm(h-h_new)/norm(h_new)    
        h = h_new   
        #print norm_h 
    return h
        
def prognostic_to_classic(sim, q, delta, gamma):
    # Aware that this can be written in another way, 
    # but this is the best in the context 
    # (q, delta, gamma) --> (u,v,h)   
    # 1. solve the eqn for h
    # 2. compute zeta with this h
    # 3. solve the eqns for u and v         
           
    for ii in range(sim.Nz): 
        if sim.system == "Real":
            ubar = sim.soln.u[:,:,ii].mean()
            vbar = sim.soln.v[:,:,ii].mean()
        elif sim.system == "Ramped":
            ubar = sim.ramp_soln.u[:,:,ii].mean()
            vbar = sim.ramp_soln.v[:,:,ii].mean()     
                 
        # find h via Helmholtz-type equation
        hhat = IHelm_solver(sim, q, gamma)  # mean free
        # compute zeta from q
        zeta = sim.Filt( (q*(sim.h0 + hhat) - 1.)/sim.eps )
        # find u,v by Helmholtz decomposition
        psi = sim.ILap(zeta)      # stream function
        phi = sim.ILap(delta)     # velocity potential
        uhat = sim.Div_perp(psi,phi) + ubar
        vhat = sim.Div(psi,phi) + vbar
        
    return (uhat, vhat, hhat)


def classic_to_prognostic(sim):
    # (u,v,h) --> (q, delta, gamma)
    
    for ii in range(sim.Nz): 
        if sim.system == "Real":
            u = sim.soln.u[:,:,ii]
            v = sim.soln.v[:,:,ii]
            h = sim.soln.h[:,:,ii]
        elif sim.system == "Ramped":
            u = sim.ramp_soln.u[:,:,ii]
            v = sim.ramp_soln.v[:,:,ii]
            h = sim.ramp_soln.h[:,:,ii]
        
        # relative vorticity zeta
        zeta = sim.Div_perp(u,v)
        # potential vorticity q: materially conserved
        q = sim.Filt((sim.eps*zeta + 1.)/(sim.h0+h))
        # velocity divergence delta
        delta = sim.Div(u,v)
        # ageostrophic vorticity gamma 
        gamma = zeta - sim.Lap(h)
        
        if sim.prognostic:
            if sim.system == "Real":
                sim.solnp.q[:,:,ii]     = q
                sim.solnp.delta[:,:,ii] = delta
                sim.solnp.gamma[:,:,ii] = gamma
            elif sim.system == "Ramped":
                sim.ramp_solnp.q[:,:,ii]     = q
                sim.ramp_solnp.delta[:,:,ii] = delta
                sim.ramp_solnp.gamma[:,:,ii] = gamma
        else:
            return  (q, delta, gamma)
            
#############################################
# Function to test PV-inversion equations

def prog_classic_test(sim, q, delta, gamma):
    # test PV and PV_inversion equations
    (u_inv, v_inv, h_inv) = prognostic_to_classic(sim, q, delta, gamma)
    if sim.system == "Real":
        print ("")
        print ("Comparison for real system:")
        print (absrel_norm(sim, sim.soln.u, u_inv[:,:,np.newaxis]))
        print (absrel_norm(sim, sim.soln.v, v_inv[:,:,np.newaxis]))
        print (absrel_norm(sim, sim.soln.h, h_inv[:,:,np.newaxis]))
    elif sim.system == "Ramped":
        print ("")
        print ("Comparison for ramped system:")
        print (absrel_norm(sim, sim.ramp_soln.u, u_inv[:,:,np.newaxis]))
        print (absrel_norm(sim, sim.ramp_soln.v, v_inv[:,:,np.newaxis]))
        print (absrel_norm(sim, sim.ramp_soln.h, h_inv[:,:,np.newaxis]))
    return (u_inv, v_inv, h_inv)
    
#############################################


def project_on_physicalspace(sim):
	
	for ii in range(sim.Nz):
		sim.soln.u[:,:,ii] = sim.ramp_soln.u[:,:,ii]
		sim.soln.v[:,:,ii] = sim.ramp_soln.v[:,:,ii]
		sim.soln.h[:,:,ii] = sim.ramp_soln.h[:,:,ii]
		# needed to save balanced states
		if sim.prognostic:
			sim.solnp.q[:,:,ii] = sim.ramp_solnp.q[:,:,ii]         
			sim.solnp.delta[:,:,ii] = sim.ramp_solnp.delta[:,:,ii]
			sim.solnp.gamma[:,:,ii] = sim.ramp_solnp.gamma[:,:,ii]


def compute_norm(sim):

	for ii in range(sim.Nz):
		if sim.prognostic:
			(q, qT) = (sim.solnp.q[:,:,ii],sim.ramp_solnp.q[:,:,ii])
			(delta, deltaT) = (sim.solnp.delta[:,:,ii],sim.ramp_solnp.delta[:,:,ii])
			(gamma, gammaT) = (sim.solnp.gamma[:,:,ii],sim.ramp_solnp.gamma[:,:,ii])
		else:            
			sim.system = "Real"
			(q, delta, gamma) = classic_to_prognostic(sim)
			sim.system = "Ramped"
			(qT, deltaT, gammaT) = classic_to_prognostic(sim)
		(h,hr) = (sim.soln.h[:,:,ii],sim.ramp_soln.h[:,:,ii])
		(u,ur) = (sim.soln.u[:,:,ii],sim.ramp_soln.u[:,:,ii])
		(v,vr) = (sim.soln.v[:,:,ii],sim.ramp_soln.v[:,:,ii])
		
		sim.norm_h = absrel_norm(sim, h, hr)
		sim.norm_u = absrel_norm(sim, u, ur)
		sim.norm_v = absrel_norm(sim, v, vr)
		sim.norm_uv = absrel_norm(sim, sqrt(u**2+v**2), sqrt(ur**2+vr**2))
		sim.norm_q = absrel_norm(sim, q, qT)
		sim.norm_delta = absrel_norm(sim, delta, deltaT)
		sim.norm_gamma = absrel_norm(sim, gamma, gammaT)

		# at the end of rebalancing, to save rebalanced state assign all variables          
		# soln_(u,v,h) = (uT,vT,hT) at t=t' (end of forward time)
		project_on_physicalspace(sim)


def basepoint_pv(sim):  
    # run the loop until the stopping creiateria is satisfied. 
    # system = 'Ramped' at last, so back to ramping procedure.
    
    for ii in range(sim.Nz): 
        # compute soln_(u,v,h) ---> (q, delta, gamma)
        # compute ramp_soln_(u,v,h)(T) ---> (q, delta, gamma)(T)
        if sim.prognostic:
            q      = sim.solnp.q[:,:,ii] 
            qTnn   = sim.ramp_solnp.q[:,:,ii]       
            deltaT = sim.ramp_solnp.delta[:,:,ii] 
            gammaT = sim.ramp_solnp.gamma[:,:,ii]      
        else:            
            sim.system = "Real"
            (q, delta, gamma) = classic_to_prognostic(sim)
            sim.system = "Ramped"
            (qTnn, deltaT, gammaT) = classic_to_prognostic(sim)
            #prog_classic_test(sim, qTnn, deltaT, gammaT)
        
        # make the sequence
        # qTn = qT(n), qTnn = qT(n+1)
        # decide stopping criteria
        sim.norm_q  = absrel_norm(sim, q, qTnn)
        sim.norm_qT = absrel_norm(sim, qTnn, sim.qTn) 
        
        if sim.norm_qT <= sim.kappa:
            # convergence criteria is satisfied
            # stop optbal iterations
            sim.optbal_iterate = False
            
            if sim.time < sim.end_time :   
                # first cycle 
                # do not disturb the balanced state by qT -> q
                # soln_(u,v,h) = (uT,vT,hT) to evolve in real time
                project_on_physicalspace(sim)
                
                # This balanced wave is the initial data for next test.
                # take (q,dT,gT) and find (u,v,h) to save for next test cases
                # keep initial q since q is base-point
                if sim.balanced_wave:
                	(iu, iv, ih) = prognostic_to_classic(sim, q, deltaT, gammaT)
                	print ("Saving the balanced wave for upcoming tests...")
                	sim.balu = iu
                	sim.balv = iv
                	sim.balh = ih
                	
            elif sim.end_time == sim.time :  
                # second cycle
                # compute diagnosed imbalance, then attain ramp_soln to soln
                compute_norm(sim)
            
        elif sim.norm_qT > sim.kappa:
            # convergence criteria is not satisfied
            # sim.optbal_iterate = True by default
             
            # set Nonlinear BC q(T) = q and derive the new state
            (iu, iv, ih) = prognostic_to_classic(sim, q, deltaT, gammaT)
             
            # (u,v,h)_inv = ramp_soln_(u,v,h) to iterate again    
            sim.ramp_soln.u = iu[:,:,np.newaxis]
            sim.ramp_soln.v = iv[:,:,np.newaxis]
            sim.ramp_soln.h = ih[:,:,np.newaxis]
            
        return qTnn    
            

def basepoint_h(sim):
    # run the loop until the stopping creiateria is satisfied.
    # balance the system: (u,v)=(u(T),v(T)) and keep h 
    
    for ii in range(sim.Nz): 
    	# copy() is necessary here...
        h    = sim.soln.h[:,:,ii].copy()
        hTnn = sim.ramp_soln.h[:,:,ii].copy()
        
        # make the sequence
        # hTn = hT(n), hTnn = hT(n+1)
        # decide stopping criteria
        sim.norm_h  = absrel_norm(sim, h, hTnn)
        sim.norm_hT = absrel_norm(sim, hTnn, sim.hTn) 
        
        if sim.norm_hT <= sim.kappa:
            # convergence criteria is satisfied
            # stop optbal iterations
            sim.optbal_iterate = False
            
            if sim.time < sim.end_time :   
                # first cycle 
                
                # This balanced wave is the initial data for next test.
                # Take it before projecting balwave on physical space
                # keep initial h (soln.h) since h is base-point. 
                if sim.balanced_wave:
                	print ("Saving the balanced wave for upcoming tests...")
                	sim.balu = sim.ramp_soln.u[:,:,ii]
                	sim.balv = sim.ramp_soln.v[:,:,ii]
                	sim.balh = sim.soln.h[:,:,ii]
                
                # take the balance wave on physical space
                project_on_physicalspace(sim)

            elif sim.end_time == sim.time:  
                # second cycle
                # compute diagnosed imbalance, then attain ramp_soln to soln
                compute_norm(sim)
        
        elif sim.norm_hT > sim.kappa:
            # convergence criteria is not satisfied
            # sim.optbal_iterate = True by default
        
            # now take (u(T),v(T),h(T)) = (u(T),v(T),h) in artificial time evolution
            sim.ramp_soln.h[:,:,ii] = sim.soln.h[:,:,ii]
            
        return hTnn


def set_nonlinear_bc(sim):

    if sim.base_point == "h":
        
        print ("")
        print ("setting nonlinear_bc for base-point h...")
 
        # produce a sequence of hT = (hT0,...,hTn,hTnn,...)
        hTnn = basepoint_h(sim)
        
        print ("")
        if sim.norm_type == 'absolute':
            print (' ||h-h(T)||  =', sim.norm_h)
            print ('||hTnn-hTn|| =', sim.norm_hT)
        elif sim.norm_type == 'relative':
            print ('   ||h-h(T)||/||h||   =', sim.norm_h)
            print ('||hTnn-hTn||/||hTnn|| =', sim.norm_hT)
        
        # update hTn
        sim.hTn[:,:] = hTnn[:,:]
        
    elif sim.base_point == "q":
        
        print ("")
        print ("setting nonlinear_bc for base-point q...")
       
        # produce a sequence of qT = (qT0,...,qTn,qTnn,...)
        qTnn = basepoint_pv(sim) 
        
        print ("")
        if sim.norm_type == 'absolute':
            print (' ||q-q(T)||  =', sim.norm_q)
            print ('||qTnn-qTn|| =', sim.norm_qT        )
        elif sim.norm_type == 'relative':
            print ('   ||q-q(T)||/||q||   =', sim.norm_q)
            print ('||qTnn-qTn||/||qTnn|| =', sim.norm_qT)
        
        # update qTn
        sim.qTn[:,:] = qTnn[:,:]

