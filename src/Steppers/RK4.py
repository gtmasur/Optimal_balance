import numpy as np
import sys

# y^{n+1} = y^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
# k1 = flux(y^n)
# k2 = flux(y^n + (dt/2)*k1)
# k3 = flux(y^n + (dt/2)*k2)
# k4 = flux(y^n + dt*k3)

def RK4_real(sim):

    def compute_k(sim):
        sim.flux()   # gives curr_flux
        k_u = sim.curr_flux.u
        k_v = sim.curr_flux.v
        k_h = sim.curr_flux.h
        return (k_u,k_v,k_h)

    dt = sim.dt

    # Compute k1
    k1_u, k1_v, k1_h = compute_k(sim)
    
    # Compute k2
    sim.soln.u += k1_u*dt/2.
    sim.soln.v += k1_v*dt/2.
    sim.soln.h += k1_h*dt/2.
    k2_u, k2_v, k2_h = compute_k(sim)

    # Compute k3
    sim.soln.u += -k1_u*dt/2. + k2_u*dt/2.
    sim.soln.v += -k1_v*dt/2. + k2_v*dt/2.
    sim.soln.h += -k1_h*dt/2. + k2_h*dt/2.
    k3_u, k3_v, k3_h = compute_k(sim)

    # Compute k4
    sim.soln.u += -k2_u*dt/2. + k3_u*dt
    sim.soln.v += -k2_v*dt/2. + k3_v*dt
    sim.soln.h += -k2_h*dt/2. + k3_h*dt
    k4_u, k4_v, k4_h = compute_k(sim)

    sim.soln.u += -k3_u*dt + (dt/6.)*(k1_u + 2.*k2_u + 2.*k3_u + k4_u)
    sim.soln.v += -k3_v*dt + (dt/6.)*(k1_v + 2.*k2_v + 2.*k3_v + k4_v)
    sim.soln.h += -k3_h*dt + (dt/6.)*(k1_h + 2.*k2_h + 2.*k3_h + k4_h)
    
def RK4_ramp(sim):

    def compute_k(sim):
        sim.ramp_flux()   # gives curr_flux
        k_u = sim.ramp_curr_flux.u
        k_v = sim.ramp_curr_flux.v
        k_h = sim.ramp_curr_flux.h
        return (k_u,k_v,k_h)

    dt = sim.ramp_dt
    
    # Compute k1
    k1_u, k1_v, k1_h = compute_k(sim)
    
    # Compute k2
    sim.ramp_soln.u += k1_u*dt/2.
    sim.ramp_soln.v += k1_v*dt/2.
    sim.ramp_soln.h += k1_h*dt/2.
    k2_u, k2_v, k2_h = compute_k(sim)

    # Compute k3
    sim.ramp_soln.u += -k1_u*dt/2. + k2_u*dt/2.
    sim.ramp_soln.v += -k1_v*dt/2. + k2_v*dt/2.
    sim.ramp_soln.h += -k1_h*dt/2. + k2_h*dt/2.
    k3_u, k3_v, k3_h = compute_k(sim)

    # Compute k4
    sim.ramp_soln.u += -k2_u*dt/2. + k3_u*dt
    sim.ramp_soln.v += -k2_v*dt/2. + k3_v*dt
    sim.ramp_soln.h += -k2_h*dt/2. + k3_h*dt
    k4_u, k4_v, k4_h = compute_k(sim)

    sim.ramp_soln.u += -k3_u*dt + (dt/6.)*(k1_u + 2.*k2_u + 2.*k3_u + k4_u)
    sim.ramp_soln.v += -k3_v*dt + (dt/6.)*(k1_v + 2.*k2_v + 2.*k3_v + k4_v)
    sim.ramp_soln.h += -k3_h*dt + (dt/6.)*(k1_h + 2.*k2_h + 2.*k3_h + k4_h)
    
def RK4(sim):

    if sim.system == 'Real':
        return RK4_real(sim)
    elif sim.system == 'Ramped':
        return RK4_ramp(sim)
        
