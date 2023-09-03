import numpy as np

def Euler_real(sim):

    # Compute the flux
    sim.flux()
    
    # Evolve the system
    sim.soln.u += sim.curr_flux.u*sim.dt
    sim.soln.v += sim.curr_flux.v*sim.dt
    sim.soln.h += sim.curr_flux.h*sim.dt
    
    # Store the appropriate history
    if sim.nfluxes > 0:
        sim.fluxes.u = [sim.curr_flux.u.copy()]
        sim.fluxes.v = [sim.curr_flux.v.copy()]
        sim.fluxes.h = [sim.curr_flux.h.copy()]
        sim.dts    = [sim.dt]      
        
def Euler_ramp(sim):

    # Compute the flux
    sim.ramp_flux(sim)
    
    sim.ramp_soln.u += sim.ramp_curr_flux.u*sim.ramp_dt
    sim.ramp_soln.v += sim.ramp_curr_flux.v*sim.ramp_dt
    sim.ramp_soln.h += sim.ramp_curr_flux.h*sim.ramp_dt
    
    # Store the appropriate history
    if sim.ramp_nfluxes > 0:
        sim.ramp_fluxes.u = [sim.ramp_curr_flux.u.copy()]
        sim.ramp_fluxes.v = [sim.ramp_curr_flux.v.copy()]
        sim.ramp_fluxes.h = [sim.ramp_curr_flux.h.copy()]
        sim.ramp_dts    = [sim.ramp_dt]
        
def Euler(sim):

    if sim.system == 'Real':
        return Euler_real(sim)
    elif sim.system == 'Ramped':
        return Euler_ramp(sim)
    
