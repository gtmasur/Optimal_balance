from .Spectral_xy import spectral_operators
from .Rampfun import rho_val
from pylab import * 
        
        
def spectral_ramp_sw_flux(sim):

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):
    
        h = sim.ramp_soln.h[:,:,ii]
        u = sim.ramp_soln.u[:,:,ii]
        v = sim.ramp_soln.v[:,:,ii]
        
        t = sim.ramp_time
        T = sim.T_ramp
        
        if sim.ramp_direction == 'Backward':
            c = -1.
            x = (T - t)/T   #(T_ramp-ramp_time)/T_ramp
            rho = rho_val(sim, x) 
        elif sim.ramp_direction == 'Forward':
            c = 1.
            x = t/T
            rho = rho_val(sim, x)
            
        # Compute derivatives
        Dx_u, Dx_v, Dx_h = sim.Dx(u), sim.Dx(v), sim.Dx(h)
        Dy_u, Dy_v, Dy_h = sim.Dy(u), sim.Dy(v), sim.Dy(h)
        
        # Add Coriolis and derivatives
        sim.ramp_curr_flux.u[:,:,ii] = c*(- rho*sim.eps*(u*Dx_u + v*Dy_u) + v - Dx_h)
        sim.ramp_curr_flux.v[:,:,ii] = c*(- rho*sim.eps*(u*Dx_v + v*Dy_v) - u - Dy_h)
        sim.ramp_curr_flux.h[:,:,ii] = c*(- rho*(u*Dx_h+v*Dy_h) - rho*h*(Dx_u+Dy_v) - sim.h0*(Dx_u+Dy_v) )
        
        # the viscosity term is on the right hand side
        if sim.viscosity: 
            sim.ramp_curr_flux.u[:,:,ii] += c*rho*(sim.eps/sim.Re)* sim.Lap(u)
            sim.ramp_curr_flux.v[:,:,ii] += c*rho*(sim.eps/sim.Re)* sim.Lap(v)
            
        # divide by sim.eps
        sim.ramp_curr_flux.u[:,:,ii] /= sim.eps
        sim.ramp_curr_flux.v[:,:,ii] /= sim.eps
        

def apply_ramp_filter(sim):

    for ii in range(sim.Nz):
        ue = sim.ramp_soln.u[:,:,ii].copy()
        ve = sim.ramp_soln.v[:,:,ii].copy()
        he = sim.ramp_soln.h[:,:,ii].copy()
       
        # Project on physical space
        sim.ramp_soln.u[:,:,ii] = sim.Filt(ue)
        sim.ramp_soln.v[:,:,ii] = sim.Filt(ve)
        sim.ramp_soln.h[:,:,ii] = sim.Filt(he)
        
def spectral_ramp_sw(sim):
           
    sim.spectral_operators = spectral_operators

    sim.spectral_ramp_sw_flux = spectral_ramp_sw_flux
    sim.apply_ramp_filter = apply_ramp_filter


