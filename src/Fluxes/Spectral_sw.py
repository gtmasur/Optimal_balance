from .Spectral_xy import spectral_operators
from pylab import *

# GTM: general comment:
# assign name A if not doing any change in A
# copy to name A if some changes in A possible. 
# a[:,:]=b[:,:] means copying b values inside a, not assigning.

def spectral_sw_flux(sim):

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = sim.soln.h[:,:,ii]
        u = sim.soln.u[:,:,ii]
        v = sim.soln.v[:,:,ii]
        
        # Compute derivatives
        Dx_u, Dx_v, Dx_h = sim.Dx(u), sim.Dx(v), sim.Dx(h)
        Dy_u, Dy_v, Dy_h = sim.Dy(u), sim.Dy(v), sim.Dy(h)

        # Add Coriolis and derivatives
        sim.curr_flux.u[:,:,ii] = - sim.eps*(u*Dx_u + v*Dy_u) + v - Dx_h
        sim.curr_flux.v[:,:,ii] = - sim.eps*(u*Dx_v + v*Dy_v) - u - Dy_h
        sim.curr_flux.h[:,:,ii] = - (u*Dx_h+v*Dy_h) - (sim.h0+h)*(Dx_u +Dy_v)
        
        # the viscosity term is on the right hand side
        if sim.viscosity: 
            sim.curr_flux.u[:,:,ii] += (sim.eps/sim.Re)* sim.Lap(u)
            sim.curr_flux.v[:,:,ii] += (sim.eps/sim.Re)* sim.Lap(v)
        
        # divide by sim.eps
        sim.curr_flux.u[:,:,ii] /= sim.eps
        sim.curr_flux.v[:,:,ii] /= sim.eps
        
        
def apply_filter(sim):
        
    for ii in range(sim.Nz):
        ue = sim.soln.u[:,:,ii].copy()
        ve = sim.soln.v[:,:,ii].copy()
        he = sim.soln.h[:,:,ii].copy()
        
        # Project on physical space
        sim.soln.u[:,:,ii] = sim.Filt(ue)
        sim.soln.v[:,:,ii] = sim.Filt(ve)
        sim.soln.h[:,:,ii] = sim.Filt(he)
        
def spectral_sw(sim):

    sim.spectral_operators = spectral_operators

    sim.spectral_sw_flux = spectral_sw_flux
    sim.apply_filter = apply_filter
    
