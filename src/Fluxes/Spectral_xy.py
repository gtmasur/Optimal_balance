"""
GTM: All spectral operators to be used all over the simulation.
"""

from pylab import * 
import sys
from scipy.fftpack import fftfreq, fft2, ifft2


def spectral_operators(sim):       # Set the differentiation operators
    
    def Dx (f):
        return ifft2(sDx*fft2(f)).real    
            
    def Dy (f):
        return ifft2(sDy*fft2(f)).real
   
    def D2x (f):
        return ifft2(sDx**2*fft2(f)).real
   
    def D2y (f):
        return ifft2(sDy**2*fft2(f)).real 
        
    def Filt (f):
        return ifft2(sFilt*fft2(f)).real

    def Div (u,v):
        return ifft2(sDx*fft2(u) + sDy*fft2(v)).real
        
    def Div_perp (u,v):
        return ifft2(-sDy*fft2(u) + sDx*fft2(v)).real

    def Lap (h):
        return ifft2(sLap*fft2(h)).real

    def ILap (h):
        return ifft2(sILap*fft2(h)).real
        
    def IHelm (h):
        # sim.sIHelm is defined in Nonlinear_BC
        return ifft2(sim.sIHelm*fft2(h)).real

    
    # 2*pi/sim.Lx cancels out: but in case of dimensional run
    kx = 2*pi/sim.Lx*fftfreq(sim.Nx)*sim.Nx
    sDx = 1j*kx[:,newaxis]*ones((sim.Nx, sim.Ny))  # Spectral x-Derivative
    
    ky = 2*pi/sim.Ly*fftfreq(sim.Ny)*sim.Ny
    sDy = 1j*ky[newaxis,:]*ones((sim.Nx,sim.Ny)) # Spectral y-Derivative 

    # spectral Laplacian  
    # Do not mix - and +: spectral laplacian has - entries
    sLap = (sDx**2 + sDy**2).real
    
    # spectral filter
    Nalias = sim.Nx/3.0
    sFilt = array(Nalias**2+sLap > 0, dtype=float)

    # inverse Laplacian
    sILap = 1.0/(sLap + (sLap==0.0))
    sILap[0,0] = 0.0
    
    sim.sDx = sDx
    sim.sDy = sDy
    sim.Dx = Dx
    sim.Dy = Dy
    sim.D2y = D2y
    sim.D2x = D2x
    sim.sFilt = sFilt  # needed in Helmholtz solver
    sim.Filt = Filt
  
    sim.sLap = sLap 
    sim.Div = Div
    sim.Div_perp = Div_perp
    sim.Lap = Lap    
    sim.ILap = ILap
    sim.IHelm = IHelm
    
    
