#! /usr/bin/env python3

from pylab import *
from numpy.fft import *

# This random field is written for 2pi doubly periodic domain.
# Need to be adjusted for L periodic domain

# Use TeX to typeset figure annotations:
rc('text', usetex=True)
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

N = 255  # grid points
d = 6    # spectrum decay
ak0 = 6  # maximum wave number
# Number of grid points in each direction

Nalias = N/3
# This is an example of the 2/3 rule.  In Sergey's example code, there
# is a more sophisticated version on a half-staggered grid which
# requires less dealiasing.  Note: In this file, we don't actually
# require dealiasing - I just leave it in to show the code pattern.

rks = N*fftfreq(N)        # wave numbers derived from fftfreq
Lap = (rks**2)[:,newaxis] + (rks**2)[newaxis,:] # Spectral Laplacian
Filt = array(Nalias**2-Lap>=0, dtype=float)     # Aliasing Filter

kk = arange(N, dtype=float)

def k (ix,iy):
    return int(sqrt(Lap[ix,iy]))

def spectrum(u):
    s = zeros(N)
    uft = fft2(u)
    for i in range(N):
        for j in range(N): 
            s[k(i,j)] += abs(uft[i,j])**2 
    return s*(2*pi/N)**2        

def rand_field (d, ak0, s=4258463, ham=0.2):
    """Produce random field which an asymptotic spectral decay of
       k^-(d) and a spectral maximum at ak0"""
    
    # First, produce random field with absolutely flat spectrum
    seed(s) # Seed random number generator
    ur = rand(N,N)
    sr = spectrum(ur)
    
    uft = fft2(ur)
    for i in range(N):
        for j in range(N):
            uft[i,j] /= (sqrt(sr[k(i,j)]) + 1e-12) 
    
    n = (7+d)/4.0     # Controls the scaling in the decaying part of the spectrum
    a = (2*n-3)/3.0   # Ensures that spectrum peaks at k=ak0
    ss = (Lap)**(3.0/2.0)/(Lap+a*ak0**2)**n
    hh = ifft2(ss*uft*Filt).real
    hh *= ham/max( abs(amax(hh)),abs(amin(hh)))  # this gives |hh|<0.2
    return hh

uu = rand_field(d,ak0)  
lapuu = ifft2(fft2(uu)*Lap).real
uu2 = ifft2(Filt*fft2(uu**2)).real
specuu = spectrum(uu)


# The following shows how to draw a reference line which touches the
# spectrum at wavenumber kref with spectral slope duu
kref = 50
duu = - d
Cuu = specuu[kref]/kk[kref]**duu

loglog (kk, specuu, '-r',
        kk, spectrum(lapuu), '-g',
        kk, spectrum(uu2), '-b',
        kk[1:], Cuu*kk[1:]**duu, ':r')

legend ((r'$u$',
         r'$\Delta u$',
         r'$u^2$',
         r'$O(k^{' + str(duu) + '})$ reference'))
xlim(0,N)
ylim(1e-6,1e9)
#savefig('initialspectra_d{0:g}_k{1:g}.pdf'.format(d,ak0))

# Generate PDF file with the image
#savefig('spectrum_{0:d}_{1:d}.pdf'.format(d,N))

figure()

imshow(uu, vmin = -0.2, vmax = 0.2, cmap = 'seismic')
colorbar()
title('The file $u$ in physical space')
#savefig('initframe_d{0:g}_k{1:g}.pdf'.format(d,ak0))

show()

# Save data
fname = './randomh_{0:d}_{1:d}'.format(d,N)
np.savez_compressed(fname, h = uu )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# the eps range we work on
m = linspace(2,11,60)
eeps = 2**(-m/2)

# I am planning to run it for only 20 values, linspace(2,11,20),
# for time being to obtain quick results





