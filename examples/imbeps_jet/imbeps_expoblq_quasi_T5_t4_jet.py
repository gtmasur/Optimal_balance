from pylab import *
import matplotlib.pyplot as plt
import sys, os

sys.path.append('../../src')
import Simdata
import Steppers as Step
from PyRsw import Simulation

# Nonlinear SWE is solved in a doubly periodic 
# domain by spectral method
sim = Simulation()  # Create a simulation object

#============= Simulation settings ===================

"""
GTM: Aware that describing all these parameters became 
complicated; but, no other choice to run all test runs
depending on one src code.
"""

# Initial condition and naming
sim.d = 6                 # order of spectrum decay
sim.ak0 = 6
sim.run_name = 'Randomh_ord{0:g}'.format(sim.d)  

# General settings
sim.stepper      = Step.AB3           # Time-stepping algorithm: AB3, RK4  
sim.adaptive     = True 
sim.tstep        = 'dtu'              # dtmin or dtu
sim.viscosity    = False    
sim.scalinglimit = 'quasi-geostrophic' # semi or quasi geostrophic
sim.initialsetup = 'jet'            # geost or L1 initializations and 'zero' velocity  (or here jet)
# Prognostic variables
# (q,delta,gamma) can be also computed while plotting
sim.prognostic   = True               # compute (q,delta,gamma) in each steps

# Optimal balance parameters
sim.optimal_balance = True
if sim.optimal_balance: 
    sim.save_diagimb     = True          # save only diagnosed imbalances
    sim.rampfun          = 'exp'          # exp, quartic or cos
    sim.balanced_wave    = False          # use balanced waves for upcoming tests
    sim.optbal_algorithm = 'algorithm1'  # optimal balance algorithm1, 2, 3, L1, geo
    if sim.optbal_algorithm == 'algorithm3':
        sim.rebaltime = 'T'              # T or 2T for algorithm3
    # Boundary conditions of ramping at linear and nonlinear ends
    sim.linear_bc  = 'obliqpr'           # obliqpr, orthpr or altgeost, geost
    sim.base_point = 'q'                 # q or h
    sim.norm_type  = 'relative'          # relative or absolute

#============= Domain parameters ===================

# Specify parameters regarding Rossby number
sim.Lx  = 2*pi         # Domain extent               
sim.Ly  = 2*pi         # Domain extent                   
sim.Nx  = 255          # Grid points in x
sim.Ny  = 255          # Grid points in y
sim.Nz  = 1            # Number of layers        
sim.cfl = 0.5          # CFL coefficient       
sim.Re = 1e4           # Reynolds number         

#============= Plotting parameters ===================

sim.plot_vars = ['h','u','v']
sim.clims = [[0.24,-0.24], [], []]  

# Only output is separated for real and ramp cases.
sim.animate = 'None'        # 'Save' to save video frames,
                            # 'Anim' to animate,
                            # 'None' otherwise   
sim.output = False           # True or False in real time

if sim.optimal_balance:

    sim.ramp_output = False  # in ramp time
    sim.ramp_spec_output = False
    
    if sim.save_diagimb:
        # only save diagimbs
        sim.animate = 'None'
        sim.output = False
        sim.ramp_output = False
        sim.ramp_spec_output = False
  

#================ Functions =======================

def set_parameters(sim):

    # fixed time step with cgw = 1./sim.eps
    sim.fixed_dt = pi/(2*sim.Nx)*sim.eps  # Dritschel's code
 
    # plotting parameters
    if sim.animate == 'Anim' :  sim.plott = 50*sim.fixed_dt
    elif sim.animate == 'Save': sim.plott = sim.t_forward 
    sim.savet = sim.t_forward         # output savet in real time
    if sim.optimal_balance: 
        sim.ramp_savet = sim.T_ramp   # in ramp time
        
    print ("")
    if sim.viscosity == True:
        print ("************************* eps, Re =", sim.eps, sim.Re)
    else:
        print ("************************* eps =", sim.eps)
    if sim.optimal_balance:  
        print ("")
        print ('t_prime, T:', sim.tprime, sim.T)
        print ('t_forward, T_ramp, dt:', sim.t_forward, sim.T_ramp, sim.fixed_dt)
    else: 
        print ("")
        print ('t_prime:', sim.tprime)
        print ('t_forward, dt:', sim.t_forward, sim.fixed_dt)
    print ("")
       

def eps_run(sim):
    # run the simulation for a given epsilon value
    # initialize parameters for each eps
    sim.__init__()         
    sim.initialize()
    set_parameters(sim) 
    sim.initialize_fields(sim)
    sim.run()    
    if sim.optimal_balance:
    	return (sim.imb, sim.iters) 

"""GTM: Separate diagimb funcs are needed."""    
def diagimb_Tramp(sim,TT):  
    Simdata.initialize_imbtramp_dir(sim)
    if not os.path.isfile(sim.imb_dir):
        Tramp = [T for T in TT]
        Imb, Iters = [], []
        for sim.T in Tramp:
        	(imb, iters) = eps_run(sim)
        	Imb += [imb]
        	Iters += [iters]
        np.savez(sim.imb_dir, 
                 Tramp = Tramp,     tprime = sim.tprime, 
                 epsilon = sim.eps, kappa = sim.kappa,
                 imb = Imb,         iters = Iters)     
 
def diagimb_epsilon(sim,eeps):
    Simdata.initialize_imbeps_dir(sim)
    if not os.path.isfile(sim.imb_dir):
        epsilon = [eps for eps in eeps ]
        Imb, Iters = [], []
        for sim.eps in eeps:
            (imb, iters) = eps_run(sim)
            Imb += [imb]
            Iters += [iters]
        savez(sim.imb_dir, 
                 Tramp = sim.T,     tprime = sim.tprime,
                 epsilon = epsilon, kappa = sim.kappa,
                 imb = Imb,         iters = Iters)  
 
def diagimb_tfor(sim,ttprime):
    Simdata.initialize_imbtfor_dir(sim)
    if not os.path.isfile(sim.imb_dir):
        ttprime = [tprime for tprime in ttprime]
        Imb, Iters = [], []
        for sim.tprime in ttprime:
        	(imb, iters) = eps_run(sim)
        	Imb += [imb]
        	Iters += [iters]
        np.savez(sim.imb_dir, 
                 Tramp = sim.T,     tprime = ttprime,
                 epsilon = sim.eps, kappa = sim.kappa,
                 imb = Imb,         iters = Iters)       
                              
    
#==============  run Simulation  =====================

# fixed times: converted to T_ramp and t_forward
# depending on the scaling limit
sim.T      = 5.
sim.tprime = 4.

# put different kappa for base-points
sim.n = 4
# relative norm threshold for opt balance
sim.kappa= float("1e-{0:g}".format(sim.n))   

if not sim.optimal_balance:
    print ("")
    print (" ----- Normal run ---- ")
    sim.eps = 0.1
    eps_run(sim)

elif sim.optimal_balance:

    if not sim.save_diagimb:
        print ("")
        print (" ----- FixedData run ---- ")
        eeps = [0.1]
        for sim.eps in eeps:
            eps_run(sim)

    elif sim.save_diagimb:
        print ("")
        print (" ----- Imbalances run ---- ")
        
        # diagnosed imbalance as a function of epsilon 
        m = linspace(2,11,20)
        eeps = 2**(-m/2)
        diagimb_epsilon(sim,eeps)
                     
        # diagnosed imbalance as a function of Tramp  
#        sim.eps = 0.2
#        sim.tprime  = 0.001
#        n = 80
#        TT = 10**(np.linspace(-1.3,1.2,n))
#        TT = [0.001]           
#        diagimb_Tramp(sim, TT)
        
        # diagnosed imbalance as a function of tfor
#        sim.eps = 0.2
#        sim.T   = 0.001
#        n = 45
#        ttprime = 10**(np.linspace(-6.3,3,n))
#        ttprime = [0.001]
#        diagimb_tfor(sim, ttprime)
        
        


