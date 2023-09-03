"""
GTM: Brief explanation to data saving:
     i) For real system, the initial and end data is saved in save_additional_data,
        while other datas are saved in save_data. Then naming of datas:
        initial, 0, 1, ....., (the last output_counter), end
     ii) For artificial system, we only save the initial and end data
        for each forward and backward ramping.
        If needed then save other datas too.
     iii) diag_imb gives output for only one single epsilon value
"""

from pylab import * 
import os, sys, shutil
from scipy.fftpack import fftfreq, fft2, ifft2
                 
         
def initialize_mainpath(sim):
    # creating a folder tree of the form:
    # Outputs -- Scaling limit
    #         -- Run_name 
    #         -- Normal_run or Optimal_balance
    # make main folder 
    # add algorithm names in 
    
    # shorten scaling names
    if sim.scalinglimit=='quasi-geostrophic': slim='quasi-geo'
    elif sim.scalinglimit=='semi-geostrophic': slim='semi-geo'  
    
    if sim.optimal_balance:
        sim.path = 'Outputs/{0:s}/{1:s}_{2:s}/Optimal_balance/{3:s}_{4:s}_{5:s}'\
                    .format(slim,sim.run_name, sim.tstep,
                    sim.linear_bc,sim.base_point,sim.rampfun)
    else: 
        sim.path = 'Outputs/{0:s}/{1:s}_{2:s}/Normal_run'\
                    .format(slim,sim.run_name, sim.tstep)
    if not(os.path.isdir(sim.path)):
        os.makedirs(sim.path) 
    
    # shorten algorithm name to use in the coming name
    if sim.optimal_balance:
    
        global alg, kappan, balwave, ntype
        
        if sim.norm_type == 'relative': ntype = 'rel'
        elif sim.norm_type == 'absolute': ntype = 'abs'
        kappan = '_k{0:g}'.format(sim.n) 
        if sim.balanced_wave:         balwave = '_balwave'
        elif (not sim.balanced_wave): balwave = ''
        
        if sim.optbal_algorithm   == 'algorithm1':   alg = 'alg1'
        elif sim.optbal_algorithm == 'algorithm2':   alg = 'alg2'
        elif sim.optbal_algorithm == 'algorithm3':   alg = 'alg3'
        elif sim.optbal_algorithm == 'algorithmL1':  alg = 'algL1'
        elif sim.optbal_algorithm == 'algorithmGeo': alg = 'algGeo'
        

def initialize_imbeps_dir(sim):
    # try: no error if n is not defined!
    # n is not defined if we are not checking kappa sensitivity
    initialize_mainpath(sim)
    sim.imb_path = sim.path +'/Imbalances/Tramp{0:g}_tfor{1:g}_{2:s}'\
                    .format(sim.T,sim.tprime,alg)
    if not(os.path.isdir(sim.imb_path)):
        os.makedirs(sim.imb_path)    
        
    sim.imb_dir = sim.imb_path+'/imb_{0:s}{1:s}{2:s}{3:s}.npz'\
                .format(sim.initialsetup,balwave,kappan,ntype)
                
def initialize_imbepsinit_dir(sim):
    # I(eps) for different initial data
    initialize_mainpath(sim)
    sim.imb_path = sim.path +'/Imbalances/Tramp{0:g}_tfor{1:g}_{2:s}'\
                    .format(sim.T,sim.tprime,alg)
    if not(os.path.isdir(sim.imb_path)):
        os.makedirs(sim.imb_path)    
        
    # kmax is added for initial data test cases    
    sim.imb_dir = sim.imb_path+'/imb_{0:s}{1:s}_kmax{2:g}{3:s}{4:s}.npz'\
                .format(sim.initialsetup,balwave,sim.ak0,kappan,ntype)
                      
def initialize_imbtramp_dir(sim):
    # file of diagnosed imbalances as a function of Tramp
    initialize_mainpath(sim)
    sim.imb_path = sim.path +'/Imbalances/Tramp_tfor{0:g}_{1:s}'\
                    .format(sim.tprime, alg)
    if not(os.path.isdir(sim.imb_path)):
        os.makedirs(sim.imb_path)
    sim.imb_dir = sim.imb_path+'/imb_eps{0:g}_{1:s}{2:s}{3:s}{4:s}.npz'\
                  .format(sim.eps,sim.initialsetup,balwave,kappan,ntype)         
 
# TODO : This name is too long... Test it... 
def initialize_imbtrampinit_dir(sim):
    ## I(T) for different initial data
    initialize_mainpath(sim)
    sim.imb_path = sim.path +'/Imbalances/Tramp_tfor{0:g}_{1:s}'\
                    .format(sim.tprime,alg)
    if not(os.path.isdir(sim.imb_path)):
        os.makedirs(sim.imb_path)
    sim.imb_dir = sim.imb_path+'/imb_eps{0:g}_{1:s}{2:s}_kmax{3:g}{3:s}{4:s}.npz'\
                .format(sim.eps,sim.initialsetup,balwave,sim.ak0,kappan,ntype)          
         
         
def initialize_imbtfor_dir(sim):
    # file of diagnosed imbalances as a function of tprime
    initialize_mainpath(sim)
    sim.imb_path = sim.path +'/Imbalances/Tramp{0:g}_tfor_{1:s}'\
                    .format(sim.T,alg)
    if not(os.path.isdir(sim.imb_path)):
        os.makedirs(sim.imb_path)
    sim.imb_dir = sim.imb_path+'/imb_eps{0:g}_{1:s}{2:s}{3:s}{4:s}.npz'\
                  .format(sim.eps,sim.initialsetup,balwave,kappan,ntype)     
                    
def initialize_saving(sim):
    # creating a folder tree of the form:
    # Outputs -- Scaling limit 
    #         -- Run_name 
    #         -- (1) Normal_run           or (2) Optimal_balance
    #         -- (1) tfor(value)          or (2) linear_bc & base_point
    #         -- (1) eps(value)_initsetup or (2) FixedData or Imbalances  
    #         -- (1) data is written      or (2) T_ramp(value)_tfor(value)_alg  
    #         -- and so on...
    initialize_mainpath(sim)
    
    # TODO: this is set for proper kappa naming....
    if sim.optimal_balance:
        kastr = str(sim.kappa)
        kan = len(kastr)-1
        # for k3,k3h and k4
        if kastr[-1] == '5':
            kaval = "k{0:g}h".format(kan-2) 
        elif kastr[-1] == '1':
            kaval = "k{0:g}".format(kan-1) 
        else:
            # for eps=0.5
            kaval = "k{0:g}".format(sim.kappa)
    else: 
        kaval = '' 
    
    if sim.optimal_balance:
        # add kappa in the name
        
        # TODO: naming has changed....
        sim.path = sim.path+'/FixedData/Tramp{0:g}_tfor{1:g}_{2:s}/eps{3:g}_{4:s}{5:s}'\
                     .format(sim.T,sim.tprime,alg,sim.eps,sim.initialsetup,kaval)
    else:
        sim.path = sim.path + '/tfor{0:g}/eps{1:g}_{2:s}'\
                    .format(sim.tprime,sim.eps,sim.initialsetup)
        
    if os.path.isdir(sim.path):
       print('Output directory {0:s} already exists. '.format(sim.path) + \
             'Warning, deleting everything in the directory.') 
       # TODO: deal with deleting the directory later on
       shutil.rmtree(sim.path)    # delete an entire directory tree
    os.makedirs(sim.path)
    
    if sim.output:
        os.makedirs(sim.path+'/Realt_Data/Data/')
        if sim.optimal_balance and sim.ramp_output:
            os.makedirs(sim.path+ '/Artit_Data/Forward/Data/')
            os.makedirs(sim.path+ '/Artit_Data/Backward/Data/')
            if sim.ramp_spec_output:
                os.makedirs(sim.path+ '/Artit_Data/Forward/Specs/')
                os.makedirs(sim.path+ '/Artit_Data/Backward/Specs/')
    
    if sim.animate == 'Save':
        os.makedirs(sim.path+'/Realt_Data/Frames')
        # turn on to save ramp_frames
        if sim.optimal_balance:
            # forward and backward ramping scheme
            os.makedirs(sim.path+ '/Artit_Data/Forward/Frames')
            os.makedirs(sim.path+ '/Artit_Data/Backward/Frames')
            
    # Initialize saving stuff
    save_info(sim)
    if sim.output:
        # save data differently in both case
        if sim.optimal_balance:
            save_additional_state(sim)   # to save the initial data 
        else:
            save_state(sim)
            #sim.save_specs() # no need for initial data            
        # do not save grid at every instance..    
        #save_grid(sim)
        
        
def save_additional_state(sim):
    # to save initial data before balancing and end data after rebalancing.
    # do not increase sim.output_counter
    if sim.output_counter == 0.:
        name = 'initial'   # save initial data
    else:
        name = 'end'       # save data after rebalancing
    if sim.prognostic:
        savez_compressed(sim.path + '/Realt_Data/Data/{0:s}'.format(name),
                            u = sim.soln.u, v = sim.soln.v, h = sim.soln.h,
                            q = sim.solnp.q, delta = sim.solnp.delta, gamma = sim.solnp.gamma,
                            t = sim.time)
    else:
        savez_compressed(sim.path + '/Realt_Data/Data/{0:s}'.format(name),
                            u = sim.soln.u, v = sim.soln.v, h = sim.soln.h,
                            t = sim.time)
    
def save_state(sim):
    # Save current state
    if sim.prognostic:
        savez_compressed(sim.path + '/Realt_Data/Data/{0:04d}'.format(sim.output_counter),
                            u = sim.soln.u, v = sim.soln.v, h = sim.soln.h,
                            q = sim.solnp.q, delta = sim.solnp.delta, gamma = sim.solnp.gamma,
                            t = sim.time)
    else:
        savez_compressed(sim.path + '/Realt_Data/Data/{0:04d}'.format(sim.output_counter),
                            u = sim.soln.u, v = sim.soln.v, h = sim.soln.h,
                            t = sim.time)
    sim.output_counter += 1

def save_specs(sim):
    h = sim.soln.h[:,:,0]
    Spec = sum(abs(fft2(h))**2)
    if sim.time < sim.end_time:
        sim.Specs += [Spec]
    elif sim.time == sim.end_time:
        sim.Specs += [Spec]
        savez_compressed(sim.path + '/Realt_Data/Spec', Spec = array(sim.Specs))   

def ramp_save_state(sim):
    if sim.ramp_direction == 'Forward':
        ramp_dir = 'Forward' 
    elif sim.ramp_direction == 'Backward':
        ramp_dir = 'Backward' 

    def save(sim):
        if sim.prognostic:
            fname = sim.path + '/Artit_Data/'+ ramp_dir +'/Data/{0:04d}_cycle{1:d}_iter{2:d}'.format(
                    sim.ramp_output_counter,sim.nramping,sim.iter)
            savez_compressed(fname,
                                u = sim.ramp_soln.u, v = sim.ramp_soln.v, h = sim.ramp_soln.h, 
                                q = sim.ramp_solnp.q, delta = sim.ramp_solnp.delta, 
                                gamma = sim.ramp_solnp.gamma,
                                t = sim.ramp_time)
        else:
            fname = sim.path + '/Artit_Data/'+ ramp_dir +'/{0:04d}_cycle{1:d}_iter{2:d}'.format(
                    sim.ramp_output_counter,sim.nramping,sim.iter)
            savez_compressed(fname,
                                u = sim.ramp_soln.u, v = sim.ramp_soln.v, 
                                h = sim.ramp_soln.h, t = sim.ramp_time)
    
    # if needed, all states can be saved easily  
    #save(sim)
    #sim.ramp_output_counter += 1
    
    # save only first and end states  
    if sim.ramp_output_counter == 0:
        save(sim)
        sim.ramp_output_counter += 1
    elif sim.ramp_time == sim.T_ramp:
        save(sim)


def ramp_save_specs(sim):
    if sim.ramp_direction == 'Forward':
        ramp_dir = 'Forward' 
    elif sim.ramp_direction == 'Backward':
        ramp_dir = 'Backward' 
    h = sim.ramp_soln.h[:,:,0]
    Spec = sum(abs(fft2(h))**2)
    if sim.ramp_time < sim.T_ramp:
        sim.ramp_Specs += [Spec]
    elif sim.ramp_time == sim.T_ramp:
        sim.ramp_Specs += [Spec]
        # and save    
        fname = sim.path + '/Artit_Data/'+ ramp_dir +'/Specs/Spec_cycle{0:d}_iter{1:d}'.format(
                sim.nramping, sim.iter)
        savez_compressed(fname, Spec = array(sim.ramp_Specs))   
        

# Save grid 
def save_grid(sim):
    fname = sim.path + '/grid'
    savez_compressed(fname, x = sim.x, y = sim.y, X = sim.X, Y = sim.Y)


# Helper to write arrays          
def array2str(fp,arr):
    s = '[' + reduce(lambda s1,s2: s1+', '+s2, map(str,arr)) + ']'
    return s
        
        
# Save information
def save_info(sim):
    # Time being the following setting are fixed:
    # ========================
    # stepper      = AB3
    # method       = spectral 
    # dynamics     = nonlinear 
    # Coriolis     = f-plane
    # ========================
    fname = sim.path + '/info.txt'
    fp = open(fname, 'w')
    fp.write('Lx           = {0:g}\n'.format(sim.Lx))
    fp.write('Ly           = {0:g}\n'.format(sim.Ly))
    fp.write('Nx           = {0:d}\n'.format(sim.Nx))
    fp.write('Ny           = {0:d}\n'.format(sim.Ny))
    fp.write('Nz           = {0:d}\n'.format(sim.Nz))
    fp.write('savet        = {0:g}\n'.format(sim.savet))
    if sim.animate != 'None':
        fp.write('plott        = {0:g}\n'.format(sim.plott)) 
    fp.write('eps          = {0:g}\n'.format(sim.eps))
    fp.write('initialsetup = {0:s}\n'.format(sim.initialsetup))
    if sim.viscosity:
        fp.write('nu           = {0:g}\n'.format(sim.nu))    
#    if sim.optimal_balance:
#        fp.write('linear_bc    = {0:s}\n'.format(sim.linear_bc))
#        fp.write('base_point   = {0:s}\n'.format(sim.base_point))
#        fp.write('tprime       = {0:g}\n'.format(sim.tprime)) 
#        fp.write('T            = {0:g}\n'.format(sim.T)) 
#        fp.write('kappa        = {0:g}\n'.format(sim.kappa))  
    fp.close()
    return 

def save_imbalances(sim):
    # save diagnosed imbalance
	#savez_compressed(sim.path +'/diag_imb', imb = sim.imb)
	savez(sim.path +'/diag_imb', 
             tprime = sim.tprime,
             T = sim.T,
             epsilon = sim.eps,
             imb = sim.imb) 
             
	# write in additional info in file
	fname = sim.path + '/info.txt'
	fp = open(fname, 'a+')
	# TODO: this part has changed...!!!!
	# write kappa value after it is updated in set_kappa
	if sim.optimal_balance:
	    fp.write('linear_bc    = {0:s}\n'.format(sim.linear_bc))
	    fp.write('base_point   = {0:s}\n'.format(sim.base_point))
	    fp.write('tprime       = {0:g}\n'.format(sim.tprime)) 
	    fp.write('T            = {0:g}\n'.format(sim.T)) 
	    fp.write('kappa        = {0:g}\n'.format(sim.kappa))  
	fp.write("\n")
	if sim.base_point == "q":   
	    fp.write('norm_qT      = {0:g}\n'.format(sim.norm_qT)) 
	if sim.base_point == "h":   
	    fp.write('norm_hT      = {0:g}\n'.format(sim.norm_hT))  
	fp.write('norm_q       = {0:g}\n'.format(sim.norm_q)) 
	fp.write('norm_delta   = {0:g}\n'.format(sim.norm_delta)) 
	fp.write('norm_gamma   = {0:g}\n'.format(sim.norm_gamma))
	fp.write('norm_u       = {0:g}\n'.format(sim.norm_u))
	fp.write('norm_v       = {0:g}\n'.format(sim.norm_v))
	fp.write('norm_h       = {0:g}\n'.format(sim.norm_h))
	fp.write('iterations   = {0:s}\n'.format(str(sim.iters)))  
	fp.close()
	return       
    
