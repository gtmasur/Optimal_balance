# GTM: Ramp_conditions
#     This module provides boundary condition
#     for ramping procedure

#from Geostrophy import set_geostrophy
#from L1balancedinit import set_L1u
#from Randominit import randfield

from .Linear_BC import set_linear_bc
from .Linear_BC import build_proj_mats
from .Linear_BC import wave_separation

from .Nonlinear_BC import set_nonlinear_bc
from .Nonlinear_BC import classic_to_prognostic
from .Nonlinear_BC import prog_classic_test
from .Nonlinear_BC import compute_norm
from .Nonlinear_BC import project_on_physicalspace

# add optimal balance algorithm as a black box
from .Optimal_balance import set_kappa # TODO: newly added...
from .Optimal_balance import initialize_ramp
from .Optimal_balance import balance_wave
from .Optimal_balance import reach_linearend
from .Optimal_balance import balance_wave_L1orGeo
from .Optimal_balance import imbalances_nonlinear
from .Optimal_balance import imbalances_linear


#from Initializations import L1_initialdata
from .Initializations import initialize_fields 
from .Initializations import set_conditions
