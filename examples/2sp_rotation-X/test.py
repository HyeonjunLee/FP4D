# FP4D_reader_v3_2
## 2025.02.12
# - Designed based on FP4D_v1.1_mpi_oo_250210.
# - To use iopt_NEO, NEO must be run based on out.fp4d.neo_input.
# - For simulation, Z number and mass number are provided as input, but for plotting, the proton number must be input to determine the species.

import sys
sys.path.append('/home/youknow/python/FP4D_reader_v3_2')
from import_FP4D_reader_v3_2 import *


# Start Main: Input the sim_path for FP4D and NEO, and the proton number of the species.

# e,D,W+35,He-3
num_proton = [-1,1,74,2]

folder_path = '/home/youknow/kairos/FP4D_v1.1.a/n20_keV1_new/4sp-He3W35_rotation-O_nW-1e-5-Nt3_flow0'

iopt_NEO = True

neo_path = '/home/youknow/Codes/Works/NEO_for_FP4D/n20_keV1/4sp-He3W35_rotation-O_nW-1e-5'


# START

# READ DATA
data, data_path = load_data(folder_path)

# Set comman, global, vars
FP4D = global_FP4D(data, data_path, folder_path)
FP4D()

PLOT = PLOT_FP4D(data, data_path, num_proton)
PLOT()

FP4D.print_time()

Full, f0, f1, g1, Drift, PS, Qf, Cf, NB, RHS, expand_f0_M, expand_drift     = call_PostProcessing(data, FP4D.input['OPT_TIME_SOLVER'], FP4D.input['OPT_TIME_SOLVER'])

FP4D.init_drift()
g1_flux = FP4D.drift_psi[np.newaxis,:,:,:,:,:] * data['g1']

g1_flux = PostProcessing(data,g1_flux)
g1_flux()

if FP4D.input['OPT_ROTATION'] != 0:
    FP4D.init_drift_rot()
    g1_flux_rot = FP4D.drift_psi_rot[np.newaxis,:,:,:,:,:] * data['g1']
    g1_flux_rot = PostProcessing(data,g1_flux_rot)
    g1_flux_rot()

# friction compare w analytic [helander] -eq(8.29)
col_flux = PostProcessing(data,FP4D.col_flux)
col_flux()


# Compare the FSA_UparB
# 1. NEO way
# 2. youknow way

if (iopt_NEO):
    Full.exe_UparB_NEO()
    FSA_UparB_NEO = Full.FluxSurface_average(Full.UparB_NEO)

    # compare FSA_UparB & FSA_UparB_NEO
    print('ratio: UparB_FP4D of UparB_NEO')
    print(Full.FSA_UparB/FSA_UparB_NEO)
    # ==> 0.99918

# # EXTRACT_NEO
import sys
sys.path.append('/home/youknow/python/tool_NEO-EXTRACT')

from tool_NEO    import *


if (iopt_NEO):
    # Double-check n_norm. It must be compared with the actual NEO input.
    # To prevent issues, NEO_input was created in FP4D, and using it as is should be fine.
    # Don't forget that electrons must be processed first.
    n_norm_NEO = FP4D.n0[0,0]
    m_norm_NEO = 3.3452e-27 # mD in FP4D : 3.3440e-27
    T_norm_NEO = FP4D.T0[0,0]
    v_norm_NEO = np.sqrt(T_norm_NEO*data['eV2J']/m_norm_NEO)

    B_norm_NEO = FP4D.B0 # Need to modify ... in general case ...
    a_norm_NEO = FP4D.a0

    jparB_norm_NEO = data['eV2J'] * n_norm_NEO * v_norm_NEO * B_norm_NEO


# n_norm = 1e20; T_norm = 1000; B_norm = 1.8; a=0.5
if (iopt_NEO):
    NEO = extract_NEO(neo_path, n_norm_NEO, T_norm_NEO, B_norm_NEO, a_norm_NEO)
    NEO()

if (iopt_NEO):
    print('NEO kpar', NEO.klittle_upar)
    print('NEO jparB', NEO.jparB)
    print('NEO normed jparB',NEO.jparB/(NEO.n_norm*1.602e-19*NEO.v_norm*NEO.B_norm))


print(Full.kpar[:,0,-1])

# current
for ipsi in range(FP4D.Npsi):
    if (iopt_NEO):
        PLOT.jpar(Full, ipsi, iopt_NEO, jparB_norm_NEO, NEO.B_norm, t_index='step')   
    else:
        PLOT.jpar(Full, ipsi)


# kpar
isp_index = range(FP4D.Nsp)

for ipsi in range(FP4D.Npsi):
    PLOT.kpar(Full, ipsi, t_index='step',
             isp_index = isp_index)


# moment over time for theta

# isp_index = range(FP4D.Nsp)

# ith = 0

# for ipsi in range(FP4D.Npsi):
#     PLOT.moment(Full, ipsi, ith, t_index = 'step',
#                isp_index = isp_index)

# # PLOT.moment(f1, ipsi, ith, t_index = 'step',
# #            isp_index = isp_index)


# (FSA) moment over time

isp_index = range(FP4D.Nsp)

for ipsi in range(FP4D.Npsi):
    PLOT.FSA_moment(Full, ipsi, t_index = 'step',
                   isp_index = isp_index)

# PLOT.FSA_moment(f1, ipsi, t_index = 'step',
#                isp_index = isp_index)


# moment versus theta
isp_index = range(FP4D.Nsp)

it = -1

for ipsi in range(FP4D.Npsi):
    PLOT.moment_theta(Full, ipsi, it,
                      isp_index = isp_index)

# PLOT.moment_theta(f1, ipsi, it, t_index = 'step',
#                   isp_index = isp_index)


# (FSA) temperature over time
isp_index = range(FP4D.Nsp)

for ipsi in range(FP4D.Npsi):
    PLOT.FSA_temperature(Full, ipsi, t_index = 'step',
                         isp_index = isp_index)

# Flux Check
ipsi=0
it=-1

for ipsi in range(FP4D.Npsi):
    print('----------------------------')
    print('* FP4D flux (dr)')        
    print('----------------------------')
    FP4D.print_flux_FP4D(g1_flux, FP4D.inv_dpsidr, ipsi, it)
    print("")

    print('----------------------------')
    print('* collisional FP4D flux (dr)')        
    print('----------------------------')
    FP4D.print_flux_FP4D(col_flux, FP4D.inv_dpsidr, ipsi, it)
    print("")

    print('----------------------------')
    print('* NEO flux (dr)')               
    print('----------------------------')
    FP4D.print_flux_NEO(NEO, 1.0, ipsi)

# flux over time for theta

# isp_index = range(FP4D.Nsp)

# ith = 1

# for ipsi in range(FP4D.Npsi):
#     PLOT.flux(g1_flux, ipsi, ith, t_index = 'step',
#              isp_index = isp_index)


# (FSA) flux over time

isp_index = range(FP4D.Nsp)

for ipsi in range(FP4D.Npsi):
    PLOT.FSA_flux(g1_flux, ipsi, t_index = 'step',
                 isp_index = isp_index)


isp_index = [1,2]

for ipsi in range(FP4D.Npsi):
    PLOT.FSA_flux(g1_flux, ipsi, t_index = 'step',
                 isp_index = isp_index)

# flux versus theta

# isp_index = range(FP4D.Nsp)

# it = -1

# for ipsi in range(FP4D.Npsi):
#     PLOT.flux_theta(g1_flux, ipsi, it, 
#                     isp_index = isp_index)


# volume integral terms versus theta

ipsi=0
isp=0
it1=0; it2=-1

# for isp in range(FP4D.Nsp):
#     print("isp",isp)
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, Drift.dens, 'Drift.dens')
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, PS.   dens, 'PS.   dens')
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, RHS.  dens, 'RHS.  dens')
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, f1.   dens, 'f1.   dens')
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, g1.   dens, 'g1.   dens')
#     print(" ")
    
# for isp in range(FP4D.Nsp):
#     print("isp",isp)
#     PLOT.plt_theta_int_vol(isp, ipsi, it1, it2, Qf.   eni, 'Qf.eni')
#     print(" ")


# collisional power deposition for theta

jsp_index = range(FP4D.Nsp)

ith = 0

for ipsi in range(FP4D.Npsi):
    for isp in range(FP4D.Nsp):
        PLOT.plt_col_power(FP4D.col_power, ipsi, ith, isp, t_index = 'step',
                           jsp_index = jsp_index)


# (FSA) collisional power deposition
jsp_index = range(FP4D.Nsp)

for ipsi in range(FP4D.Npsi):
    for isp in range(FP4D.Nsp):
        PLOT.plt_FSA_col_power(FP4D.FSA_col_power, ipsi, isp, t_index = 'step',
                               jsp_index = jsp_index)


# Potential
ipsi = 0
it = -1
    
PLOT.Phi_1(ipsi=ipsi,it=it,theta_mode = '0_to_2pi')
PLOT.Phi_1(ipsi=ipsi,it=it,theta_mode = '-pi_to_pi')


# Comparison of the potential FP4D with NEO

# PHI : FP4D vs NEO
ipsi = 0; it = -1

if (iopt_NEO):
    for ipsi in range(FP4D.Npsi):
        PLOT.Phi_1_NEO(ipsi, it, NEO.GRID_THETA, NEO.PHI)


if FP4D.input['OPT_ROTATION'] != 0 and (iopt_NEO):
    for ipsi in range(FP4D.Npsi):
        PLOT.Phi_0_NEO(ipsi, it, NEO.GRID_THETA, NEO.phi_rot)


if FP4D.input['OPT_ROTATION'] != 0 and (iopt_NEO):
    for ipsi in range(FP4D.Npsi):
        PLOT.Phi_tot(ipsi, it)


# 0th order density w/ rotation

it=0

if FP4D.input['OPT_ROTATION'] != 0 and (iopt_NEO):
    for ipsi in range(FP4D.Npsi):
        for isp in range(FP4D.Nsp):
            PLOT.n0_NEO(ipsi, isp, it, NEO, f0)

# n0+n1=>total density using f1
it=-1

if FP4D.input['OPT_ROTATION'] != 0 and (iopt_NEO):
    for ipsi in range(FP4D.Npsi):
        for isp in range(FP4D.Nsp):             
            PLOT.n0_n1_tot(ipsi, isp, it, f0, f1)


# n0+n1=>total density using g1

# it=-1

# if FP4D.input['OPT_ROTATION'] != 0 and (iopt_NEO):
#     for ipsi in range(FP4D.Npsi):
#         for isp in range(FP4D.Nsp):             
#             PLOT.n0_n1_tot(ipsi, isp, it, f0, g1)


# Coutour plot of distribution function

# expand_drift = np.expand_dims(data['dfdt_d'], axis=0)
# expand_drift = np.tile(expand_drift, (FP4D.Nt,1,1,1,1,1))

RHS_sum = data['dfdt_Q'] + data['dfdt_col']/data['step_Col'] + data['dfdt_d'][np.newaxis,:,:,:,:,:]


plt_list = [
    Full.dist,
    f0.dist,
    f1.dist,
    g1.dist,
    g1.dist-f1.dist,
    data['dfdt_Q'],
    data['dfdt_col']/data['step_Col'],
    data['dfdt_par'],
    expand_drift,
    RHS_sum]

titles_list= [
    r'$f_{0a}+f_{1a}$',
    r'$f_{0a}$',
    r'$f_{1a}$', 
    r'$g_{1a}$',
    r'$g_{1a}-f_{1a}=f_{0a} \frac{Z_a e \Phi_1}{T_{0a}}$',
    r'$Q_a(f_{0a}) \times \Delta t$',
    r'$(C+Q_a^{implicit})(f_{0a}^*+f_{1a}^*) \times \Delta t$',
    r'$-v_\parallel \hat{\mathbf{b}}\cdot \mathbf{\nabla|_\mathit{\mu,W}}\,g_a \times \Delta t$',
    r'$-\mathbf{v_{ma}}\cdot \mathbf{\nabla|_\mathit{\mu,W}}\,f_{Ma} \times \Delta t$',
    r'$f_{1a}^{n} - f_{1a}^{n-1} = RHS \times \Delta t$'] 


# True : contourf, Fals : contour
iopt_contourf = True
# iopt_contourf = False

n_levels = 100

it=0
isp=0
ith=8
ipsi=0

PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)

it=-1
isp=0
ith=FP4D.Nth//4
ipsi=0

PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)
    


# # True : contourf, Fals : contour
# iopt_contourf = True
# # iopt_contourf = False


# n_levels = 100


# it=0
# isp=0
# ith=8
# ipsi=0

# PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)

# it=-1
# isp=0
# ith=FP4D.Nth//4
# ipsi=0

# PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)
    


# # True : contourf, Fals : contour
# iopt_contourf = True
# iopt_contourf = False


# n_levels = 100


# it=-1
# isp=1
# ith=FP4D.Nth//4
# ipsi=0

# PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)


# # True : contourf, Fals : contour
# iopt_contourf = True
# iopt_contourf = False


# n_levels = 100


# it=-1
# isp=2
# ith=FP4D.Nth//4
# ipsi=0

# PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)


# # True : contourf, Fals : contour
# iopt_contourf = True
# iopt_contourf = False

# n_levels = 100


# it=-1
# isp=3
# ith=FP4D.Nth//4
# ipsi=0

# PLOT.contour(iopt_contourf, plt_list,titles_list, isp,ipsi,ith,it, n_levels)



plt_dist_list = [f0.dist,
                 f1.dist,
                 f0.dist+f1.dist,
                 g1.dist
                ]
plt_name_list = ['f0',
                 'f1',
                'f0+f1',
                'g1']


it=-1;isp=0;ith=0;ipsi=0;

ir_list=range(0, FP4D.Nvperp//8)

PLOT.plt_1d_vpar(plt_dist_list, plt_name_list, 
                 isp, ipsi, ith, it, ir_list)


it=-1;isp=0;ith=0;ipsi=0;

iz_list=range(FP4D.Nvpar//8, FP4D.Nvpar//8*2)

PLOT.plt_1d_vperp(plt_dist_list, plt_name_list, 
                 isp, ipsi, ith, it, iz_list)


plt_dist_list_2 = [f0.dist,
                 f1.dist,
                ]
plt_name_list_2 = ['f0',
                 'f1'
                ]


it=-1;isp=0;ith=FP4D.Nth//4;ipsi=0;

ir_list=range(0, FP4D.Nvperp)

PLOT.plt_1d_vpar(plt_dist_list_2, plt_name_list_2, 
                 isp, ipsi, ith, it, ir_list)

