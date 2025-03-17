/* --------------------------------------------------------------- */  
/*    Header file for constant varialbes which are                 */
/*    frequently used in gyrokinetic tokamak simulation.           */
/*                                                                 */
/*    Note that all varialbes in this file have double physical      */
/*    dimension.                                                   */
/* --------------------------------------------------------------- */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>
#include <memory.h>
#include <assert.h>
#include <mpi.h>
#include <fftw3.h>
#ifdef FLAG_KAIROS
  #include <umfpack.h>
#else // hnuc
  #include <suitesparse/umfpack.h>
#endif

// HYPER
/*#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_ls.h"
*/
/* ------------------------------------------------- */
/* switch this definition to change double and float */
#define double_sin(x) sin((x))
#define double_cos(x) cos((x))
#define double_tan(x) tan((x))
#define double_atan(x) atan((x))
#define double_acos(x) acos((x))
#define double_cosh(x) cosh((x))
#define double_exp(x) exp((x))
#define double_log(x) log((x))
#define double_log10(x) log10((x))
#define double_fabs(x) fabs((x))
#define double_sqrt(x) sqrt((x))
#define double_erf(x) erf((x))
#define double_erfc(x) erfc((x))
#define double_pow(x, y) pow((x), (y))
#define double_tanh(x) tanh((x))
#define double_atan2(x,y) atan2((x),(y))
#define MPI_MYREAL MPI_DOUBLE
// typedef double double;
typedef double complex dcomplex;
/* ------------------------------------------------- */


/* definition for option ---------------------------- */
// uncomment for debug with more printing
//#define PR_DEBUG 1

// uncomment for GAM simulation
//#define GAM_SIM 1

// uncomment to include FLR correction
#define DEF_FLR_COR 1

// uncomment for 3D field diagnostics
//#define DEF_DIAG_3DF 1

// uncomment for fluid electron
//#define EFLUID 1

// uncomment for de-aliasing
//#define DEF_DE_AL 1

// uncomment to remove asymmetic All
//#define DEF_ZERO_All0 1

// uncomment to turn-on electron parallel Landau damping
//#define DEF_LANDAU 1

// uncomment to apply Neumann BC at core for n=0 mode
#define DEF_CBC_NEU 1

// uncomment to turn-on resistivity term in Ohm's law
//#define EF_RESIS 1

/* uncomment to turn-off fluid moment based passing electron weigth reset
         and to turn-on passing electron weight evolution along trajectory */
#define DEF_EPASS_WGT 1

// uncomment for bounced averaged drift-kinetic electrons for TEM
//#define BADKE 1
#define N_BTH 128
#define BADK_NBA_MAX 60
#define BADK_NBA_MIN 10
#define N_BADK_LB_FBP 80
#define BADK_VP2_MIN 1.0e-10
#define BADK_VL2_MIN 1.0e-10
#define BADK_VV_CUT 0.1

// uncoment for qausi-linear background evoltuion uding DG
//#define DG_FQL 1

// uncomment for ETG simulation with electrons as main species
//#define DEF_ETG_SIM 1
/* ------------------------------------------------- */

#define NUM_SP_MAX  3
/*#define NUM_ISP     2
#define SP_IDX(ip) ((ip) < g_ptl.np[0] ? SP_ion : \
                   ((ip) < g_ptl.np_ions ? SP_Imp : SP_ele))*/

#define NUM_PVAR 8

/* particle flag related definitions */
#define BPF_NULL 0x00000000 // null flag
#define init_pflag(flag) (flag = BPF_NULL)

#define BPF_LOST 0x00000001 // lost particle flag
#define BPF_TEST 0x00000002 // test particle flag
#define BPF_XTH  0x00000004 // pth-coordinate flag
#define BPF_OUT  0x00000008 // self-consistent domain out flag
#define BPF_RCYL 0x00000010 // recycling flag
#define BPF_ADW  0x00000020 // full adiabatic response flag
#define BPF_TRAP 0x00000040 // trapped particle flag

#define chk_pflag_notlost(flag) ((flag & BPF_LOST) != BPF_LOST)
#define chk_pflag_lost(flag)    ((flag & BPF_LOST) == BPF_LOST)
#define chk_pflag_test(flag)    ((flag & BPF_TEST) == BPF_TEST)
#define chk_pflag_RZ(flag)      ((flag & BPF_XTH)  != BPF_XTH)
#define chk_pflag_xth(flag)     ((flag & BPF_XTH)  == BPF_XTH)
#define chk_pflag_in(flag)      ((flag & BPF_OUT)  != BPF_OUT)
#define chk_pflag_out(flag)     ((flag & BPF_OUT)  == BPF_OUT)
#define chk_pflag_rcyl(flag)    ((flag & BPF_RCYL) == BPF_RCYL)
#define chk_pflag_adw(flag)     ((flag & BPF_ADW)  == BPF_ADW)
#define chk_pflag_trap(flag)    ((flag & BPF_TRAP) == BPF_TRAP)

#define set_pflag_lost(flag)    (flag = flag | BPF_LOST)
#define set_pflag_test(flag)    (flag = flag | BPF_TEST)
//#define set_pflag_RZ(flag)      (flag = flag ^ BPF_XTH)
#define set_pflag_xth(flag)     (flag = flag | BPF_XTH)
#define set_pflag_out(flag)     (flag = flag | BPF_OUT)
#define set_pflag_rcyl(flag)    (flag = flag | BPF_RCYL)
#define set_pflag_adw(flag)     (flag = flag | BPF_ADW)
#define set_pflag_trap(flag)    (flag = flag | BPF_TRAP)

#define unset_pflag_lost(flag)  (flag = flag ^ BPF_LOST)
#define unset_pflag_test(flag)  (flag = flag ^ BPF_TEST)
//#define unset_pflag_RZ(flag)    (flag = flag | BPF_XTH)
#define unset_pflag_xth(flag)   (flag = flag ^ BPF_XTH)
#define unset_pflag_out(flag)   (flag = flag ^ BPF_OUT)
#define unset_pflag_rcyl(flag)  (flag = flag ^ BPF_RCYL)
#define unset_pflag_adw(flag)   (flag = flag ^ BPF_ADW)
#define unset_pflag_trap(flag)  (flag = flag ^ BPF_TRAP)

/* finer grid for interpolation */
#define N_FTH 128

/* grid for phase space diagnostics */
#define N_DG_PHS_XX 4
#define N_DG_PHS_TH 1

/* definition for PADE poisson solver */
#define N_PADE_MAT     10
#define PADE_MID_L     0
#define PADE_MID_D     1
#define PADE_MID_D1    2
#define PADE_MID_R     3
#define PADE_MID_S     4
#define PADE_MID_N     5
#define PADE_MID_N1    6
#define PADE_MID_Nfx   7
#define PADE_MID_N1fx  8
#define PADE_MID_SB    9

#define N_GYP_MAX 4
#define CHARLEN 1000
#define N_NULLX 1

/* for diag-f */
#define N_DGF_XX  10
#define N_DGF_THE 12
#define N_DGF_NV  30
#define N_DGF_VTOT (N_DGF_NV*(2*N_DGF_NV-1))

inline int idx_dgf(int l, int p, int nv);


/* fori read-eprof.c */
#define DEF_RPF_NUM 3

/* for runtime.c */
#define N_DGRT 8


/* ------------------------------------------------- */
inline dcomplex phi_thb(int n, int i, int j);
inline dcomplex All_thb(int n, int i, int j);
inline dcomplex fld_thb(dcomplex f[], int n, int i, int j);
inline dcomplex thbf(int n, int i, int j);
inline dcomplex inv_thbf(int n, int i, int j);
inline int jmod_th(int j);
inline int idx_s3(int i, int j);
inline int idx_s5(int i, int j);
inline int idx_1d_CS(int i);
inline int idx_1d_QS(int i);


/* ------------------------------------------------ */
inline int idxc(int i, int j);
inline int idxV(int i, int j);

inline int idxn_ij(int n, int ij);
inline int idx_nV25(int n, int ij, int ijp);
inline int idxnij(int n, int i, int j);
inline int idx_thb(int n, int i);
inline int idxM00(int i, int ip);
inline int idx_nskin(int n, int ij, int ijp);
inline int idxMp(int n, int ij);
inline int idxMn(int n, int ij, int ijp);
inline int idx_pconv(int sp_id, int ic);
inline int idx_pdcprf(int sp_id, int ic, int i);
/* ------------------------------------------------- */

/* ------------------------------------------------- */
inline double modulus(double x);
double mod_th(double x); 
inline double mod_zt(double x);
double th_fn(double R, double Z);
/* ------------------------------------------------ */

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define absv(x) ((x) >= 0.0 ? (x) : -(x))
/* ------------------------------------------------- */

/* 0th-spline */
inline double ngp(double x);

/* 1st-spline */
double cic(double x);

/* 2nd-spline */
double tsc(double x);
inline double dtsc(double x);

/* 3rd spline */
inline double csc(double x);

/* ------------------------------------------------- */
#define psi_fn(R,Z) (intp_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
                                g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.psi_CS))
#define I_fn(R,Z) (intp_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
			g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.I_CS))
#define B_fn(R,Z) (intp_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
                              g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.B_CS))

#define intp_dpsi_RZ(R,Z,pt_val) \
	intp_dfn_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
	               g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.psi_CS, pt_val)
#define intp_dxx_RZ(R,Z,pt_val) \
	intp_dfn_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
	               g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.xx_CS, pt_val)
#define intp_dI_RZ(R,Z,pt_val) \
	intp_dfn_2d_CS(R, Z, g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR,\
	               g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.I_CS, pt_val)

#define intp_RZ_xfth(xx, th, pt_R, pt_Z) \
        intp2_2d_CS(xx, mod_th(th), g_xth.xx0, g_xth.xx1, g_xth.nxx, \
	            0.0, g_c.twopi, N_FTH,\
                    g_xth.R_CS, g_xth.Z_CS, pt_R, pt_Z)

#define epm_fn(xx) intp_1d_CS(xx, g_xth.xx0, g_xth.xx1, g_xth.nxx, g_xth.epm_CS)
/* ------------------------------------------------- */

/* ------------------------------------------------- */
inline void notnan(double x);
/* ------------------------------------------------- */


/* for QUADPACK ------------------------------------------------------- */
double dqag(double f(),double a,double b,double epsabs,
            double epsrel,int irule,double *abserr,int *neval,int *ier);
#define my_G_K(f,a,b) dqag(f,a,b, g_qdpk.epsabs, g_qdpk.epsrel, g_qdpk.irule,\
                           &(g_qdpk.r1),&(g_qdpk.neval),&(g_qdpk.iter))
/* for QUADPACK ------------------------------------------------------- */


/* for phase space diagnostics ---------------------------------------- */
inline int idx_dg_phsx(int i, int j);
inline int idx_dg_phsv(int l, int p, int k);
/* for phase space diagnostics ---------------------------------------- */


/* for fluid passing electron ----------------------------------------- */
#define EF_IMEX_MAX 4

#define N_EF_TVECTOR   2
#define EF_TVID_dAldt  0
#define EF_TVID_dnedt  1

#define N_EF_VECTOR  4
#define EF_VID_All   0
#define EF_VID_dne   1
#define EF_VID_phi   2
#define EF_VID_dje   3

#define N_EF_MATRIX   20
#define EF_MID_CJ     0
#define EF_MID_CA     1
#define EF_MID_CP     2
#define EF_MID_CD     3
#define EF_MID_AD     4
#define EF_MID_AA     5
#define EF_MID_AP     6
#define EF_MID_NN     7
#define EF_MID_LA     8
#define EF_MID_SS     9
#define EF_MID_DD     10
#define EF_MID_DDl    11
#define EF_MID_DDp    12
#define EF_MID_AJ     13
#define EF_MID_CPT    14
#define EF_MID_APT    15
#define EF_MID_RR     16
#define EF_MID_LDA    17
#define EF_MID_DDl_LD 18
#define EF_MID_SS_LD  19

#define N_EF_TENSOR  4
#define EF_TID_CAJ   0
#define EF_TID_CPD   1
#define EF_TID_AAD   2
#define EF_TID_AAP   3

#define N_EFV 4

//#define DEF_HYPRE 1

int idx_ef_vec(int n, int ij);
int idx_ef_mat(int myc_nn, int myc_ij, int dip, int djp);
int idx_ef_ten(int myc_nn, int myc_ij,
               int nnp, int dip, int djp, int dipp, int djpp);
int idx_ef_big_vec(int ij, int p, int myc_nn);
/* for fluid passing electron ----------------------------------------- */

/* End of cnst-srt1.h */
/* START of gvar-srt.h ----------------------------------------------------------------------- */
/* enum type:           -1       0       1       2       3      4       */
enum species_type { sp_non=-1, sp_ion, sp_ele, sp_imp, sp_ep, sp_eps };
typedef enum species_type sp_t;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* simulation step and time */
struct simulation_time_type
{
  int flag0;

  double time, dt;
  int niter, nstep;
  int pstep, pstep0, pstep1;

  /* skipping frequency */
  int skip_coll, skip_snsp;
  int skip_dg_phi, skip_dg_eng;
  int skip_nrdc;
};
typedef struct simulation_time_type stime_t;

extern stime_t g_stime;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* particle variables */
struct phase_variable_type
{ 
  double x[NUM_PVAR];
};
typedef struct phase_variable_type phsv_t;

struct particle_variable_type
{
  int flag;
  double x[NUM_PVAR];
  double B, xx, Iv, tcm, vv, Ik;
  double xg[N_GYP_MAX], tg[N_GYP_MAX];
};
typedef struct particle_variable_type ptlv_t;

struct particle_type
{
  int flag0;

  // number of particles for each species
  int np[NUM_SP_MAX];
  int nummk[NUM_SP_MAX];
  int np_tot, np_ion, np_imp, np_ele, np_max;

  // number of gyro-points etc.
  int ng[NUM_SP_MAX];
  double Zs_ng[NUM_SP_MAX];

  // geometric information for particle loading
  double V_load[NUM_SP_MAX], xx_load[NUM_SP_MAX];

  ptlv_t *array[NUM_SP_MAX], *array0[NUM_SP_MAX];

  double rhoR[NUM_SP_MAX][N_GYP_MAX], rhoZ[NUM_SP_MAX][N_GYP_MAX];

  phsv_t *k0[NUM_SP_MAX], *k1[NUM_SP_MAX], *k2[NUM_SP_MAX];
};
typedef struct particle_type pl_t;

extern pl_t g_ptl;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
struct xxp_theta_grid_point_type
{
  double R, Z;
};
typedef struct xxp_theta_grid_point_type grid_pt_t;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
struct integral_info_qdpk_type
{
  double xx, th1, th2;
  double r1, r2, r3;
  double dlam, lam0, Mep;
  int is, n, i, ip, j, jp;

  double epsabs, epsrel;
  int irule, neval, iter;
};
typedef struct integral_info_qdpk_type qdpk_t;
extern qdpk_t g_qdpk;

struct integral_info_cuba_type
{
  char mtype;
  double xx1, xx2, th1, th2, B, lB1, lB2, vv;
  int n, np, i, ip, ipp, j, jp, jpp, m, mp;
  int mat_type, pmat_type, ir_part;
};
typedef struct integral_info_cuba_type cuba_t;
extern cuba_t g_cuba;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* distribution function diagnostics in phase space */
struct diagf_type
{
  int flag0;

  int sp_id, wid;

  double xx1, xx2, dxx;
  double th1, th2, dth;
  double dvol;

  double den, T;

  int nv;
  double vmax;
  double vs, dv;

  double *f;
};
typedef struct diagf_type diagf_t;

extern diagf_t g_fms0[NUM_SP_MAX][N_DGF_XX];
extern diagf_t g_fms1[NUM_SP_MAX][N_DGF_XX];
extern diagf_t g_fls0[NUM_SP_MAX][N_DGF_THE];
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
// UMFPACK complex type
struct umfz_type
{
  int m, nnz;
  int *Mp, *Mi;
  double *Mx, *Mz;
  void *Mn;
};
typedef struct umfz_type umfz_t;

// UMFPACK double type
struct umfr_type
{
  int m, nnz;
  int *Mp, *Mi;
  double *Mx;
  void *Mn;
};
typedef struct umfr_type umfr_t;

// HYPER IJ-type
struct hypi_type
{
  int m, m_loc, fst_row, lst_row, nnz;
  int *Mp, *Mi;
  double *Mx;

  int *rows_loc;
/*  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b, x;
  HYPRE_ParVector par_b, par_x;
  HYPRE_Solver solver, precond;*/
};
typedef struct hypi_type hypi_t;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* physical constants */
struct physical_const_type
{
  int flag0;

  /* physical values for normalization
  double ni00, Ti0, TI0, Te0;*/

  // physical values at plasma center
  double vti0, Omi0, rhoi0, rh0n, ompi0, lamD0, 
       rhs0n[NUM_SP_MAX], la0n, ls0, ls0n, beta0, alfv0;

  // physical values at mid-minor radius
  double vthi_mid, tau_t_mid, rhoi_mid, tau_b_mid,
       nuc_ii_mid, tc_ii_mid;

  /* ratios
  double mi_mI, mI_mi, mI_me, mi_me, me_mi, me_mI, 
       ZI_Zi, Zi_ZI, TI0_Ti0, Ti0_TI0;*/

  /* mass and charge number
  double Ms[NUM_SP_MAX], Zs[NUM_SP_MAX], 
       Ms_Zs[NUM_SP_MAX], Zs_Ms[NUM_SP_MAX], Zs_ng[NUM_SP_MAX];*/

  /* ion mass normalized by proton mass and charge numbers
  double mi, mI, Zi, ZI;*/

  // reference values for minimal density and temperature
  double den_min, tem_min;

  // basic constant
  double c, e, me, mp, eV, keV;

  // conversion factor to dimensional value (cgs)
  double cf_t, cf_x, cf_v, cf_Dif, cf_pfx, cf_tfx, cf_lfx,
       cf_Efld, cf_Pot, cf_pwrw, cf_res, cf_res1, cf_j, cf_j1;

};
typedef struct physical_const_type phc_t;

extern phc_t g_phc;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* mathematical and numerical constants */
struct const_type
{
  int flag0;

  double Zj[3], dZj[3];
  double epsilon, pi, twopi, twopi_inv, sqrt_pi;
};
typedef struct const_type c_t;

extern c_t g_c;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* mpi parallelization */
struct parallel_type
{ 
  int flag0;

  int ipe, Npe;
  int ipe_nlocal, Npe_nlocal;
  int ipe_fld, Npe_fld;

  MPI_Group grp_world;
  MPI_Group grp_nlocal;
  MPI_Group grp_fld;

  MPI_Comm  com_nlocal;
  MPI_Comm  com_fld;

  char fld_type;
};
typedef struct parallel_type pall_t;

extern pall_t g_pall;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* system variables */
struct system_type
{
  int flag0;

  // directory names
  char work_dir[2*CHARLEN], 
       snsp_dir[2*CHARLEN];

  // global file pointers
  FILE *fp_con, *fp_eng, *fp_traj;
  FILE *fp_phi00_his;
  FILE *fp_zonal_his, *fp_turbf_his;

  // seed for random number generator
  long rngs_seed;
};
typedef struct system_type sys_t;

extern sys_t g_sys;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* limiter information */
struct line_seg_type
{ 
  int p;
  double R1, Z1, R2, Z2, th1, th2;
};
typedef struct line_seg_type line_seg_t;

struct my_lim_type
{
  int flag0;

  line_seg_t *seg;
  int *idx;
};
typedef struct my_lim_type my_lim_t;

extern my_lim_t g_my_lim;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* GEQDSK data  */
struct eqdsk_type
{
  int flag0;

  int nw, nh, limitr, nbbbs;
  int NR, NZ, Npsi, ix[N_NULLX];

  double Rmin, Rmax, dR, Zmin, Zmax, dZ, psi0, psi1, dpsi, simag;
  double a0, Rc, Zc, psic, Bc, Ic, Rx[N_NULLX], Zx[N_NULLX], psix[N_NULLX];
  double *I_CS, *psi_2d_CS, *I_2d_CS, *rbbbs, *zbbbs, *rlim, *zlim;
  double rbbbs_max, rbbbs_min, zbbbs_max, zbbbs_min;
};
typedef struct eqdsk_type eqd_t;

extern eqd_t g_eqd;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* simple edge model */
struct simple_edge_type
{
  int flag0;

  double Rx_CS[N_FTH+3], Zx_CS[N_FTH+3],
       RI_CS[N_FTH+3], ZI_CS[N_FTH+3],
       epx_CS[N_FTH+3], epI_CS[N_FTH+3],
       nmfp_CS[N_FTH+3], Tn_CS[N_FTH+3];

  double nii_pdf_of_th_CS[N_FTH+3], 
       th_of_nii_Pf_CS[N_FTH+3];
};
typedef struct simple_edge_type edge_t;

extern edge_t g_edge;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* (R,Z) geometry and grid system */
struct RZ_grid_type
{
  int flag0;

  int NR, NZ, my_limitr;
  double Rmax, Rmin, Zmax, Zmin, Rx[N_NULLX], Zx[N_NULLX];
  double Rbmax, Rbmin, Zbmax, Zbmin;
  double dR, dZ, dR_inv, dZ_inv;
  double *psi_CS, *xx_CS, *I_CS, *B_CS, *bdB_CS, *J_B_CS, *I_R2B_CS;
  double lim_Rmin, lim_Rmax, lim_Zmin, lim_Zmax, *lim_R, *lim_Z;
  double psi_last;
};
typedef struct RZ_grid_type RZ_grid_t;

extern RZ_grid_t g_RZ;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* circular geometry */
struct cir_eq_type
{
  int flag0;

  double a0, R0, B0, q[5];
};
typedef struct cir_eq_type cir_t;

extern cir_t g_cir;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* analytic shaped geometry */
struct ase_eq_type
{
  int flag0;

  double R0, B0, a0, q0, kappa, del, shift;
};
typedef struct ase_eq_type ase_t;

extern ase_t g_ase;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* equilibrium profile */
struct eq_prof_type
{
  int flag0;

  // number of species for which profiles are set by input
  int nsp, op_qn, is_ion, is_ele;

  // indicator for multiple ionic species (imp, ep are counted)
  int nn_imp, *is_imp;

  sp_t sp[NUM_SP_MAX];

  double mass[NUM_SP_MAX], // mass in unit of proton or electron mass
       Ms[NUM_SP_MAX],   // normalized mass by main ion mass
       Zs[NUM_SP_MAX];   // charge in unit of e

  // central values of initial profiles
  double n00[NUM_SP_MAX], T00[NUM_SP_MAX], Pep0[NUM_SP_MAX];

  // initial profiles
  double *n0_CS[NUM_SP_MAX],  *T0_CS[NUM_SP_MAX],
       *Ll0_CS[NUM_SP_MAX], *Vl0_CS[NUM_SP_MAX];
  double *Pep0_CS[NUM_SP_MAX], *Cn0_CS[NUM_SP_MAX];
  double *Zeff0_CS;

  // time evolving profiles
  double *n_CS[NUM_SP_MAX],  *T_CS[NUM_SP_MAX], 
       *Ml_CS[NUM_SP_MAX], *Ll_CS[NUM_SP_MAX],
       *Vl_CS[NUM_SP_MAX], *Vlm_CS[NUM_SP_MAX];
  double *Pep_CS[NUM_SP_MAX], *Cn_CS[NUM_SP_MAX];
  double *Zeff_CS;
};
typedef struct eq_prof_type eprof_t;

extern eprof_t g_eprof;
/* ----------------------------------------------------------------------- */

  
/* ----------------------------------------------------------------------- */
/* (x, theta) geometry and grid system */
struct xth_grid_type
{
  int flag0;

  // radial direction
  int nxx,
      nxx1, // nxx + 1
      nxx2, // nxx + 2
      nxx3; // nxx + 3
  double xx0, xx1, *xx, dxx, dxx_inv, dxx2_inv;

  // poloidal direction
  int nth;
  double dth, dth_inv, dth2_inv;

  // grids on poloidal plain
  int nv;

  // toroidal direction
  int n1, n2,        // range of toroidal mode numbers [n1,n2]
      dn,            // increment of toroidal mode numbers
      nn;            // total number of toroidal mode numbers
  int my_n1, my_n2,  // range of my toroidal mode numbers [my_n1,my_n2]
      my_nn;         // increment of my toroidal mode numbers
  double zt1, zt1_inv;

  double xx_fd0, xx_fd1, xx_pd0[NUM_SP_MAX], xx_pd1[NUM_SP_MAX];
  double xx_out0, xx_out1;
  double xx_cchk0, xx_cchk1, dV_cchk;
  double ep_min, ep_max, epc, xxc;

  int tss_i;

  double *S2pi, *A, V_tot, nV_tot[NUM_SP_MAX], kapa;
  double edge_A_epxI, edge_V_epxI;
  double A_iL, V_iL;
  double *q_CS;
  dcomplex *thb;

  double xbry_R[N_FTH], xbry_Z[N_FTH];
  double dfth, dfth_inv, dfth2_inv;

  double *S_CS, *Sval, *epm_CS, *R_CS, *Z_CS, *eqR_CS, *eqZ_CS;
  double *Bmod_CS, *Bthmod_CS, *Bmin_CS, *Bmax_CS;
  double *Bth_B_xth_CS, *I_R2B_xth_CS;
  double **Qj_fx_CS, *Qj_fx;
  double *trho_CS, *xx_trho_CS;

  double xx_mid, ep_mid, q_mid, B_mid;
};
typedef struct xth_grid_type xth_t;

extern xth_t g_xth;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* field solvers and variables */
struct field_type
{
  int flag0;

  /* flag for adibatic response correction
     by substracting flux surface averaged quantity */
  int flag_adc;

  // axisymmetric field and source
  double *phi00, *den00;
  double *All00, *jll00;
  dcomplex *phi, *den;
  dcomplex *All, *jll;

  umfz_t *p_sol;    // Poisson solver
  umfz_t *r_sol;    // FRL correction matrix
  umfz_t *s_sol;    // Smoothing solver

  umfz_t zp_sol2;   // 2D zonal solver

  umfr_t zp_sol1;   // 1D zonal solver     (double matrices !!!)  
  umfr_t zr_sol1;   // 1D zonal FRL matrix (double matrices !!!)  
  umfr_t zs_sol1;   // 1D zonal smoother   (double matrices !!!)  

#ifdef POL_FFT
  int op_pfft_on;
#endif
};
typedef struct field_type field_t;

extern field_t g_fld;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* temporary variables */
struct tmp_type
{
  int flag0;

  int nn;
  double *buf_re, *buf_im, *tot_re, *tot_im;
  double *bx, *bz, *xx, *xz;
  double dgf[N_DGF_VTOT];
};
typedef struct tmp_type tmp_t;

extern tmp_t g_tmp;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* my-quadrature types */
struct my_qd5_type
{
  double t[5], w[5];
  double x2D[25], y2D[25], w2D[25];
};
typedef struct my_qd5_type my_qd5_t;

extern my_qd5_t g_my_qd5;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* damping profiles */
struct damp_type
{
  int flag0;
  double *edp_prof;
};
typedef struct damp_type damp_t;

extern damp_t g_damp;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* Coulomb collision */
struct coll_type
{
  int flag0;
  
  double *nuc_ii, *vcs0, *cprof;
  double nu0[NUM_SP_MAX][NUM_SP_MAX];
  double *den0[NUM_SP_MAX], *tem0[NUM_SP_MAX],
       *vth0[NUM_SP_MAX], *Vll0[NUM_SP_MAX];

  double *Vl0_2d[NUM_SP_MAX];

  double *cln[NUM_SP_MAX][NUM_SP_MAX];

  int *cnth, *acnth;

  double *dth, *dth_inv, *ep_cir, v_cmin0, *ndv;

  double st12, st1pi, stf2, ap[5];

  double *buf_tmp, *buf_tot;

  int itcd0, itcd1, iccd0, iccd1;

  double *vl0, *vv0, *wt0;
  double *dvl[NUM_SP_MAX][NUM_SP_MAX];
  double *dvv[NUM_SP_MAX][NUM_SP_MAX];
  double *dwt[NUM_SP_MAX];

  double *s_dvl_Ie, *s_dvl_Ii, *s_dvv_Ie, *s_dvv_Ii;
};
typedef struct coll_type coll_t;

extern coll_t g_coll;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* common diagnostics */
struct diag_type
{
  int flag0;

  int flag_Er;

  // geometric infomration for diagnostics
  double thm1, thm2;
  double *ep, *xx, *q, *R, *dV, *dVm, *fsa,
       *dp2, *dpdR, *dxdR, *dpdx, *dxdp,
       *Bm, *B2, *B2inv, *Bp_B, *epd, *R_om, *B_om, *Bp_om;

  // reference values for diagnostic
  double *n0[NUM_SP_MAX], *T0[NUM_SP_MAX], *Ll0[NUM_SP_MAX], *Vl0[NUM_SP_MAX],
       *dndp0[NUM_SP_MAX], *dTdp0[NUM_SP_MAX],
       *dLldp0[NUM_SP_MAX], *dVldp0[NUM_SP_MAX];
  double *XiCH0, *XgB0[NUM_SP_MAX],
       *XB0[NUM_SP_MAX], *cs0[NUM_SP_MAX];

  // diagnostic variables
  double *pfx[NUM_SP_MAX], *tfx[NUM_SP_MAX], *lfx[NUM_SP_MAX], *vfx[NUM_SP_MAX],
       *T4fx[NUM_SP_MAX], *V4fx[NUM_SP_MAX],
       *n[NUM_SP_MAX], *T[NUM_SP_MAX], *Ll[NUM_SP_MAX], *Vl[NUM_SP_MAX],
       *Vlm[NUM_SP_MAX], *vis_xx_ll[NUM_SP_MAX], *ttorq[NUM_SP_MAX];

  // 2d diagnostics for collision
  double *Vl_2d[NUM_SP_MAX];

  // instantaneous field
  double *E00, *P00, *P00_rhk;
  double *A00, *dA00;
  // time averaged field
  double *phi00neo, *E00neo;

  // file pointers
  FILE *fp1[NUM_SP_MAX];
  FILE *fpfx[NUM_SP_MAX];
  FILE *fp10[NUM_SP_MAX];
  FILE *fp10_2d[NUM_SP_MAX];

  // counter for average 
  int cnt_t1,   max_t1;
  int cnt_fx,   max_fx;
  int cnt_t10,  max_t10;
  int cnt_coll, max_coll;
  int cnt_Eneo, max_Eneo;
  int cnt_rf,   max_rf;
  int cnt_fmid, max_fmid;
  int cnt_fms[NUM_SP_MAX], max_fms[NUM_SP_MAX];
  int cnt_fls[NUM_SP_MAX], max_fls[NUM_SP_MAX];
  int cnt_phs[NUM_SP_MAX], max_phs[NUM_SP_MAX];
};
typedef struct diag_type diag_t;

extern diag_t g_dg;
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* bounce-averaged kinetic electrons */
struct badk_tb_type
{
  double xx0, xx1, dxx, dxx_inv;

  double *vte, *v0, *v1, xi0, xi1;
  double *dv, *dv_inv, dxi, dxi_inv;

  /*double *vte, *vl0, *vl1, *vp0, *vp1;
  double *dvl, *dvp, *dvl_inv, *dvp_inv;*/

  double *tb, *trp, *thb;
};
typedef struct badk_tb_type badk_tb_t;
struct badk_type
{
  int flag0;

  int op_pp, pp_ns_tot, op_pass_wgt;

  /*int op_nbadk;*/
  int op_wbadk;

  double dbth, dbth_inv, dtp;

  double ptl_baw[2*BADK_NBA_MAX+2], ptl_Bnc[2*BADK_NBA_MAX+2],
       ptl_thb[2*BADK_NBA_MAX+2], ptl_ztb[2*BADK_NBA_MAX+2];

  double ptl_tb, ptl_tb_dg_phs;

  double *nu_CS;

  double *trap_frac0_CS, *trap_frac_CS,
       *Bom_CS, *Rom_CS, *ep_CS, *q_CS, *ff_CS, 
       *Cpne_Te0_CS, *Cpne_Te_CS;
#if defined(SEO_MOD)
  int trap_frac_2D_intp_bd, trap_frac_2D_intp_nth;
  double trap_frac_2D_intp_th0, trap_frac_2D_intp_th1;
  double *trap_frac_2D_CS, *Cpne_Te_2D_CS;
#endif

  dcomplex *den;

  double *den00;

  double *thb_up_CS;
  double *thb_down_CS;

  badk_tb_t tb;

  ptlv_t *ele_array;
};
typedef struct badk_type badk_t;

extern badk_t g_badk;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* phase space diagnostics */
struct phase_type
{
  double xx0, xx1, dxx;
  double th0, th1, dth;
  double dvol, den, T, Jx;
  double vth, dvl, dvp, vl0, vl1, vp0, vp1;
  double *fv;
};
typedef struct phase_type phs_t;

struct diag_phase_type
{
  int flag0;

  double *tmp_buf;
  phs_t phs[NUM_SP_MAX][N_DG_PHS_XX*N_DG_PHS_TH];
  double xx0, xx1, dxx, dxx_inv;
  double dth, th0[N_DG_PHS_TH], th1[N_DG_PHS_TH];
  double dzt, zt0, zt1;
};
typedef struct diag_phase_type diag_phs_t;

extern diag_phs_t g_dg_phs;
/* ----------------------------------------------------------------------- */

#ifdef MOD_QD
int g_my_qd_order;
double *g_my_qd_t, *g_my_qd_w, *g_my_qd_x2D, *g_my_qd_y2D, *g_my_qd_w2D;
#endif

#ifdef SEO_MOD
int g_add_init_pot_on, g_add_init_den00_on;
double g_add_init_pot_wave_num, g_add_init_pot_x_st, g_add_init_pot_x_ed;
double g_add_init_pot_xx;
double g_add_init_pot_dpotdx_ref, g_add_init_pot_dpotdx, g_add_init_pot_del_den_ref;
double g_add_init_pot_ep_st, g_add_init_pot_ep_ed, g_add_init_pot_ep;
#endif

/* ----------------------------------------------------------------------- */
/* 3D field diagnostics */
struct diag_3df_type
{
  int flag0;

  dcomplex *den[NUM_SP_MAX], *tem[NUM_SP_MAX], *vll[NUM_SP_MAX];
};
typedef struct diag_3df_type diag_3df_t;

extern diag_3df_t g_diag_3df;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* zonal preserving noise reduction */
struct nrdc_zpk_type
{
  int flag0;

  int nxx, nE, nthfs, nv, nl, nlnv, nv_1, nl_1;
  double dxx, dxx_inv, dthfs, dthfs_inv, dE, dE_inv, E0, E1;
  double v0[NUM_SP_MAX], v1[NUM_SP_MAX], dv[NUM_SP_MAX], dv_inv[NUM_SP_MAX];
  double vl0[NUM_SP_MAX], vl1[NUM_SP_MAX], dvl[NUM_SP_MAX], dvl_inv[NUM_SP_MAX];
  double *gs, *gz, *gn, *ge, *gh;
  double *vll_B[NUM_SP_MAX], *src, *mat, *tsrc, *tmat, *hc1;
  double *sum_df, *sum_f0, *tsum_df, *tsum_f0;
  double *itgl_df_f0, *itgl_df, *itgl_f0;
  double *df3_fsa, *df3_num, *df3_tmp;
#ifdef MOD_IMP
  double *w_ratio[NUM_SP_MAX], *v_th[NUM_SP_MAX];
#endif
};
typedef struct nrdc_zpk_type zpk_t;

extern zpk_t g_zpk;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* fluid passing electrons */
struct efluid_type
{
  int flag0;

  /* time integration algorithm for hybrid scheme */
  int op_tint;

  /* further parallelization parameters */
  int nnpe, myc_nxy, myc_ij0, myc_ij1;
  int nidx_nnp, nidx_ten1, nidx_ten2, nidx_ten, nn_bvec, nn_vv;

  /* parallel diffusion parameter */
  double Dll;

  /* Te vector for pressure calculation */
  double *Te0;

  /* dAlldt for vll formulation */
  dcomplex *dAlldt;

  dcomplex *vec0[N_EF_VECTOR], *vec[N_EF_VECTOR], *vecp[N_EF_VECTOR];
  dcomplex *tmp_vec[N_EF_VECTOR];

  /* dAll_Ld for landau damping */
  dcomplex *dAll_Ld;

  /* ex rk4 */
  dcomplex *dvdt0[N_EF_TVECTOR], *dvdt1[N_EF_TVECTOR],
	   *dvdt2[N_EF_TVECTOR], *dvdt3[N_EF_TVECTOR];

  /* coupling with bounce-averaged trapped electrons */
  dcomplex *rden, *dne_t, *dPel_t, *dPep_t;

  /* imex */
  // fluid vectors
  dcomplex *U0[EF_IMEX_MAX][N_EF_TVECTOR], 
	   *U1[EF_IMEX_MAX][N_EF_TVECTOR];
  // particle vectors
  double *F[EF_IMEX_MAX];
  // IMEX table
  double aa[EF_IMEX_MAX][EF_IMEX_MAX], 
       at[EF_IMEX_MAX][EF_IMEX_MAX],
       ww[EF_IMEX_MAX], wt[EF_IMEX_MAX];

  /* M[n][i,j][di,dj] */
  dcomplex *mat[N_EF_MATRIX];

  /* T[n][i,j][np][di,dj][dii,djj] */
  dcomplex *ten[N_EF_TENSOR];

  /* Big matrix for implicit time integration */
  double *big_sol_xx, *big_sol_xz, *big_src_xx, *big_src_xz;
  int *BM_nnz, **BMp, **BMi;
  double **BMx, **BMz;
  void **BMn;

  /* HYPRE types for Big matrix solver */
  int sid, Nh_tot, Nh_loc, Ih0, Ih1;
  double *b_hyp, *x_hyp, *sol_hyp;
  hypi_t BM_hyp;

  /* Maximum matrix size for diagnostic purpose */
  dcomplex mat_max[N_EF_MATRIX], ten_max[N_EF_TENSOR], pmat_max[N_PADE_MAT];
  double mat_amax[N_EF_MATRIX], ten_amax[N_EF_TENSOR], pmat_amax[N_PADE_MAT];

  double *res_CS, *j0_CS;
};
typedef struct efluid_type ef_t;

extern ef_t g_ef;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* read experimental profiles */
struct read_eprof_type
{
  int flag0;

  // flag to indicate which profiles to read
  int flag[NUM_SP_MAX][DEF_RPF_NUM];

  // conversion factor to be multiplied
  double cf[NUM_SP_MAX][DEF_RPF_NUM];

  // conversion factor for xx-coordinate
  double cf_xx[NUM_SP_MAX]; 

  // number of data points
  int nep[NUM_SP_MAX][DEF_RPF_NUM];

  double ep0[NUM_SP_MAX][DEF_RPF_NUM],
       ep1[NUM_SP_MAX][DEF_RPF_NUM],
       *ep_CS[NUM_SP_MAX][DEF_RPF_NUM];
};
typedef struct read_eprof_type read_eprof_t;

extern read_eprof_t g_read_eprof;
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */
/* input variables */
struct inpv_type
{
  int flag0;

  /* simulatin kind */
  int op_sim_kind;

  /* physics options */
  int op_neo, op_zonal, op_turbf, op_em_sim;

  /* radial domain */
  int nxx, nth;
  double xx0, xx1;

  /* range of self-consistent field */
  int ifd0, ifd1;

  /* boundary option */
  int op_lim, op_recycl;
  double rrl_cir;

  /* collision options */
  int op_coll, skip_coll, sub_coll, sub_colle, op_coll_prof;
  double coll_prof_A0;

  /* noise control and source */
  int op_nrdc_zpk, nthfs_zpk;
  double gs_zpk, gh_zpk;

  /* damping */
  int op_edp, op_hypv, op_wdp, op_Dll;
  double edp_dxx, edp_alpha, edp0, hypL, wdpx, Dll;

  /* smoothing zonal */
  int nxx_sm00;
  double xx0_sm00, xx1_sm00;

  /* neoclassical Er */
  int op_neo_Er;
  double efac_neo_Er;

  /* toroidal spectrum */
  int tor_nn, tor_dn, tor_n1, tor_n2;

  /* time integration */
  int op_tint, tot_nstep;
  double dt;

  /* diagnostic freq */
  int dg_phi_cnt, dg_t1_cnt, dg_fx_cnt, dg_t10_cnt;
  int op_dg_wgt, skip_dg_wgt;

  /* -------------------------------------------------------------------- */
  /* input variables for profiles and marker particles for each species */
  int nsp, op_qn, op_eqp_model;

  int op_load_ef;
  int op_load_ptl[NUM_SP_MAX];

  sp_t sp[NUM_SP_MAX];
  double mass[NUM_SP_MAX], charge[NUM_SP_MAX];
  int ng[NUM_SP_MAX];

  int nptl[NUM_SP_MAX];
  int ipd0[NUM_SP_MAX], ipd1[NUM_SP_MAX], ipd2[NUM_SP_MAX];
  double kv[NUM_SP_MAX];

  int op_eqp_n[NUM_SP_MAX];
  double eqp_n00[NUM_SP_MAX], eqp_Ln[NUM_SP_MAX], 
       eqp_rc_n[NUM_SP_MAX], eqp_dr_n[NUM_SP_MAX];

  int op_eqp_T[NUM_SP_MAX];
  double eqp_T00[NUM_SP_MAX], eqp_LT[NUM_SP_MAX], 
       eqp_rc_T[NUM_SP_MAX], eqp_dr_T[NUM_SP_MAX];

  int op_eqp_V[NUM_SP_MAX];
  double eqp_LV[NUM_SP_MAX], eqp_rc_V[NUM_SP_MAX], eqp_dr_V[NUM_SP_MAX],
       eqp_V0[NUM_SP_MAX], eqp_V1[NUM_SP_MAX];

  int op_eqp_Pep[NUM_SP_MAX];
  double eqp_Pep0[NUM_SP_MAX], eqp_L_Pep[NUM_SP_MAX],
       eqp_rc_Pep[NUM_SP_MAX], eqp_dr_Pep[NUM_SP_MAX];

  double eqp_lambda0[NUM_SP_MAX], eqp_del_lambda[NUM_SP_MAX],
       eqp_Emax[NUM_SP_MAX];
  /* -------------------------------------------------------------------- */

  /* geometry option */
  int op_geo_info;
  double cir_q0, cir_q1, cir_q2, cir_q3, cir_q4, cir_a0, cir_R0, cir_B0;
  double ase_R0, ase_B0, ase_a0, ase_q0, ase_kappa, ase_del, ase_shift;
  
  /* option to load auxilary info and snapshot */
  int op_RZ_xfth_file,
      op_mtrxf, op_ef_matf, op_ef_tenf, 
      op_ef_BMf, op_ef_LDf,
      op_test_output,
      op_snsp;



  /* -------------------------------------------------------------------- */
  /* input variables for experimental equilibrium and profiles */
  int eqd_nw, eqd_nh;
  int flag_rpf[NUM_SP_MAX][DEF_RPF_NUM];
  double cf_rpf[NUM_SP_MAX][DEF_RPF_NUM];
  double cf_rpf_xx[NUM_SP_MAX];
  /* -------------------------------------------------------------------- */



  /* -------------------------------------------------------------------- */
  /* internal input variables set by combinations of external ones */
  int op_1d_zonal;
  int mixcoord[NUM_SP_MAX];
  int op_dg_fms[NUM_SP_MAX];
  int op_dg_phs[NUM_SP_MAX];
  int op_dg_fls[NUM_SP_MAX];
  /* -------------------------------------------------------------------- */

#ifdef SEO_ADD 
  double delr_ep_den, delr_ep_temp;
  int fft_half_m_width;
#endif

};
typedef struct inpv_type inpv_t;

extern inpv_t g_inpv;
/* ----------------------------------------------------------------------- */





/* ----------------------------------------------------------------------- */
/* external constants */
extern const double BADK_DTHB_MAX;
extern const double BADK_DTHB_MIN;

extern const int WT_Imp;

extern const int N_EINFO;

extern const int N_CS_ITER;

extern const double Kll_SIGN;
extern const double CUR_DIR;
extern const double BP_SIGN;
extern const double BT_SIGN;
extern const double BPT_SIGN;

extern const int DEF_CUBA_VERBOSE;
extern const int DEF_CUBA_LAST;
extern const int DEF_CUBA_MINEVAL;
extern const int DEF_CUBA_MAXEVAL;
extern const double DEF_CUBA_REL_ERR;
extern const double DEF_CUBA_ABS_ERR;
extern const int DEF_CUBA_KEY;

extern const int PHS_VARS_F0;
extern const int PHS_VARS_F1;
extern const int PHS_VARS_PFX;
extern const int PHS_VARS_TFX;
extern const int PHS_VARS_MFX;
extern const int N_DG_PHS_VAR;
extern const double DG_PHS_VMAX;
extern const int N_DG_PHS_VL;
extern const int N_DG_PHS_VP;
extern const int N_DG_PHS_FTOT;

extern const int IDX_RPF_DEN;
extern const int IDX_RPF_TEM;
extern const int IDX_RPF_VLL;

extern const int DGRT_STEP; // 0
extern const int DGRT_ION;  // 1
extern const int DGRT_TELE; // 2
extern const int DGRT_PELE; // 3
extern const int DGRT_FLD;  // 4
extern const int DGRT_COLL; // 5
extern const int DGRT_EFT;  // 6
extern const int DGRT_EFM;  // 7
/* END of gvar-srt.c ----------------------------------------------------------------------- */



/* START of func-srt.h ================================================================== */

/* dxdt-RZ.c */
void dxdt_RZ_vll(stime_t *stime, ptlv_t *pp, int sp_id, double *dxdt);

/* dxdt-xth.c */
void dxdt_xth_vll(stime_t *stime, ptlv_t *pp, int sp_id, double *dxdt);


/* fem-func.c */
void cal_S2c(double x, double x0, double dx, int *i, double *w);
void cal_QSw_xth(double xx, double th, int *i, int *j, double *wi, double *wj);
double Qx(double xx, int i);
double Qt(double th, int j);
double Lx_fn(double xx, int i);
double Lt_fn(double th, int j);
double S0x(double xx, int j);
double S0t(double th, int j);
double dQx(double xx, int i);
double ddQx(double xx, int i);
double dQt(double th, int j);
double ddQt(double th, int j);


/* func.c */
void set_my_constants(void);
void reset_flag0(void);
void set_basic_paras(void);
void set_tor_spectrum(void);
void select_simulation_kind(void);
void set_mid_radius_info(void);
void set_stime(stime_t *stime);
void allocate_tmp_array(void);
void check_markers(int nntime, int ppstep);


/* geom-ase.c */
void set_ase_efit_data(void);
double ase_psi_fn_RZ(double R, double Z);


/* geom-cir.c */
void set_cir_efit_data(void);
double cir_q_fn_ep(double ep);
double cir_q_afn_ep(double ep);
double cir_dqde_afn_ep(double ep);
double cir_dpde_afn_ep(double ep);
double cir_psit_afn_ep(double ep);
double cir_psi_fn_ep(double ep);

/* geqdsk.c */
void read_geqdsk_and_test(void);
void read_geqdsk_file(FILE *fp, int eqd_nw, int eqd_nh);
void my_geqdsk_to_cgs(void);
void refine_my_geqdsk(void);
void write_my_geqdsk(void);
int comp_grid_pt(const void *x, const void *y);
void sort_grid_pt(grid_pt_t *curve, int num);
double nrc_wrap_eqd_psi_fn(double *x);
void nrc_wrap_eqd_dpsi_fn(double *x, double *dy);
double nrc_wrap_eqd_gdp2_fn(double *x);
void nrc_wrap_eqd_dgdp2_fn(double *x, double *dy);
void find_min_gdp2_along_bbbs(int istart, int num,
                              int *nmin, double *psi, double *gdp2);

/* limiter.c */
void simplify_RZ_curve(double *cv_R, double *cv_Z, int num, double err, int *rnum);
void init_my_limiter(void);
int check_limiter(ptlv_t *pp, double *Rh, double *Zh);

/* RZ-grid.c */
void set_RZ_simulation_domain(void);
void set_limiter(void);

/* xth-grid.c */
void setup_xx_theta_grid(void);
void set_plasma_boundary(void);
void set_epm_R_Z_xfth(void);
void find_psi_RZ(double psi, double R_s, double Z_s, double err, double *R, double *Z);
void S_fth_CS(int i, double th, double *S, double *dSdt);
void S_fth_CS2(int i, double th, double *Si, double *dSi, double *ddSi);
void init_S_fth_CS(double *data);
void init_Qj_fx(void);
void intp_all_Qj_fx(double x, double *Qjfx);
void test_grid_setup(void);

#ifdef POL_FFT
void init_fft_mode_filter(void);
void pfft_set_alpha_intp(void);
void pfft_al_to_th_raw(double xx_in, double al_in, double *th_out);
void set_pfft_mat(void);

double cal_Z_ii_mat_fn(int i_in, int ii_in);
double fn_Z_ii_mat_arg(double xx_in);

void cal_ptou_mat(int ii_in, int jj_in, int i_in, int m_in, int loc_n, double *val_re, double *val_im);
void Integrand_func_ptou_mat(const int *ndim, const double *zz, const int *ncomp, double *ff);
void cal_utos_mat(int ii_in, int jj_in, int i_in, int m_in, int loc_n, double *val_re, double *val_im);
void Integrand_func_utos_mat(const int *ndim, const double *zz, const int *ncomp, double *ff);
void cal_ptop_mat(int i_in, int j_in, int jj_in, int loc_n, double *val_re, double *val_im);
void Integrand_func_ptop_mat(const int *ndim, const double *zz, const int *ncomp, double *ff);

void th_al_conversion_output(void);
void test_pol_filter_mod(void);
void fft_filter(dcomplex *src, int flag_input);

void pfft_theta_to_dal(double xx_in, double th_in, double *daldx_out, double *daldth_out, int op_in);

double fn_theta_to_alpha(double xx, double th_in, int op_in);
double fn_alpha_to_theta(double xx, double al_in, int op_in);

void S_CS_ij_at_th(int i, int j_in, double th_in, double *S_out, int op_in);
double Qt_pol_filter(double th_in, int j_in);

void Landau_damping_fft(double dtCL, dcomplex *jll, dcomplex *All);
void fft_kpara_diag(dcomplex *All);
#endif

/* nrc-bessj.c */
double bessj0(double x);
double bessj1(double x);
double bessj(int n, double x);


/* nrc-miz.c */
void frprmn(double *p, int n, double ftol, int *iter, double *fret,
            double (*func)(double []), void (*dfunc)(double [], double []));


/* nrc-mtdc.c */
void my_choldc(dcomplex *m, int n, double *p, int fop);
void my_cholsl(dcomplex *m, int n, double *p, dcomplex *b, dcomplex *x);
void cludcmp(dcomplex *a, int n, int *indx, double *d);
void clubksb(dcomplex *a, int n, int *indx, dcomplex *b);
void ludcmp(double *a, int n, int *indx, double *d);
void lubksb(double *a, int n, int *indx, double *b);
/* nrc ===================================================================== */

/* my-spline =============================================================== */
/* cs-1d-util.c */
double intp_1d_CS(double x, double x1, double x2, int nx, double *array);
void intp_dfn_1d_CS(double x, double x1, double x2, int nx, double *array, double *val);
void init_intp_1d_CS(double x1, double x2, int nx, double *data, int bc);

/* cs-2d-util.c */
double intp_2d_CS(double x, double y,
                double x1, double x2, int nx,
                double y1, double y2, int ny, double *array);
void intp2_2d_CS(double x, double y,
                 double x1, double x2, int nx,
                 double y1, double y2, int ny, 
		 double *array1, double *array2, double *val1, double *val2);
void intp_dfn_2d_CS(double x, double y,
                    double x1, double x2, int nx,
                    double y1, double y2, int ny,
                    double *array, double *val);
void init_intp_2d_CS(double x1, double x2, int nx,
                     double y1, double y2, int ny,
                     double *data, int bc_x, int bc_y);


/* cs-xth-util.c */
double intp_xth_CS(double x, double y,
                 double x1, double x2, int nx,
                 double y1, double y2, int nt,
                 double *array);
void intp_dfn_xth_CS(double x, double y,
                     double x1, double x2, int nx,
                     double y1, double y2, int nt,
                     double *array, double *val);


/* qs-1d-util.c */
double intp_1d_QS(double x, double x1, double x2, int nn, double *data);
void init_intp_1d_QS(double x1, double x2, int nn, double *data, int bc);


/* spline.c */
void cal_Qss(double x, double x0, double dx_inv, 
             int *i, double *Q, double *dQ, double *ddQ);
void cal_Css(double x, double x0, double dx_inv, 
             int *i, double *C, double *dC, double *ddC);
void cal_Qs(double t, double dx_inv, double *Qs);
void cal_Cs(double t, double dx_inv, double *Cs);


/* umfpack-util.c */
void constr_umfpackM_1d_CS(int nn, int **Ap, int **Ai, double **Ax, void **Num);
void constr_umfpackM_1d_QS(int nn, int **Ap, int **Ai, double **Ax, void **Num);
void init_intp_1d_CS_umfpack(double x1, double x2, int nn, double *data);
void init_intp_1d_QS_umfpack(double x1, double x2, int nn, double *data);
void constr_umfpackM_2d_CS(int nx, int ny,
                           int **Ap, int **Ai, double **Ax, void **Num);
void init_intp_2d_CS_umfpack(double x1, double x2, int nx,
                             double y1, double y2, int ny, double *data);
void constr_umfpackM_xth_CS(int nx, int nt,
                            int **Ap, int **Ai, double **Ax, void **Num);
void init_intp_xth_CS_umfpack(double x1, double x2, int nx,
                              double y1, double y2, int nt, double *data);
/* my-spline =============================================================== */

/* diag-gen.c */
void set_dg0_vars(stime_t *stime);
void set_dg0_SOL(stime_t *stime);
void my_reduce_dg_rbuf(double *buf);
void my_allreduce_dg_rbuf(double *buf);

/* inpv.c // */
void check_input(void);
void read_input_file(void);

/* START of cquadpak.h  */


#define uflow 	DBL_MIN
#define oflow 	DBL_MAX
#define epmach 	DBL_EPSILON
#define LIMIT   10000	
/*#define LIMIT 	500*/
#define MAXP1 	21
#define Pi      M_PI
#define COSINE 	1
#define SINE	2

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif

/* Integration routines */
/* Gauss-Kronrod for integration over finite range. */
double G_K15(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K21(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K31(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K41(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K51(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K61(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);

/* Gauss-Kronrod for integration over infinite range. */

/* Gauss-Kronrod for integration of weighted function. */
double dqext(int *n,double epstab[],double *abserr,
	double res3la[],int *nres);
void dqsort(int limit,int last,int *maxerr,double *ermax,
	double elist[],int iord[],int *nrmax);
double dqagi(double f(),double bound,int inf,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqags(double f(),double a,double b,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqagp(double f(),double a,double b,int npts2,double *points,
	double epsabs,double epsrel,double *abserr,int *neval,int *ier);
double dqng(double f(),double a,double b,double epsabs,double epsrel,
	double *abserr,int *neval,int *ier);
double dqag(double f(),double a,double b,double epsabs,double epsrel,
	int irule,double *abserr,int *neval,int *ier);
double dqage(double f(),double a,double b,double epsabs,double epsrel,
    int irule,double *abserr,int *neval,int *ier,int *last);
double dqwgtc(double x,double c,double p2,double p3,double p4,
	int kp);
double dqwgto(double x,double omega,double p2,double p3,double p4,
	int integr);
double dqwgts(double x,double a,double b,double alpha,double beta,
	int integr);
void dqcheb(double *x,double *fval,double *cheb12,double *cheb24);
double dqc25o(double f(),double a,double b,double omega,int integr,
	int nrmom,int maxp1,int ksave,double *abserr,int *neval,
	double *resabs,double *resasc,int *momcom,double **chebmo);		
double dqfour(double f(),double a,double b,double omega,int integr,
    double epsabs,double epsrel,int icall,int maxp1,
    double *abserr,int *neval,int *ier,int *momcom,
	double **chebmo);
double dqawfe(double f(),double a,double omega,int integr,double epsabs,
    int limlst,int maxp1,double *abserr,int *neval,int *ier,
    double *rslst,double *erlist,int *ierlst,double **chebmo);
double dqawf(double f(),double a,double omega,int integr,double epsabs,
	double *abserr,int *neval,int *ier);
double dqawo(double f(),double a,double b,double omega,int integr,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqaws(double f(),double a,double b,double alfa,double beta,int wgtfunc,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier);
double dqawse(double f(),double a,double b,double alfa,double beta,
    int wgtfunc,double epsabs,double epsrel,double *abserr,
    int *neval,int *ier);
void dqmomo(double alfa,double beta,double ri[],double rj[],double rg[],
    double rh[],int wgtfunc);
double dqc25s(double f(),double a,double b,double bl,double br,double alfa,
    double beta,double ri[],double rj[],double rg[],double rh[],
    double *abserr,double *resasc,int wgtfunc,int *nev);
double dqc25c(double f(),double a,double b,double c,double *abserr,
    int *krul,int *neval);
double dqawc(double f(),double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier);
double dqawce(double f(),double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier);

double G_B15(double f(),double a,double b,double *abserr,
	double *resabs, double *resasc);



