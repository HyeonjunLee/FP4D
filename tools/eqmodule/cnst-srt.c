#include "cnst-srt.h"
#include <math.h>


const double BADK_DTHB_MAX = 0.3;
const double BADK_DTHB_MIN = 0.05;

// set this to 5(6) for df(ff) impurity
const int WT_Imp = 5; 

/* number of edge-interpolation variables */
const int N_EINFO = 8;

const int N_CS_ITER = 3000;

const double Kll_SIGN = 1.0;
const double CUR_DIR  = 1.0;
const double BP_SIGN = 1.0;
const double BT_SIGN  = 1.0;
const double BPT_SIGN = 1.0;

const int DEF_CUBA_VERBOSE = 2;
const int DEF_CUBA_LAST = 4;
const int DEF_CUBA_MINEVAL = 1;
const int DEF_CUBA_MAXEVAL = 10000;
/*const int DEF_CUBA_MINEVAL = 100;
  const int DEF_CUBA_MAXEVAL = 100000;*/
#ifdef DEF_ETG_SIM
const double DEF_CUBA_REL_ERR = 1.0e-8;
const double DEF_CUBA_ABS_ERR = 1.0e-15;
#else
const double DEF_CUBA_REL_ERR = 1.0e-7;
const double DEF_CUBA_ABS_ERR = 1.0e-12;
#endif
const int DEF_CUBA_KEY = 0;

/* for diag-phs.c */
const int PHS_VARS_F0  = 0;
const int PHS_VARS_F1  = 1;
const int PHS_VARS_PFX = 2;
const int PHS_VARS_TFX = 3;
const int PHS_VARS_MFX = 4;
const int N_DG_PHS_VAR = 5;
const double DG_PHS_VMAX  = 3.5;
const int N_DG_PHS_VL = 80;
const int N_DG_PHS_VP = 40;
/* this should be set as N_DG_PHS_VAR*(N_DG_PHS_VL+1)*(N_DG_PHS_VP+1) */
const int N_DG_PHS_FTOT = 16605;

/* for read-eprof.c */
const int IDX_RPF_DEN = 0;
const int IDX_RPF_TEM = 1;
const int IDX_RPF_VLL = 2;

/* for runtime.c */
const int DGRT_STEP = 0;
const int DGRT_ION  = 1;
const int DGRT_TELE = 2;
const int DGRT_PELE = 3;
const int DGRT_FLD  = 4;
const int DGRT_COLL = 5;
const int DGRT_EFT  = 6;
const int DGRT_EFM  = 7;

/* START of gvar-srt.c */
/* input variable group */
inpv_t g_inpv;


/* mpi parallelization */
pall_t g_pall;


/* physical constants */
phc_t g_phc;


/* numerical constants */
c_t g_c;


/* system variables */
sys_t g_sys;


/* integration */
qdpk_t g_qdpk;
cuba_t g_cuba;


/* marker */
pl_t g_ptl;


/* time steps */
stime_t g_stime;


/* limiter information */
my_lim_t g_my_lim;


/* geqdsk data */
eqd_t g_eqd;


/* simple edge model */
edge_t g_edge;


/* (R,Z)-grid system */
RZ_grid_t g_RZ;


/* circular geometry */
cir_t g_cir;


/* analytic shaped geometry */
ase_t g_ase;


/* useful mid-radius information 
   real g_xx_mid, g_ep_mid, g_q_mid, g_B_mid,
   g_vthi_mid, g_tau_t_mid, g_rhoi_mid, g_tau_b_mid, 
   g_nuc_ii_mid, g_tc_ii_mid;*/

/* equilibrium profile */
eprof_t g_eprof;


/* (xx,theta)-grid system */
xth_t g_xth;


/* field and solver */
field_t g_fld;


/* temporary storage */
tmp_t g_tmp;


/* my-quadratures */
my_qd5_t g_my_qd5;


/* damping profiles */
damp_t g_damp;


/* collision related */
coll_t g_coll;


/* common diagnostics */
diag_t g_dg;


/* diag-f */
diagf_t g_fms0[NUM_SP_MAX][N_DGF_XX];
diagf_t g_fms1[NUM_SP_MAX][N_DGF_XX];
diagf_t g_fls0[NUM_SP_MAX][N_DGF_THE];


/* BADK module */
badk_t g_badk;


/* phase space diagnostics */
diag_phs_t g_dg_phs;


/* 3D diagnosis */
diag_3df_t g_diag_3df;


/* noise reduction */
zpk_t g_zpk;


/* fluid passing electrons */
ef_t g_ef;

/* read experimental profiles */
read_eprof_t g_read_eprof;


/* START of func.c */


void set_my_constants_(void)
{
  /* general constants --------------------------------------------------- */
  g_c.epsilon = 1.0e-34;

  g_c.pi        = M_PI;
  g_c.twopi     = 2.0*M_PI;
  g_c.twopi_inv = 1.0/g_c.twopi;
  g_c.sqrt_pi   = double_sqrt(g_c.pi);

  g_c.Zj[0]  = 0.0;
  g_c.Zj[1]  = 1.0;
  g_c.Zj[2]  = 0.0;
  g_c.dZj[0] = 0.0;
  g_c.dZj[1] = 0.0;
  g_c.dZj[2] = 0.0;
  /* general constants --------------------------------------------------- */


  /* physical constants -------------------------------------------------- */
  g_phc.c = 2.9979e10;
  g_phc.e = 4.8032e-10;

  g_phc.mp = 1.6726e-24;
  g_phc.me = 9.1094e-28; // original value

  g_phc.eV  = 1.6022e-12;
  g_phc.keV = 1.0e3*1.6022e-12;

  g_phc.den_min = 1.0e-5;
  g_phc.tem_min = 1.0e-5;

  /* set QUADPACK parameters for circular case */
  g_qdpk.epsabs = 1.0e-10;
  g_qdpk.epsrel = 1.0e-5;
  g_qdpk.irule = 5;

  /*g_phc.mi = 2.0;
    g_phc.mI = 12.0;
    g_phc.Zi = 1.0;
    g_phc.ZI = 6.0;*/
  /* physical constants -------------------------------------------------- */
}


/* reset flag0's of all global variables */
void reset_flag0_(void)
{
  int i, j, is;

  g_stime.flag0 = 0;
  g_ptl.flag0 = 0;
  g_phc.flag0 = 0;
  g_c.flag0 = 0;
  g_pall.flag0 = 0;
  g_sys.flag0 = 0;
  g_my_lim.flag0 = 0;
  g_eqd.flag0 = 0;
  g_RZ.flag0 = 0;
  g_cir.flag0 = 0;
  g_ase.flag0 = 0;
  g_eprof.flag0 = 0;
  g_xth.flag0 = 0;
  g_fld.flag0 = 0;
  g_tmp.flag0 = 0;
  g_damp.flag0 = 0;
  g_coll.flag0 = 0;
  g_dg.flag0 = 0;
  g_badk.flag0 = 0;
  g_dg_phs.flag0 = 0;
  g_diag_3df.flag0 = 0;
  g_zpk.flag0 = 0;
  g_ef.flag0 = 0;
  g_read_eprof.flag0 = 0;
  g_inpv.flag0 = 0;

  for(is = 0; is < NUM_SP_MAX; is++) 
  {
    for(j = 0; j < N_DGF_THE; j++) 
    {
      g_fls0[is][j].flag0 = 0;
    }
  }

  for(is = 0; is < NUM_SP_MAX; is++)
  {
    g_eprof.sp[is] = sp_non;
  }
}

void set_equilibrium_geometry_(void)
{
  switch(g_inpv.op_geo_info)
  {
    /* analytic concentric circular geometry */
    case 0: set_cir_efit_data();    
            break;
            /* general geometry given in GEQDSK format */
    case 1: read_geqdsk_and_test(); 
            break;
            /* analytic shaped geometry with specified central q-value */
    case 2: set_ase_efit_data();    
            break;
            /* analytic shaped geometry with q-profile (circular) */
    case 3: set_ase_efit_data();    
            break;

            /* otherwise default is set as circular */
    default: set_cir_efit_data();  
             break;
  }
}



/* skip characters until nn-'\n' are met. */
void skip_line(FILE *fp, int nn)
{
  int n=0;
  char ch='a';

  while(n < nn)
  {
    fscanf(fp, "%c", &ch);
    if(ch == '\n') n++;
  }
}

/* set global mpi-barrier with message */
void my_mpi_barrier(char *message)
{
  if(g_pall.ipe == 0)
  {
    printf("%s\n\n", message);
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/* inline.c */

inline int idx_1d_QS(int i)
{
        return i+1;
}

inline int idx_1d_CS(int i)
{
  return i+1;
}

/*inline double absv(double x)
{
	if(x >= 0.0)
		return x;
	else
		return -x;
}*/

double mod_th(double x)
{
  if(x >= 0.0 && x < g_c.twopi)
    return x;
  else
    if(x >= g_c.twopi)
      return x - ((int)(x*g_c.twopi_inv))*g_c.twopi;
    else
      return x - ((int)(x*g_c.twopi_inv)-1)*g_c.twopi;
}

/* 0th-spline */
inline double ngp(double x)
{
	if(x < 0.5)
		return 1.0;
	else
		return 0.0;
}

/* 1st-spline */
double cic(double x)
{
	if(x < 1.0)
		return 1.0 - x;
	else
		return 0.0;
}


/* 2nd-spline */
double tsc(double x)
{
	if(x < 0.5)
		return 0.75 - x*x;
	else
		if(x < 1.5)
			return 0.5*(1.5-x)*(1.5-x);
		else
			return 0.0;
}


inline double dtsc(double x)
{
	return cic(x+0.5) - cic(x-0.5);
}

/* 3rd spline */
inline double csc(double x)
{
	if(x < 2.0)
		if(x < 1.0)
			return (1.0+3.0*(1.0-x)*(1.0 + (1.0-x)*(1.0 - (1.0-x))))/6.0;
		else
			return (2.0-x)*(2.0-x)*(2.0-x)/6.0;
	else
		return 0.0;
}

double th_fn(double R, double Z)
{
  if(R == 1.0)
    if(Z == 0.0)
      return 0.0;
    else
      if(Z < 0.0)
        return 0.75*g_c.twopi;
      else
        return 0.25*g_c.twopi;
  else
    return double_atan2(Z,R-1.0);
}

inline int idx_thb(int n, int i)
{
  return i + (g_xth.nxx+1)*(n-g_xth.n1)/g_xth.dn;
}

inline int idxV(int i, int j)
{
  return j + g_xth.nth*i;
}

inline int jmod_th(int j)
{
  if(j < g_xth.nth)
    if(j >= 0)
      return j;
    else
      return j + g_xth.nth;
  else
    return j - g_xth.nth;
}

inline int idxnij(int n, int i, int j)
{
  return j + g_xth.nth*i + g_xth.nv*(n-g_xth.n1)/g_xth.dn;
}

inline complex fld_thb(dcomplex f[], int n, int i, int j)
{
  if(j < g_xth.nth)
    if(j >= 0)
      return f[idxnij(n,i,j)];
    else
      return f[idxnij(n,i,j+g_xth.nth)]*g_xth.thb[idx_thb(n,i)];
  else
    return f[idxnij(n,i,j-g_xth.nth)]/g_xth.thb[idx_thb(n,i)];
}

inline double mod_zt(double x)
{
  if(x >= 0.0 && x < g_xth.zt1)
    return x;
  else
    if(x >= g_xth.zt1)
      return x - ((int)(x*g_xth.zt1_inv))*g_xth.zt1;
    else
      return x - ((int)(x*g_xth.zt1_inv)-1)*g_xth.zt1;
}

inline int idx_dg_phsx(int i, int j)
{
  return j + N_DG_PHS_TH*i;
}

inline int idx_dg_phsv(int l, int p, int k)
{
  return k + N_DG_PHS_VAR*p + N_DG_PHS_VAR*(N_DG_PHS_VP+1)*l;
}

inline int idx_dgf(int l, int p, int nv)
{
  return p + nv*(l+nv-1);
}

inline int idxn_ij(int n, int ij)
{
  return ij + g_xth.nv*(n - g_xth.n1)/g_xth.dn;
}


