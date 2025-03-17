#include "cnst-srt.h"

#define idx_1d_CS(i) ((i)+1)
#define idx_2d_CS(i,j,nx,ny) ((j)+1 + (ny+3)*((i)+1))
#define idx_RZ_CS(i,j) ((j)+1 + (g_RZ.NZ+3)*((i)+1))


void set_cir_efit_data(void)
{
  int i, j;
  double R, Z, rr, ma0, mR0, B0, rrl, rrm, pv[6], psi;


  /* take inputs ------------------------------------------------- */
  ma0 = g_inpv.cir_a0;
  mR0 = g_inpv.cir_R0;
  B0  = g_inpv.cir_B0;

  g_cir.q[0] = g_inpv.cir_q0;
  g_cir.q[1] = g_inpv.cir_q1;
  g_cir.q[2] = g_inpv.cir_q2;
  g_cir.q[3] = g_inpv.cir_q3;
  g_cir.q[4] = g_inpv.cir_q4;

  /*// Cyclon base case 
    ma0 = 0.48; mR0 = 1.3; B0  = 1.9;
    g_cir.q[0] = 0.84; 
    g_cir.q[1] = 0.0; 
    g_cir.q[2] = 2.18*mR0*mR0/(ma0*ma0);
    g_cir.q[3] = 0.0;
    g_cir.q[4] = 0.0;
  // Waltz standard case 
  ma0 = 0.5; mR0 = 1.5; B0  = 2.0;
  g_cir.q[0] = 1.0; 
  g_cir.q[1] = 0.0; 
  g_cir.q[2] = 3.0*3.0/(0.5*0.5);
  g_cir.q[3] = 0.0;
  g_cir.q[4] = 0.0;
  // Waltz standard case2 (linear q) 
  ma0 = 0.5; mR0 = 1.5; B0  = 2.0;
  g_cir.q[0] = 0.0; 
  g_cir.q[1] = 0.0; 
  g_cir.q[2] = 12.0;
  g_cir.q[3] = 0.0;
  g_cir.q[4] = 0.0;
  // TFTR L-mode case 
  ma0 = 0.94; mR0 = 2.7; B0  = 3.1;
  g_cir.q[0] = 1.117; 
  g_cir.q[1] = 0.0;
  g_cir.q[2] = 36.0705;
  g_cir.q[3] = 0.0;
  g_cir.q[4] = 0.0;*/
  /* take inputs ------------------------------------------------- */


  /* circular parameters in cgs units */
  g_cir.a0 = 1.0e2*ma0;
  g_cir.R0 = 1.0e2*mR0;
  g_cir.B0 = 1.0e4*B0;


  /* ------------------------------------------------------------- */
  g_eqd.NR = 129-1;
  g_eqd.NZ = 129-1;
  g_eqd.Npsi = 129-1;

  /* this parameter indicates the position of circular limiter
     in the unit of the minor radius */
  rrm = rrl = g_inpv.rrl_cir;
  //rrm = rrl*1.2;

  rrl *= ma0;
  rrm *= ma0;

  g_eqd.Rmin = mR0 - rrm;
  g_eqd.Rmax = mR0 + rrm;
  g_eqd.Zmin = -rrm;
  g_eqd.Zmax =  rrm;
  g_eqd.psi0 = 0.0;
  g_eqd.psi1 = cir_psi_fn_ep(ma0/mR0)*mR0*mR0*B0;
  g_eqd.dR = (g_eqd.Rmax - g_eqd.Rmin)/g_eqd.NR;
  g_eqd.dZ = (g_eqd.Zmax - g_eqd.Zmin)/g_eqd.NZ;
  g_eqd.dpsi = (g_eqd.psi1 - g_eqd.psi0)/g_eqd.Npsi;

  /* allocate variables */
  g_eqd.I_CS = (double*)malloc(sizeof(double)*(g_eqd.Npsi+3));
  g_eqd.psi_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));
  g_eqd.I_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));

  /* initialize interpolation */
  for(i = 0; i <= g_eqd.Npsi; i++) g_eqd.I_CS[idx_1d_CS(i)] = mR0*B0; 
  init_intp_1d_CS_umfpack(g_eqd.psi0, g_eqd.psi1, g_eqd.Npsi, g_eqd.I_CS);

  for(i = 0; i <= g_eqd.NR; i++) for(j = 0; j <= g_eqd.NZ; j++)
  {
    R = g_eqd.Rmin + i*g_eqd.dR;
    Z = g_eqd.Zmin + j*g_eqd.dZ;
    rr = (R-mR0)*(R-mR0) + Z*Z;
    rr /= mR0*mR0;

    psi = cir_psi_fn_ep(double_sqrt(rr));
    g_eqd.psi_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] = psi*mR0*mR0*B0;
  }
  init_intp_2d_CS_umfpack(g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS);


  /* set plasma center */
  g_eqd.a0 = ma0;
  g_eqd.Rc = mR0;
  g_eqd.Zc = 0.0;
  intp_dfn_2d_CS(g_eqd.Rc, g_eqd.Zc,
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);
  g_eqd.psic = pv[0];
  g_eqd.Ic = intp_1d_CS(pv[0],g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  g_eqd.Bc = double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+g_eqd.Ic*g_eqd.Ic)/g_eqd.Rc;

  /* initialize global variables describing edge geometry */
  g_eqd.nbbbs = 72;
  g_eqd.rbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  g_eqd.zbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    g_eqd.rbbbs[j] = ma0*double_cos(j*g_c.twopi/g_eqd.nbbbs) + mR0;
    g_eqd.zbbbs[j] = ma0*double_sin(j*g_c.twopi/g_eqd.nbbbs);
  }

  /* simple rectangular limiter */
  g_eqd.limitr = 4;
  g_eqd.rlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_eqd.zlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_eqd.rlim[0] = g_eqd.Rmax; g_eqd.zlim[0] = g_eqd.Zmax;
  g_eqd.rlim[1] = g_eqd.Rmin; g_eqd.zlim[1] = g_eqd.Zmax;
  g_eqd.rlim[2] = g_eqd.Rmin; g_eqd.zlim[2] = g_eqd.Zmin;
  g_eqd.rlim[3] = g_eqd.Rmax; g_eqd.zlim[3] = g_eqd.Zmin;

  /* simple circular limiter 
     g_eqd.limitr = 32;
     g_eqd.rlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
     g_eqd.zlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
     for(j = 0; j < g_eqd.limitr; j++)
     {
     g_eqd.rlim[j] = rrl*double_cos(j*g_c.twopi/g_eqd.limitr) + mR0;
     g_eqd.zlim[j] = rrl*double_sin(j*g_c.twopi/g_eqd.limitr);
     }*/

  /* set x-point */
  assert(N_NULLX == 1);
  g_eqd.Rx[0] = mR0;
  g_eqd.Zx[0] = -ma0;
  g_eqd.psix[0] = g_eqd.psi1;


  /* initialize I_2d_CS interpolation. */
  for(i = 0; i <= g_eqd.NR; i++) for(j = 0; j <= g_eqd.NZ; j++)
  {
    R = g_eqd.Rmin + i*g_eqd.dR;
    Z = g_eqd.Zmin + j*g_eqd.dZ;

    psi = intp_2d_CS(R, Z, g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
        g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ, g_eqd.psi_2d_CS);

    if(psi <= g_eqd.psi1 && psi >= g_eqd.psi0 && Z > g_eqd.Zx[0])
      g_eqd.I_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] =
        intp_1d_CS(psi,g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
    else
      g_eqd.I_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] =
        intp_1d_CS(g_eqd.psi1,g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  }
  init_intp_2d_CS_umfpack(g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.I_2d_CS);

  my_geqdsk_to_cgs();
  /* ------------------------------------------------------------- */

  //my_mpi_barrier("set_cir_efit_data()");
}


/* return numerical interpolated q-profile in circular geometry */
/* ---------------------------------------------------------------------- */
double cir_q_afn_ep(double ep)
{
  double val;

  val = g_cir.q[0] + g_cir.q[1]*ep + g_cir.q[2]*ep*ep 
    + g_cir.q[3]*ep*ep*ep + g_cir.q[4]*ep*ep*ep*ep;
  return val;
}
double cir_dpde_afn_ep(double ep)
{
  double q;

  if(ep >= 1.0)
  {
    printf("ep in cir_dpde_afn exceeds 1.0 : %e\n", ep);
    exit(1);
  }

  q = cir_q_afn_ep(ep);

  if(q <= 0.001)
  {
    printf("q-value in cir_dpde_afn is too small : %e %e\n", q,ep);
    exit(1);
  }

  return ep/(q*double_sqrt(1.0-ep*ep));
}
/* ---------------------------------------------------------------------- */


double cir_psi_fn_ep(double ep)
{
  double psi;

  psi = my_G_K(cir_dpde_afn_ep, 0.0, ep);
  return max(psi, 0.0);

  /*// Cyclon base case 
    psi = 0.132837 - 0.0609567*atanh(0.974726*double_sqrt(1.0-ep*ep));
  // Waltz standard case 
  psi = 0.0682743 - 0.0273998*atanh(0.986394*double_sqrt(1.0-ep*ep));
  // Waltz standard case2 (linear q) 
  psi = 0.0833333*asin(ep); 
  // TFTR L-mode case 
  psi = 0.0665729 - 0.0273039*atanh(0.984867*double_sqrt(1.0-ep*ep));*/
}

/* START of geom-ase.c */
/* set analytic shaped equilibrium */
void set_ase_efit_data(void)
{
  int i, j;
  double R, Z, ep, A, B, cs, ss, ma0, mR0, B0, rrl, rrm, pv[6], psi, rbx, zbx;
  char fname[CHARLEN]; 
  FILE *fp;

  /* -------------------------------------------------------------------- */
  /* read and set shaping parameters from input */
  g_ase.R0     = g_inpv.ase_R0;   // R at magnetic axis
  g_ase.B0     = g_inpv.ase_B0;   // B at magnetic axis
  g_ase.a0     = g_inpv.ase_a0;   // minor radius on outer mid-plane
  g_ase.q0     = g_inpv.ase_q0;   // q at magnetic axis

  g_ase.kappa  = g_inpv.ase_kappa;  // elongation
  g_ase.del    = g_inpv.ase_del;    // triangularity
  g_ase.shift  = g_inpv.ase_shift;  // Shafranov shift
  assert(g_ase.shift < 1.0);

  mR0 = g_ase.R0;
  ma0 = g_ase.a0;
  B0  = g_ase.B0;

  /* convert to cgs unit */
  g_ase.R0 *= 1.0e2;
  g_ase.a0 *= 1.0e2;
  g_ase.B0 *= 1.0e4;

  /* q-profile parameters for concentric circular case */
  g_cir.q[0] = g_inpv.cir_q0; 
  g_cir.q[1] = g_inpv.cir_q1;
  g_cir.q[2] = g_inpv.cir_q2;
  g_cir.q[3] = g_inpv.cir_q3;
  g_cir.q[4] = g_inpv.cir_q4;
  /* -------------------------------------------------------------------- */


  /* -------------------------------------------------------------------- */
  /* set equilibrium interpolation resolution */
  g_eqd.NR = 129-1;
  g_eqd.NZ = 129-1;
  g_eqd.Npsi = 129-1;
  /* -------------------------------------------------------------------- */


  /* -------------------------------------------------------------------- */
  /* this parameter indicates the position of circular limiter
     in the unit of the minor radius */
  rrm = rrl = g_inpv.rrl_cir; //rrm = rrl*1.2;

  rrl *= ma0;
  rrm *= ma0;

  g_eqd.Rmin = mR0*(1.0 - g_ase.shift) - rrm;
  g_eqd.Rmax = mR0*(1.0 - g_ase.shift) + rrm;
  g_eqd.Zmin = -g_ase.kappa*rrm;
  g_eqd.Zmax =  g_ase.kappa*rrm;
  g_eqd.psi0 = 0.0;
  g_eqd.psi1 = ase_psi_fn_RZ(1.0 - g_ase.shift + ma0/mR0, 0.0)*mR0*mR0*B0;
  g_eqd.dR = (g_eqd.Rmax - g_eqd.Rmin)/g_eqd.NR;
  g_eqd.dZ = (g_eqd.Zmax - g_eqd.Zmin)/g_eqd.NZ;
  g_eqd.dpsi = (g_eqd.psi1 - g_eqd.psi0)/g_eqd.Npsi;

  /* allocate variables */
  g_eqd.I_CS = (double*)malloc(sizeof(double)*(g_eqd.Npsi+3));
  g_eqd.psi_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));
  g_eqd.I_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));

  /* initialize interpolation */
  for(i = 0; i <= g_eqd.Npsi; i++) g_eqd.I_CS[idx_1d_CS(i)] = mR0*B0; 
  init_intp_1d_CS_umfpack(g_eqd.psi0, g_eqd.psi1, g_eqd.Npsi, g_eqd.I_CS);

  if(g_pall.ipe == 0) 
  {
    sprintf(fname, "test-ase1.txt", g_sys.work_dir);
    fp = fopen(fname, "w"); assert(fp != NULL);
  }

  for(i = 0; i <= g_eqd.NR; i++) 
  {
    R = g_eqd.Rmin + i*g_eqd.dR;

    for(j = 0; j <= g_eqd.NZ; j++)
    {
      Z = g_eqd.Zmin + j*g_eqd.dZ;

      psi = ase_psi_fn_RZ(R/mR0, Z/mR0);
      g_eqd.psi_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] = psi*mR0*mR0*B0;

      if(g_pall.ipe == 0) 
        fprintf(fp, "%d %d %e %e %e\n", i, j, R/mR0, Z/mR0, psi);
    }
    if(g_pall.ipe == 0) fprintf(fp, "\n");
  }
  init_intp_2d_CS_umfpack(g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS);

  if(g_pall.ipe == 0) fclose(fp);
  /* -------------------------------------------------------------------- */


  /* -------------------------------------------------------------------- */
  /* set plasma center */
  g_eqd.a0 = ma0;
  g_eqd.Rc = mR0*(1.0 - g_ase.shift);
  g_eqd.Zc = 0.0;
  intp_dfn_2d_CS(g_eqd.Rc, g_eqd.Zc,
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);
  g_eqd.psic = pv[0];
  g_eqd.Ic = intp_1d_CS(pv[0],g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  g_eqd.Bc = double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+g_eqd.Ic*g_eqd.Ic)/g_eqd.Rc;
  /* -------------------------------------------------------------------- */


  /* ---------------------------------------------------------------------- */
  /* initialize and set global variables describing edge geometry */
  g_eqd.nbbbs = 72;
  g_eqd.rbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  g_eqd.zbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);

  ep = ma0/mR0; 
  A = g_ase.shift/(ep*ep);
  B = g_ase.shift/(ep*ep) + g_ase.del/ep;
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    cs = double_cos(j*g_c.twopi/g_eqd.nbbbs);
    ss = double_sin(j*g_c.twopi/g_eqd.nbbbs);

    g_eqd.rbbbs[j] = mR0*(1.0 + ep*cs - A*ep*ep*cs*cs - B*ep*ep*ss*ss);
    g_eqd.zbbbs[j] = mR0*g_ase.kappa*ep*ss;

    /*g_eqd.rbbbs[j] = mR0 + 1.2*(g_eqd.rbbbs[j] - mR0);
      g_eqd.zbbbs[j] = 1.2*g_eqd.zbbbs[j];*/
  }

  /* set x-point */
  rbx =  mR0*(1.0 - B*ep*ep);
  zbx = -mR0*g_ase.kappa*ep;
  assert(N_NULLX == 1);
  g_eqd.Rx[0] = rbx;
  g_eqd.Zx[0] = zbx;
  g_eqd.psix[0] = g_eqd.psi1;
  /* ---------------------------------------------------------------------- */


  /* ---------------------------------------------------------------------- */
  /* simple rectangular limiter */
  g_eqd.limitr = 4;
  g_eqd.rlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_eqd.zlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_eqd.rlim[0] = g_eqd.Rmax; g_eqd.zlim[0] = g_eqd.Zmax;
  g_eqd.rlim[1] = g_eqd.Rmin; g_eqd.zlim[1] = g_eqd.Zmax;
  g_eqd.rlim[2] = g_eqd.Rmin; g_eqd.zlim[2] = g_eqd.Zmin;
  g_eqd.rlim[3] = g_eqd.Rmax; g_eqd.zlim[3] = g_eqd.Zmin;

  /* simple circular limiter 
     g_eqd.limitr = 32;
     g_eqd.rlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
     g_eqd.zlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
     for(j = 0; j < g_eqd.limitr; j++)
     {
     g_eqd.rlim[j] = rrl*double_cos(j*g_c.twopi/g_eqd.limitr) + mR0;
     g_eqd.zlim[j] = rrl*double_sin(j*g_c.twopi/g_eqd.limitr);
     }*/
  /* ---------------------------------------------------------------------- */


  /* ---------------------------------------------------------------------- */
  /* initialize I_2d_CS interpolation. */
  for(i = 0; i <= g_eqd.NR; i++) for(j = 0; j <= g_eqd.NZ; j++)
  {
    R = g_eqd.Rmin + i*g_eqd.dR;
    Z = g_eqd.Zmin + j*g_eqd.dZ;

    psi = intp_2d_CS(R, Z, g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
        g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ, g_eqd.psi_2d_CS);

    if(psi <= g_eqd.psi1 && psi >= g_eqd.psi0 && Z > g_eqd.Zx[0])
      g_eqd.I_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] =
        intp_1d_CS(psi,g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
    else
      g_eqd.I_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] =
        intp_1d_CS(g_eqd.psi1,g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  }
  init_intp_2d_CS_umfpack(g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.I_2d_CS);

  my_geqdsk_to_cgs();
  /* ---------------------------------------------------------------------- */

  //my_mpi_barrier("set_ase_efit_data()");
}

/* return psi for analytic shaped equilibrium normalized by B0*R0*R0 */
double ase_psi_fn_RZ(double R, double Z)
{
  int ll, ll_max = 100;
  double psi0, tau, E, ep, Rx2, fc, R2, R22, dR, epp, err;

  ep = double_sqrt((R-1.0+g_ase.shift)*(R-1.0+g_ase.shift)
      + Z*Z/(g_ase.kappa*g_ase.kappa));

  if(ep > 1.0e-7)
  {
    for(ll = 0; ll < ll_max; ll++)
    {
      dR = R + g_ase.shift - 1.0 + g_ase.del*Z*Z/(ep*g_ase.kappa*g_ase.kappa);
      epp = double_sqrt(dR*dR + Z*Z/(g_ase.kappa*g_ase.kappa));

      err = epp-ep; err = absv(err)/ep;
      if(err < 1.0e-7) break;
      else ep = epp;
    }
  }

  if(g_inpv.op_geo_info == 2)
  {
    psi0 = g_ase.kappa*(1.0 + 2.0*g_ase.shift)/(8.0*g_ase.q0);
    return 4.0*psi0*ep*ep;
  }
  else if(g_inpv.op_geo_info == 3)
  {
    psi0 = g_ase.kappa*(1.0 + 2.0*g_ase.shift)/8.0;
    return 8.0*psi0*cir_psi_fn_ep(ep);
  }
  else
  {
    printf("Unknown equilibrium option ???\n");
    exit(1);
  }
}

