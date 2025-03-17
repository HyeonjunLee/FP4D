#include "cnst-srt.h"


#define idx_RZ_CS(i,j) ((j)+1 + (g_RZ.NZ+3)*((i)+1))
/* 
   Set RZ simulation box by normalize and shift geqdsk information.
   */
void set_rz_simulation_domain_(void)
{
  int i, j, ij;
  double R, Z, pv[6], Ival, eqd_R, eqd_Z, psi_val;
  double R_inv, R2_inv, R2B_inv, Iv[6], Bv[3], B, B_inv, dBdR, dBdZ;
  double th, ep, si_th, cs_th, dtdR, dtdZ;

  /* number of R, Z grid points */
  g_RZ.NR = g_eqd.NR;
  g_RZ.NZ = g_eqd.NZ;

  /* normalized & shifted RZ domain box info 
     In this box, plasma center is moved to (1.0, 0.0) 
     R/Rc, (Z - Zc)/Rc 
     */
  g_RZ.Rmax = g_eqd.Rmax/g_eqd.Rc;
  g_RZ.Rmin = g_eqd.Rmin/g_eqd.Rc;
  g_RZ.Zmax = (g_eqd.Zmax - g_eqd.Zc)/g_eqd.Rc;
  g_RZ.Zmin = (g_eqd.Zmin - g_eqd.Zc)/g_eqd.Rc;
  g_RZ.dR = (g_RZ.Rmax - g_RZ.Rmin)/g_RZ.NR;
  g_RZ.dZ = (g_RZ.Zmax - g_RZ.Zmin)/g_RZ.NZ;
  g_RZ.dR_inv = 1.0/g_RZ.dR;
  g_RZ.dZ_inv = 1.0/g_RZ.dZ;

  /* normalized & shifted x-points */
  for(i = 0; i < N_NULLX; i++)
  {
    g_RZ.Rx[i] = g_eqd.Rx[i]/g_eqd.Rc;
    g_RZ.Zx[i] = (g_eqd.Zx[i] - g_eqd.Zc)/g_eqd.Rc;
  }

  /* normalized & shifted plasma boundary rectangle */
  g_RZ.Rbmax = g_eqd.rbbbs_max/g_eqd.Rc; 
  g_RZ.Rbmin = g_eqd.rbbbs_min/g_eqd.Rc;
  g_RZ.Zbmax = (g_eqd.zbbbs_max - g_eqd.Zc)/g_eqd.Rc;
  g_RZ.Zbmin = (g_eqd.zbbbs_min - g_eqd.Zc)/g_eqd.Rc; 
  /* allocation */
  g_RZ.psi_CS   = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.xx_CS    = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.I_CS     = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.B_CS     = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.bdB_CS   = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.J_B_CS   = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));
  g_RZ.I_R2B_CS = (double*)malloc(sizeof(double)*(g_RZ.NR+3)*(g_RZ.NZ+3));

  /* normalized psi value at the last closed surface */
  g_RZ.psi_last = (g_eqd.psi1-g_eqd.psi0)/(g_eqd.Bc*g_eqd.Rc*g_eqd.Rc);
  /* set normalized and shifted psi & I array */
  for(i = 0; i <= g_RZ.NR; i++)
  {
    R = g_RZ.Rmin + i*g_RZ.dR;
    eqd_R = R*g_eqd.Rc;

    for(j = 0; j <= g_RZ.NZ; j++)
    {
      Z = g_RZ.Zmin + j*g_RZ.dZ;
      eqd_Z = g_eqd.Zc + Z*g_eqd.Rc;

      ij = idx_RZ_CS(i,j);

      intp_dfn_2d_CS(eqd_R, eqd_Z, 
          g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
          g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
          g_eqd.psi_2d_CS, pv);

      if(pv[0] <= g_eqd.psi1 && pv[0] >= g_eqd.psi0 && Z > g_RZ.Zx[0]) 
        Ival = intp_1d_CS(pv[0],g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,
            g_eqd.I_CS);
      else
        Ival = intp_1d_CS(g_eqd.psi1,g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,
            g_eqd.I_CS);

      /* (psi - psic)/(Bc*Rc*Rc) */
      g_RZ.psi_CS[ij] = (pv[0] - g_eqd.psic)/(g_eqd.Bc*g_eqd.Rc*g_eqd.Rc);

      psi_val = g_RZ.psi_CS[ij];
      psi_val = max(0.0, psi_val);
      g_RZ.xx_CS[ij] = double_sqrt(psi_val/g_RZ.psi_last);

      /* I/(Rc*Bc) */
      g_RZ.I_CS[ij] = Ival/(g_eqd.Rc*g_eqd.Bc);

      /* B/Bc */
      g_RZ.B_CS[ij] = 
        double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+Ival*Ival)/(eqd_R*g_eqd.Bc);
    }
  }

  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.psi_CS);
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.xx_CS);
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.I_CS);
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.B_CS);

  /* commonly used geometric quantities */
  for(i = 0; i <= g_RZ.NR; i++)
  {
    R = g_RZ.Rmin + i*g_RZ.dR;
    R_inv = 1.0/R;
    R2_inv = R_inv*R_inv;

    for(j = 0; j <= g_RZ.NZ; j++)
    {
      Z = g_RZ.Zmin + j*g_RZ.dZ;

      ij = idx_RZ_CS(i,j);

      intp_dpsi_RZ(R, Z, pv);
      intp_dI_RZ(R, Z, Iv);

      Bv[0] =  BP_SIGN*pv[2]*R_inv;
      Bv[1] = -BP_SIGN*pv[1]*R_inv;
      Bv[2] =  BT_SIGN*Iv[0]*R2_inv;

      B = double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+Iv[0]*Iv[0])*R_inv;
      B_inv = 1.0/B;
      R2B_inv = R2_inv*B_inv;

      dBdR = R2B_inv*(pv[1]*pv[3]+pv[2]*pv[5]+Iv[0]*Iv[1]) - B*R_inv;
      dBdZ = R2B_inv*(pv[2]*pv[4]+pv[1]*pv[5]+Iv[0]*Iv[2]);

      /* b dot grad B */
      g_RZ.bdB_CS[idx_RZ_CS(i,j)] = (Bv[0]*dBdR + Bv[1]*dBdZ)*B_inv;

      th = th_fn(R, Z);
      ep = double_sqrt((R-1.0)*(R-1.0) + Z*Z);
      si_th = double_sin(th); dtdR = -si_th/(ep+1.0e-30);
      cs_th = double_cos(th); dtdZ =  cs_th/(ep+1.0e-30);

      /* J/B */
      g_RZ.J_B_CS[idx_RZ_CS(i,j)] = (pv[1]*dtdZ - pv[2]*dtdR)/(R*B);
      /* I/(R*R*B) */
      g_RZ.I_R2B_CS[idx_RZ_CS(i,j)] = Iv[0]/(R*R*B);
    }
  }
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.bdB_CS);
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.J_B_CS);
  init_intp_2d_CS_umfpack(g_RZ.Rmin, g_RZ.Rmax, g_RZ.NR, 
      g_RZ.Zmin, g_RZ.Zmax, g_RZ.NZ, g_RZ.I_R2B_CS);


  //my_mpi_barrier("set_RZ_simulation_domain()");

  set_limiter();
  simplify_RZ_curve(g_RZ.lim_R, g_RZ.lim_Z, g_eqd.limitr, 
      0.001, &g_RZ.my_limitr);
  if(g_pall.ipe == 0) 
  { 
    printf("g_my_limitr = %d (%d)\n", g_RZ.my_limitr, g_eqd.limitr); 
    fflush(stdout);
  }
  init_my_limiter();
}


void set_limiter(void)
{
  int j;

  g_RZ.lim_R = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_RZ.lim_Z = (double*)malloc(sizeof(double)*g_eqd.limitr);

  /* set limiter information (R, Z) : normalized & shifted */
  for(j = 0; j < g_eqd.limitr; j++)
  {
    g_RZ.lim_R[j] = g_eqd.rlim[j]/g_eqd.Rc;
    g_RZ.lim_Z[j] = (g_eqd.zlim[j] - g_eqd.Zc)/g_eqd.Rc;
  }

  /* find bounding rectangle of the limiter curve */
  g_RZ.lim_Rmax = g_RZ.lim_Rmin = g_RZ.lim_R[0];
  g_RZ.lim_Zmax = g_RZ.lim_Zmin = g_RZ.lim_Z[0];
  for(j = 1; j < g_eqd.limitr; j++)
  {
    g_RZ.lim_Rmax = max(g_RZ.lim_Rmax, g_RZ.lim_R[j]);
    g_RZ.lim_Rmin = min(g_RZ.lim_Rmin, g_RZ.lim_R[j]);

    g_RZ.lim_Zmax = max(g_RZ.lim_Zmax, g_RZ.lim_Z[j]);
    g_RZ.lim_Zmin = min(g_RZ.lim_Zmin, g_RZ.lim_Z[j]);
  }
}


/* limiter.c */
void simplify_RZ_curve(double *cv_R, double *cv_Z, int num, double err, int *rnum)
{
  int j;
  double *tmpR, *tmpZ, dR, dZ, len;

  tmpR = (double*)malloc(sizeof(double)*num);
  tmpZ = (double*)malloc(sizeof(double)*num);

  /* remove almost identical point */
  *rnum = 0;
  for(j = 0; j < num; j++)
  {
    if(j == 0)
    {
      tmpR[*rnum] = cv_R[j];
      tmpZ[*rnum] = cv_Z[j];
      (*rnum)++;
    }
    else
    {
      dR = cv_R[j] - tmpR[*rnum-1];
      dZ = cv_Z[j] - tmpZ[*rnum-1];
      len = double_sqrt(dR*dR + dZ*dZ);

      if(len > err)
      {
        tmpR[*rnum] = cv_R[j];
        tmpZ[*rnum] = cv_Z[j];
        (*rnum)++;
      }
    }
  }

  for(j = 0; j < *rnum; j++)
  {
    cv_R[j] = tmpR[j];
    cv_Z[j] = tmpZ[j];
  }

  free(tmpR);
  free(tmpZ);
}


void init_my_limiter(void)
{
  int j, j1, j2;
  double th1, th2, Rt, Zt, th;
  FILE *fp;
  char fname[CHARLEN];

  g_my_lim.seg = (line_seg_t*)malloc(sizeof(line_seg_t)*g_RZ.my_limitr);
  g_my_lim.idx = (int*)malloc(sizeof(int)*g_RZ.my_limitr);

  for(j = 0; j < g_RZ.my_limitr; j++)
  {
    j1 = j;
    j2 = j+1;
    if(j2 == g_RZ.my_limitr) j2 = 0;

    th1 = th_fn(g_RZ.lim_R[j1],g_RZ.lim_Z[j1]);
    th2 = th_fn(g_RZ.lim_R[j2],g_RZ.lim_Z[j2]);
    th1 = mod_th(th1);
    th2 = mod_th(th2);

    g_my_lim.seg[j].R1 = g_RZ.lim_R[j1]; 
    g_my_lim.seg[j].Z1 = g_RZ.lim_Z[j1];
    g_my_lim.seg[j].th1 = th1;

    g_my_lim.seg[j].R2 = g_RZ.lim_R[j2];
    g_my_lim.seg[j].Z2 = g_RZ.lim_Z[j2];
    g_my_lim.seg[j].th2 = th2;

    /* determine the orientation of the limiter segment 
       bugfix Jan.30.2008 : in determining the segment crossing outside
       mid-plane */
    if(((g_my_lim.seg[j].Z1)*(g_my_lim.seg[j].Z2) <= 0.0) && 
        (g_my_lim.seg[j].R1 > 1.0)) 
    {
      if(g_my_lim.seg[j].Z1 < 0.0) g_my_lim.seg[j].p = -1;
      g_my_lim.seg[j].p = -1;
    }
    else 
      g_my_lim.seg[j].p =  1;

    if(th2 < th1)
      /*if((g_my_lim.seg[j].p ==  1 && th2 < th1) ||
        (g_my_lim.seg[j].p == -1 && th2 > th1))*/
    {
      Rt = g_my_lim.seg[j].R1;
      Zt = g_my_lim.seg[j].Z1;
      th = g_my_lim.seg[j].th1;

      g_my_lim.seg[j].R1 = g_my_lim.seg[j].R2;
      g_my_lim.seg[j].Z1 = g_my_lim.seg[j].Z2;
      g_my_lim.seg[j].th1 = g_my_lim.seg[j].th2;

      g_my_lim.seg[j].R2 = Rt;
      g_my_lim.seg[j].Z2 = Zt;
      g_my_lim.seg[j].th2 = th;
    }
  }

  for(j = 0; j < g_RZ.my_limitr; j++)
  {
    assert(g_my_lim.seg[j].th1 <= g_my_lim.seg[j].th2);
    /*if(g_my_lim.seg[j].p == 1) 
      assert(g_my_lim.seg[j].th1 <= g_my_lim.seg[j].th2);
      else 
      assert(g_my_lim.seg[j].th1 >= g_my_lim.seg[j].th2);*/
  }

  if(g_pall.ipe == 0 && g_inpv.op_test_output == 1)
  {
    sprintf(fname, "test-lim1.txt", g_sys.work_dir);
    fp = fopen(fname, "w");
    assert(fp != NULL);

    for(j = 0; j < g_RZ.my_limitr; j++)
      fprintf(fp, "%e %e %d\n", g_RZ.lim_R[j], g_RZ.lim_Z[j], 
          g_my_lim.seg[j].p);

    fclose(fp);
  }
}

