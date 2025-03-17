#include "cnst-srt.h"

#define idx_S_CS(i,j)    ((j) + (N_FTH+1)*(i))
#define idx_xfth_CS(i,j) ((j)+1 + (N_FTH+3)*((i)+1))
#define idx_xth_CS(i,j,nx,nt) ((j) + (nt)*((i)+1))

void spline(double *x, double *y, int n, const double yp1, const double ypn, double *y2);
void splint(double *xa, double *ya, double *y2a, int n, const double x, double *y);

/* set basic information for uniform (g_xth.nxx,g_xth.nth) grid */
void setup_xx_theta_grid_(void)
{
  int is, i, j, ij;
  double ep, aa, R, Z;
  double x, th, dtdR, dtdZ, pv[6], Iv[6], B;

  /* radial direction -------------------------------------------------- */
  g_xth.nxx = g_inpv.nxx;
  g_xth.nxx1 = g_xth.nxx + 1;
  g_xth.nxx2 = g_xth.nxx + 2;
  g_xth.nxx3 = g_xth.nxx + 3;

  g_xth.xx0 = g_inpv.xx0;
  g_xth.xx0 = max(0.0, g_xth.xx0);

  g_xth.xx1 = g_inpv.xx1;
  g_xth.xx1 = min(g_xth.xx1, 1.0);

  g_xth.dxx = (g_xth.xx1-g_xth.xx0)/g_xth.nxx;
  g_xth.dxx_inv = 1.0/g_xth.dxx;
  g_xth.dxx2_inv = g_xth.dxx_inv*g_xth.dxx_inv;

  g_xth.xx = (double*)malloc(sizeof(double)*(g_xth.nxx+1));
  g_xth.xx[0] = g_xth.xx0;
  for(i = 1; i < g_xth.nxx; i++) g_xth.xx[i] = g_xth.xx0 + i*g_xth.dxx;
  g_xth.xx[g_xth.nxx] = g_xth.xx1;
  /* radial direction -------------------------------------------------- */


  /* poloidal direction ------------------------------------------------ */
  g_xth.nth = g_inpv.nth;

  g_xth.dth = g_c.twopi/g_xth.nth;
  g_xth.dth_inv = 1.0/g_xth.dth;
  g_xth.dth2_inv = g_xth.dth_inv*g_xth.dth_inv;
  /* poloidal direction ------------------------------------------------ */


  // number of grids on a poloidal plain
  g_xth.nv = g_xth.nxx1*g_xth.nth;

  /* uniform N_FTH grid */
  g_xth.dfth = g_c.twopi/N_FTH;
  g_xth.dfth_inv = 1.0/g_xth.dfth;
  g_xth.dfth2_inv = g_xth.dfth_inv*g_xth.dfth_inv;

  /* lookup and set plasma boundary (last flux surface) */
  set_plasma_boundary();

  /* setup epm(x), R(x,th), Z(x,th) interpolation */
  set_epm_R_Z_xfth();
  /* useful radial coordinat (for transp) : trho = sqrt(tflux/tflux_edge) */
  //setup_trho_grid();

  /* set boundary for field */
  g_inpv.ifd0 = max(1, g_inpv.ifd0);
  g_inpv.ifd1 = min(g_inpv.ifd1, g_xth.nxx-1);
  g_xth.xx_fd0 = g_xth.xx0 + g_inpv.ifd0*g_xth.dxx;
  g_xth.xx_fd1 = g_xth.xx0 + g_inpv.ifd1*g_xth.dxx;
  for(is = 0; is < NUM_SP_MAX; is++)
  {
    /* set boundary for plasma */
    g_inpv.ipd0[is] = g_inpv.ifd0;
    g_inpv.ipd1[is] = g_inpv.ifd1;

    g_xth.xx_pd0[is] = g_xth.xx0 + g_inpv.ipd0[is]*g_xth.dxx;
    g_xth.xx_pd1[is] = g_xth.xx0 + g_inpv.ipd1[is]*g_xth.dxx;
  }

  /* core grid : use theta-symmetric spline */
  g_xth.tss_i = 1;
  //g_xth.tss_i = 5;
  g_xth.tss_i = max(1, g_xth.tss_i);

  /* outter mid-plane r/R */
  g_xth.ep_min = epm_fn(g_xth.xx_fd0);
  g_xth.ep_max = epm_fn(g_xth.xx_fd1);

  /* mid-radious */
  g_xth.epc = 0.5*(g_xth.ep_min + g_xth.ep_max);
  intp_dpsi_RZ(1.0 + g_xth.epc, 0.0, pv);
  g_xth.xxc = double_sqrt(pv[0]/g_RZ.psi_last);

  /* XTH-particle domain */
  g_xth.xx_out0 = g_xth.xx_pd0[0];
  g_xth.xx_out1 = g_xth.xx_pd1[0];

  // radial domain to check conservation
  intp_dpsi_RZ(1.0 + 0.08, 0.0, pv); 
  g_xth.xx_cchk0 = double_sqrt(pv[0]/g_RZ.psi_last);
  intp_dpsi_RZ(1.0 + 0.28, 0.0, pv); 
  g_xth.xx_cchk1 = double_sqrt(pv[0]/g_RZ.psi_last);

  /* total simulation volume inside last flux surface */
  /* plasma elongation */
  intp_RZ_xfth(g_xth.xx1, 0.0, &aa, &Z);
  intp_RZ_xfth(g_xth.xx1, g_c.pi, &R, &Z);
  aa -= R;
  aa *= 0.5;
  aa *= aa;

  /* set geomeric quantities for later use -------------------------------- */
  g_xth.Bmod_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*(N_FTH+3));
  g_xth.Bthmod_CS = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*(N_FTH+3));
  g_xth.Bmax_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+3));
  g_xth.Bmin_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+3));

  // note that these are theta-periodic
  g_xth.Bth_B_xth_CS = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*N_FTH);
  g_xth.I_R2B_xth_CS = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*N_FTH);

  for(i = 0; i <= g_xth.nxx; i++) for(j = 0; j < N_FTH; j++)
  {
    ij = idx_xfth_CS(i,j);

    x = g_xth.xx[i];
    th = j*g_xth.dfth; th = mod_th(th);

    intp_RZ_xfth(x, th, &R, &Z);
    intp_dpsi_RZ(R, Z, pv); 
    Iv[0] = I_fn(R, Z);

    B = double_sqrt(pv[2]*pv[2] + Iv[0]*Iv[0] + pv[1]*pv[1])/R;

    ep = double_sqrt((R-1.0)*(R-1.0) + Z*Z);
    ep = max(ep, g_xth.ep_min*1.0e-8);
    dtdR = -double_sin(th)/ep;
    dtdZ = double_cos(th)/ep;

    if(j == 0) g_xth.Bmin_CS[idx_1d_CS(i)] = B;
    else if(j == N_FTH/2) g_xth.Bmax_CS[idx_1d_CS(i)] = B;

    g_xth.Bmod_CS[ij] = B;
    g_xth.Bthmod_CS[ij] = BP_SIGN*(pv[2]*dtdR - pv[1]*dtdZ)/R;
    

    ij = idx_xth_CS(i,j,g_xth.nxx,N_FTH);
    g_xth.Bth_B_xth_CS[ij] = BP_SIGN*(pv[2]*dtdR - pv[1]*dtdZ)/R;
    g_xth.I_R2B_xth_CS[ij] = Iv[0]/(R*R*B);
  }
  for(i = 0; i <= g_xth.nxx; i++)
  {
    g_xth.Bmod_CS[idx_xfth_CS(i,N_FTH)] = g_xth.Bmod_CS[idx_xfth_CS(i,0)];
    g_xth.Bthmod_CS[idx_xfth_CS(i,N_FTH)] = g_xth.Bthmod_CS[idx_xfth_CS(i,0)];
  }
  init_intp_2d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, 
      0.0, g_c.twopi, N_FTH, g_xth.Bmod_CS);
  init_intp_2d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, 
      0.0, g_c.twopi, N_FTH, g_xth.Bthmod_CS);
  
  init_intp_1d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, g_xth.Bmin_CS);
  init_intp_1d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, g_xth.Bmax_CS);

  init_intp_xth_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx,
      0.0, g_c.twopi, N_FTH, g_xth.Bth_B_xth_CS);
  init_intp_xth_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx,
      0.0, g_c.twopi, N_FTH, g_xth.I_R2B_xth_CS);
  /* set geomeric quantities for later use -------------------------------- */


//  my_mpi_barrier("setup_xx_theta_grid()");
}


void set_plasma_boundary(void)
{
  char fname[CHARLEN];
  FILE *fp;
  int j, jp;
  grid_pt_t *xc;
  double *xc_ep, *xc_th, *nrc_xc_ep, th, epx, R, Z, psi, depx, errx;

  xc = (grid_pt_t*)malloc(sizeof(grid_pt_t)*g_eqd.nbbbs);
  xc_th = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  xc_ep = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  nrc_xc_ep = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  assert(xc != NULL && xc_th != NULL && xc_ep != NULL && nrc_xc_ep != NULL);

  /* set x-curve : normalized & shifted */
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    xc[j].R = g_eqd.rbbbs[j]/g_eqd.Rc;
    xc[j].Z = (g_eqd.zbbbs[j] - g_eqd.Zc)/g_eqd.Rc;
  }

  /* sort curve point according to th */
  sort_grid_pt(xc, g_eqd.nbbbs);

  /* prepare ep interpolation along plasma boundary */ 
  jp = 0;
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    th = th_fn(xc[j].R, xc[j].Z);
    if(th < 0.0) th += g_c.twopi;
    else if(th >= g_c.twopi) th -= g_c.twopi;

    if(j == 0)
    {
      xc_th[jp] = th;
      xc_ep[jp] = double_sqrt((xc[j].R-1.0)*(xc[j].R-1.0) + xc[j].Z*xc[j].Z);
      jp++;
    }
    else
    {
      if(th > xc_th[jp-1])
      {
        xc_th[jp] = th;
        xc_ep[jp] = double_sqrt((xc[j].R-1.0)*(xc[j].R-1.0) + xc[j].Z*xc[j].Z);
        jp++;
      }
    }
  }

  //for(j = 0; j < g_eqd.nbbbs-1; j++) 
  for(j = 0; j < jp-1; j++) 
  {
    //assert(xc_th[j] < xc_th[j+1]);
    if(xc_th[j] > xc_th[j+1])
    {
      printf("wrong th interval : %d (%e, %e)\n", j, xc_th[j],xc_th[j+1]);
      exit(1);
    }
  }

  spline(xc_th, xc_ep, jp, 3.0e30, 3.0e30, nrc_xc_ep);


  if(N_NULLX == 1)
    depx = 0.5*absv(g_RZ.Zx[0] - g_RZ.Zmin);
  else
    depx = 0.5*min(absv(g_RZ.Zx[0]-g_RZ.Zmin), absv(g_RZ.Zx[1]-g_RZ.Zmax));

  errx = 1.0e-8;

  /* set xbry_(R,Z) : (R,Z) of plasma boundary */
  for(j = 0; j < N_FTH; j++)
  {
    th = j*g_xth.dfth;

    /* approximate ep of plasma boundary in this direction (th) */
    //if(th < xc_th[0] || th > xc_th[g_eqd.nbbbs-1])
    //  epx = 0.5*(xc_ep[0] + xc_ep[g_eqd.nbbbs-1]);
    if(th < xc_th[0] || th > xc_th[jp-1])
      epx = 0.5*(xc_ep[0] + xc_ep[jp-1]);
    else
      splint(xc_th, xc_ep, nrc_xc_ep, g_eqd.nbbbs, th, &epx);

    //epx *= 1.0+1.0e-3;
    R = epx*double_cos(th) + 1.0;
    Z = epx*double_sin(th);
    psi = psi_fn(R,Z);

    if(absv(psi-g_RZ.psi_last) < g_RZ.psi_last*errx) 
    {
      g_xth.xbry_R[j] = R;
      g_xth.xbry_Z[j] = Z;
    }
    else
    {
      find_psi_RZ(g_RZ.psi_last, R, Z, g_RZ.psi_last*errx, 
          &(g_xth.xbry_R[j]), &(g_xth.xbry_Z[j]));
    }
  }

  if(g_pall.ipe == 0)
  {
    sprintf(fname, "geo_xth-test3.txt", g_sys.work_dir);
    fp = fopen(fname, "w"); assert(fp != NULL);
    for(j = 0; j < g_eqd.nbbbs; j++)
      fprintf(fp, "%e %e\n", xc[j].R, xc[j].Z);
    fclose(fp);

    sprintf(fname, "geo_xth-test4.txt", g_sys.work_dir);
    fp = fopen(fname, "w"); assert(fp != NULL);
    for(j = 0; j < N_FTH; j++)
      fprintf(fp, "%e %e\n", g_xth.xbry_R[j], g_xth.xbry_Z[j]);

    fclose(fp);
  }

  free(xc);
  free(xc_ep);
  free(xc_th);
  free(nrc_xc_ep);


//:  my_mpi_barrier("set_plasma_boundary()");
}


void set_epm_R_Z_xfth(void)
{
  int i, j, ij, n_xx, n_fth;
  double Rx, Zx, R, Z, psi, err_min;
  char fname[CHARLEN];
  FILE *fp;

  err_min = 1.0e-10*g_RZ.psi_last;

  g_xth.epm_CS = (double*)malloc(sizeof(double)*(g_xth.nxx+3));
  g_xth.R_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*(N_FTH+3));
  g_xth.Z_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+3)*(N_FTH+3));
  g_xth.eqR_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(N_FTH+1));
  g_xth.eqZ_CS   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(N_FTH+1));
  /* initialize inverse coordinate mapping inside the plasma boundary 
     epm(psi)
     R(psi, th)
     Z(psi, th)
     */

  if(g_inpv.op_RZ_xfth_file == 0)
  {
    for(j = 0; j < N_FTH; j++)
    {
      /* set plasma boundary  // bug-fix July-17-2018
         g_xth.R_CS[idx_xfth_CS(g_xth.nxx,j)] = Rx = g_xth.xbry_R[j];
         g_xth.Z_CS[idx_xfth_CS(g_xth.nxx,j)] = Zx = g_xth.xbry_Z[j];*/

      Rx = g_xth.xbry_R[j];
      Zx = g_xth.xbry_Z[j];

      /* set nested flux surfaces inside plasma boundary */
      //for(i = g_xth.nxx - 1; i >= 0; i--)
      for(i = g_xth.nxx; i >= 0; i--)
      {
        ij = idx_xfth_CS(i,j);
        psi = g_xth.xx[i]*g_xth.xx[i]*g_RZ.psi_last;

        if(i == 0 && g_xth.xx0 < 1.0e-7) { R = 1.0; Z = 0.0; }
        else find_psi_RZ(psi, Rx, Zx, err_min, &R, &Z);

        g_xth.R_CS[ij] = R;
        g_xth.Z_CS[ij] = Z;
   // save contour
        g_xth.eqR_CS[i*(N_FTH+1)+j] = R;
        g_xth.eqZ_CS[i*(N_FTH+1)+j] = Z;
      }
    }
    /* set boundary */
    for(i = 0; i <= g_xth.nxx; i++)
    {
      g_xth.R_CS[idx_xfth_CS(i,N_FTH)] = g_xth.R_CS[idx_xfth_CS(i,0)];
      g_xth.Z_CS[idx_xfth_CS(i,N_FTH)] = g_xth.Z_CS[idx_xfth_CS(i,0)];
      g_xth.epm_CS[idx_1d_CS(i)] = g_xth.R_CS[idx_xfth_CS(i,0)]-1.0;
 
      g_xth.eqR_CS[i*(N_FTH+1)+N_FTH] = g_xth.eqR_CS[i*(N_FTH+1)];
      g_xth.eqZ_CS[i*(N_FTH+1)+N_FTH] = g_xth.eqZ_CS[i*(N_FTH+1)];
    }


    
    init_intp_1d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, g_xth.epm_CS);

    if(g_pall.ipe == 0)
    {
      sprintf(fname, "geo_RZ_xfth.dat", g_sys.work_dir);
      fp = fopen(fname, "w"); assert(fp != NULL);

      n_xx = g_xth.nxx;
      n_fth  = N_FTH;

      fwrite(&n_xx, sizeof(int), 1, fp);
      fwrite(&n_fth, sizeof(int), 1, fp);

      fwrite(g_xth.epm_CS, sizeof(double), g_xth.nxx+3, fp);

      fwrite(g_xth.R_CS, sizeof(double), (g_xth.nxx+3)*(N_FTH+3), fp);
      fwrite(g_xth.Z_CS, sizeof(double), (g_xth.nxx+3)*(N_FTH+3), fp);

      fclose(fp);

      // for test ----------------------------------------------
      sprintf(fname, "geo_RZ_xfth_test.txt", g_sys.work_dir);
      fp = fopen(fname, "w"); assert(fp != NULL);
      for(i = 0; i <= g_xth.nxx; i++)
      {
        fprintf(fp, "%i %e %e %e %e\n", i, g_xth.xx[i],
            g_xth.R_CS[idx_xfth_CS(i,0)]-1.0,
            g_xth.Z_CS[idx_xfth_CS(i,0)],
            intp_1d_CS(g_xth.xx[i], g_xth.xx0, g_xth.xx1, 
              g_xth.nxx, g_xth.epm_CS));
      }
      fclose(fp);
      // ------------------------------------------------------
    }
  }
  else
  {
    sprintf(fname, "RZ_xfth.dat", g_sys.work_dir);
    fp = fopen(fname, "r"); assert(fp != NULL);

    fread(&n_xx, sizeof(int), 1, fp);
    fread(&n_fth, sizeof(int), 1, fp);

    assert(n_xx == g_xth.nxx);
    assert(n_fth  == N_FTH);

    fread(g_xth.epm_CS, sizeof(double), g_xth.nxx+3, fp);

    fread(g_xth.R_CS, sizeof(double), (g_xth.nxx+3)*(N_FTH+3), fp);
    fread(g_xth.Z_CS, sizeof(double), (g_xth.nxx+3)*(N_FTH+3), fp);

    fclose(fp);
  }

  if (g_pall.ipe==0) {
  sprintf(fname, "geo_arr_bf_umfpack.txt", g_sys.work_dir);
  fp = fopen(fname, "w"); assert(fp != NULL);
  for(i = 0; i <= g_xth.nxx+2; i++)
  {
    for(j=0; j<=N_FTH+2;j++)
    {
      fprintf(fp, "%e %e\n",
        g_xth.R_CS[i*(N_FTH+3)+j]-1.0,
        g_xth.Z_CS[i*(N_FTH+3)+j]);
    }
  }
  fclose(fp);



  init_intp_2d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, 
      0.0, g_c.twopi, N_FTH, g_xth.R_CS);
  init_intp_2d_CS_umfpack(g_xth.xx0, g_xth.xx1, g_xth.nxx, 
      0.0, g_c.twopi, N_FTH, g_xth.Z_CS);

  sprintf(fname, "geo_arr_af_umfpack.txt", g_sys.work_dir);
  fp = fopen(fname, "w"); assert(fp != NULL);
  for(i = 0; i <= g_xth.nxx+2; i++)
  {
    for(j=0; j<=N_FTH+2;j++)
    {
      fprintf(fp, "%e %e\n",
        g_xth.R_CS[i*(N_FTH+3)+j]-1.0,
        g_xth.Z_CS[i*(N_FTH+3)+j]);
    }
  }
  fclose(fp);
  }

  //my_mpi_barrier("set_epm_R_Z_xfth()");
}



/* Find (R,Z) for given psi.
   Search is going along the line from plasma center to (R_s,Z_s) */
void find_psi_RZ(double psi, double R_s, double Z_s, double err, double *R, double *Z)
{
  int cnt = 0, cnt_max = g_xth.nxx*10000;
  double tmp_psi, R1, Z1, R2, Z2, Rt, Zt;
  double pv[6], th;

  R1 = 1.0;
  Z1 = 0.0;
  Rt = R2 = R_s; 
  Zt = Z2 = Z_s;

  tmp_psi = psi_fn(R2, Z2);

  while(double_fabs(psi-tmp_psi) > err && cnt < cnt_max)
  {
    cnt++;
    if(cnt > cnt_max)
    {
      printf("2M inter in find_psi_RZ() : %e (%e)\n", psi, tmp_psi);
      exit(1);
    }

    Rt = 0.5*(R1 + R2);
    Zt = 0.5*(Z1 + Z2);

    tmp_psi = psi_fn(Rt, Zt);
    if(tmp_psi > psi) { R2 = Rt; Z2 = Zt; }
    else { R1 = Rt; Z1 = Zt; }

    /*intp_dpsi_RZ(Rt, Zt, pv);
      tmp_psi = pv[0];
      th = th_fn(Rt, Zt);
      th = mod_th(th);
      pv[3] = double_cos(th)*pv[1] + double_sin(th)*pv[2];
      if(pv[3] >= 0.0)
      {
      if(tmp_psi > psi) { R2 = Rt; Z2 = Zt; }
      else { R1 = Rt; Z1 = Zt; }
      }
      else
      {
      if(tmp_psi > psi) { R1 = Rt; Z1 = Zt; }
      else { R2 = Rt; Z2 = Zt; }
      }*/
  }

  *R = Rt; *Z = Zt;
}



/* interpolate S-factor in fth-grid */
void S_fth_CS(int i, double th, double *S, double *dSdt)
{
  double tfj, C[4], dCdt, h, t, x1, x2;
  int tmpj, j[4], jp, ij;

#ifdef DEBUG1 
  assert(i >= 0 && i <= g_xth.nxx);
#endif

  tfj = th*g_xth.dfth_inv; tmpj = (int)tfj; 
  if(th < 0.0) tmpj--;
  h = tfj - tmpj; t = 1.0 - h;

  j[0] = tmpj - 1; 
  j[1] = tmpj;
  j[2] = tmpj + 1;
  j[3] = tmpj + 2;

  C[0] = t*t*t;
  C[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  C[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  C[3] = h*h*h;

  *S = *dSdt = 0.0;
  for(jp = 0; jp <= 3; jp++)
  {
    x1 = tfj - j[jp] + 0.5;
    x2 = tfj - j[jp] - 0.5;
    x1 = absv(x1);
    x2 = absv(x2);

    dCdt = (tsc(x1) - tsc(x2))*g_xth.dfth_inv;

    if(j[jp] >= 0 && j[jp] <= N_FTH) 
    {
      ij = idx_S_CS(i,j[jp]);

      *S += C[jp]*g_xth.S_CS[ij];
      *dSdt += dCdt*g_xth.S_CS[ij];
    }
    else if(j[jp] < 0) 
    {
      ij = idx_S_CS(i,j[jp]+N_FTH);

      *S += C[jp]*(g_xth.S_CS[ij] - g_xth.S2pi[i]);
      *dSdt += dCdt*(g_xth.S_CS[ij] - g_xth.S2pi[i]);
    }
    else 
    {
      ij = idx_S_CS(i,j[jp]-N_FTH);

      *S += C[jp]*(g_xth.S_CS[ij] + g_xth.S2pi[i]);
      *dSdt += dCdt*(g_xth.S_CS[ij] + g_xth.S2pi[i]);
    }
  }

  *S /= 6.0;
}


/* nrc-spline.c */

void spline(double *x, double *y, int n, const double yp1, const double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un, *u;

  u = (double*)malloc(sizeof(double)*n);

  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];

  free(u);
}

void splint(double *xa, double *ya, double *y2a, int n, const double x, double *y)
{
  int k, klo, khi;
  double h,b,a;

  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0)
  {
    printf("Bad xa input to routine splint");
    exit(1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
      +(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

inline int idx_1d_QS(int i)
{
        return i+1;
}

inline int idx_1d_CS(int i)
{
  return i+1;
}

// inline double mod_th(double x)
// {
//   if(x >= 0.0 && x < g_c.twopi)
//     return x;
//   else
//     if(x >= g_c.twopi)
//       return x - ((int)(x*g_c.twopi_inv))*g_c.twopi;
//     else
//       return x - ((int)(x*g_c.twopi_inv)-1)*g_c.twopi;
// }

/* 0th-spline */
inline double ngp(double x)
{
	if(x < 0.5)
		return 1.0;
	else
		return 0.0;
}

// /* 1st-spline */
// inline double cic(double x)
// {
// 	if(x < 1.0)
// 		return 1.0 - x;
// 	else
// 		return 0.0;
// }


// /* 2nd-spline */
// inline double tsc(double x)
// {
// 	if(x < 0.5)
// 		return 0.75 - x*x;
// 	else
// 		if(x < 1.5)
// 			return 0.5*(1.5-x)*(1.5-x);
// 		else
// 			return 0.0;
// }


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

// inline double th_fn(double R, double Z)
// {
//   if(R == 1.0)
//     if(Z == 0.0)
//       return 0.0;
//     else
//       if(Z < 0.0)
//         return 0.75*g_c.twopi;
//       else
//         return 0.25*g_c.twopi;
//   else
//     return double_atan2(Z,R-1.0);
// }

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

