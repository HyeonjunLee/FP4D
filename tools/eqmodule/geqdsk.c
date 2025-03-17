#include "cnst-srt.h"
#define idx_1d_CS(i) ((i)+1)
#define idx_2d_CS(i,j,nx,ny) ((j)+1 + (ny+3)*((i)+1))

void skip_line(FILE *fp, int nn);

/* read geqdsk file to setup equilibrium and make test outputs */
void read_geqdsk_and_test(void)
{
	int i, j, niter, inull1;
	double R, Z, psi, pv[6], gradp2, Xmin[2], psic, xpsi1, xgdp2;
	FILE *fp_eqd, *fp_test;
	char fname[CHARLEN];
	if(g_pall.ipe == 0)
	{
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("ALERT ALERT you should check nw & nh in read_geqdsk_file\n");
		printf("            read_geqdsk_and_test line 736               \n");
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	sprintf(fname, "geqdsk.dat", g_sys.work_dir);
	fp_eqd = fopen(fname, "r"); assert(fp_eqd != NULL);
        g_inpv.eqd_nw=65;
        g_inpv.eqd_nh=65;
	read_geqdsk_file(fp_eqd, g_inpv.eqd_nw, g_inpv.eqd_nh);

	/* test interpolation */
	if(g_pall.ipe == 0) 
	{
		/* ------------------------------------------------------------ */
		sprintf(fname, "test-eqd1.txt", g_sys.work_dir);
		fp_test = fopen(fname, "w"); assert(fp_test != NULL);
		for(i = 0; i <= g_eqd.NR; i++)
		{
			psi = g_eqd.psi0 + i*g_eqd.dpsi;

			fprintf(fp_test, "%e %e\n", psi,
					intp_1d_CS(psi, g_eqd.psi0, g_eqd.psi1, g_eqd.Npsi, g_eqd.I_CS));
		}
		fclose(fp_test);
		/* ------------------------------------------------------------ */


		/* ------------------------------------------------------------ */
		sprintf(fname, "test-eqd2.txt", g_sys.work_dir);
		fp_test = fopen(fname, "w"); assert(fp_test != NULL);
		for(i = 0; i <= g_eqd.NR; i++) 
		{
			R = g_eqd.Rmin + i*g_eqd.dR;

			for(j = 0; j <= g_eqd.NZ; j++)
			{
				Z = g_eqd.Zmin + j*g_eqd.dZ;

				/*pv[0] = intp_2d_CS(R, Z,
				  g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
				  g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
				  g_eqd.psi_2d_CS);*/

				intp_dfn_2d_CS(R, Z,
						g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
						g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
						g_eqd.psi_2d_CS, pv);

        gradp2 = pv[1]*pv[1] + pv[2]*pv[2] + 1.0e-9;

        fprintf(fp_test, "%e %e %e %e %e %e %e\n", R, Z, pv[0], pv[3], pv[4],
            intp_2d_CS(R,Z,g_eqd.Rmin,g_eqd.Rmax,g_eqd.NR,
              g_eqd.Zmin,g_eqd.Zmax,g_eqd.NZ,g_eqd.I_2d_CS),
            double_log(gradp2));
      }
      fprintf(fp_test, "\n");
    }
    fclose(fp_test);
    /* ------------------------------------------------------------ */


    /* ------------------------------------------------------------ */
    sprintf(fname, "test-lim.txt", g_sys.work_dir);
    fp_test = fopen(fname, "w"); assert(fp_test != NULL);
    for(i = 0; i < g_eqd.limitr; i++)
      fprintf(fp_test, "%e %e\n", g_eqd.rlim[i], g_eqd.zlim[i]);

    fprintf(fp_test, "\n\n");

    for(i = 0; i < g_eqd.nbbbs; i++)
      fprintf(fp_test, "%e %e\n", g_eqd.rbbbs[i], g_eqd.zbbbs[i]);
    fclose(fp_test);

    Xmin[0] = g_eqd.Rc;
    Xmin[1] = g_eqd.Zc;
    frprmn(Xmin, 2, 1.0e-30, &niter, &psic, 
        nrc_wrap_eqd_psi_fn, nrc_wrap_eqd_dpsi_fn);
    /*frprmn(Xmin, 2, 1.0e-30, &niter, &psic, 
      nrc_wrap_eqd_gdp2_fn, nrc_wrap_eqd_dgdp2_fn);*/
    printf("min psi at (%e, %e) = %e, (%d)\n", Xmin[0],Xmin[1],psic,niter);
    printf("psi at new axis = %e\n",
        intp_2d_CS(Xmin[0],Xmin[1], g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
          g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ, g_eqd.psi_2d_CS));


    find_min_gdp2_along_bbbs(0, g_eqd.nbbbs, &inull1, &xpsi1, &xgdp2);
    printf("first null : %d (%.10e, %.10e), psi = %.10e, gdp2 = %.10e\n",
        inull1, g_eqd.rbbbs[inull1], g_eqd.zbbbs[inull1],
        xpsi1, xgdp2);
    Xmin[0] = g_eqd.rbbbs[inull1];
    Xmin[1] = g_eqd.zbbbs[inull1];
    frprmn(Xmin, 2, 1.0e-30, &niter, &xgdp2, 
        nrc_wrap_eqd_gdp2_fn, nrc_wrap_eqd_dgdp2_fn);
    xpsi1 = nrc_wrap_eqd_psi_fn(Xmin);
    printf("  after refinement (%.10e, %.10e), psi = %.10e, gdp2 = %.10e\n",
        Xmin[0], Xmin[1], xpsi1, xgdp2);

    /* ------------------------------------------------------------ */

    write_my_geqdsk();
  }
}


/* ------------------------------------------------------------

   read original geqdsk file and set g_eqd.* 
   the following variables are initialized.

   geqdsk grid point number :
   g_eqd.nw, g_eqd.nh

   geqdsk box dimension and grid info :
   g_eqd.NR,   g_eqd.Rmin,   g_eqd.Rmax,   g_eqd.dR
   g_eqd.NZ,   g_eqd.Zmin,   g_eqd.Zmax,   g_eqd.dZ

   geqdsk psi grid info :
   g_eqd.Npsi, g_eqd.psi0, g_eqd.psi1, g_eqd.dpsi

   geqdsk plasma center :
   g_eqd.Rc, g_eqd.Zc, g_eqd.psic, g_eqd.Bc, g_eqd.Ic

   geqdsk plasma x-point (on the last flux surface) :
   g_eqd.Rx[N_NULLX], g_eqd.Zx[N_NULLX], g_eqd.psix[N_NULLX]

   psi & I 1d-interpolation as function of psi :
   g_eqd.I_CS[g_eqd.Npsi+3]

   psi & I 2d-interpolation as function of (R,Z) :
   g_eqd.psi_2d_CS[(g_eqd.NR+3)*(g_eqd.NZ+3)]
   g_eqd.I_2d_CS[(g_eqd.NR+3)*(g_eqd.NZ+3)]

   plasma boundary & limiter :
   g_eqd.nbbbs, g_eqd.limitr
   g_eqd.rbbbs[g_eqd.nbbbs], g_eqd.zbbbs[g_eqd.nbbbs]
   g_eqd.rlim[g_eqd.limitr], g_eqd.zlim[g_eqd.limitr] 

   -------------------------------------------------------------- */
void read_geqdsk_file(FILE *fp, int eqd_nw, int eqd_nh)
{
  int nw, nh, nhw, nbbbs, limitr, i, j, ix;
  double rdim, zdim, rcentr, rleft, zmid,
       rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  double *fpol, *pres, *ffprim, *pprime, *psirz, 
       *qpsi, *rbbbs, *zbbbs, *rlim, *zlim, pv[6];
  double R, Z, psi, ucf_len, ucf_flux, ucf_I;
  FILE *fp_test;
  char fname[CHARLEN]; //, s[MAX_LEN+1]

  /* keep MKS unit here */
  ucf_len = 1.0;
  ucf_flux = 1.0;
  ucf_I = 1.0;

  /* set manually */
  // fscanf(fp, "%s", &s [0]);
  // nw=atoi(s[51:54]);
  // nh=atoi(s[55:58]);
  g_eqd.nw = nw = eqd_nw;
  g_eqd.nh = nh = eqd_nh;
  nhw = nh*nw;


  /* allocate 1d array[nw] */
  fpol = (double*)malloc(sizeof(double)*nw);
  pres = (double*)malloc(sizeof(double)*nw);
  ffprim = (double*)malloc(sizeof(double)*nw);
  pprime = (double*)malloc(sizeof(double)*nw);
  qpsi = (double*)malloc(sizeof(double)*nw);

  /* allocate 2d array[nw*nh] */
  psirz = (double*)malloc(sizeof(double)*nw*nh);

  /* read from file ---------------------------------------------------- */
  skip_line(fp, 1); // skip one line

  fscanf(fp, "%lf %lf %lf %lf %lf",
      &rdim, &zdim, &rcentr, &rleft, &zmid);
  fscanf(fp, "%lf %lf %lf %lf %lf",
      &rmaxis, &zmaxis, &simag, &sibry, &bcentr);
  fscanf(fp, "%lf %lf %lf %lf %lf",
      &current, &simag, &xdum, &rmaxis, &xdum);
  fscanf(fp, "%lf %lf %lf %lf %lf",
      &zmaxis, &xdum, &sibry, &xdum, &xdum);

  if(g_pall.ipe == 0)
  {
    printf("RDIM = %e, ZDIM = %e, RLEFT = %e, ZMID = %e\n",
        rdim, zdim, rleft, zmid);
    printf("RMAXIS = %e, ZMAXIS = %e, SIMAG = %e, SIBRY = %e\n",
        rmaxis, zmaxis, simag, sibry);
    printf("RCENTR = %e, BCENTR = %e, CIRRENT = %e\n",
        rcentr, bcentr, current);
#ifdef PR_DEBUG
    fflush(stdout);
#endif
  }
  for(i = 0; i < nw; i++) fscanf(fp, "%lf", fpol+i);
  for(i = 0; i < nw; i++) fscanf(fp, "%lf", pres+i);
  for(i = 0; i < nw; i++) fscanf(fp, "%lf", ffprim+i);
  for(i = 0; i < nw; i++) fscanf(fp, "%lf", pprime+i);
  for(j = 0; j < nh; j++) for(i = 0; i < nw; i++)
    fscanf(fp, "%lf", psirz + (j + nh*i));

  for(i = 0; i < nw; i++) fscanf(fp, "%lf", qpsi+i);
  
  g_eqd.simag=simag;
  /* shift psi by simag */
  for(i = 0; i < nhw; i++) psirz[i] -= simag;
  sibry -= simag;
  simag = 0.0;

  /* fix fpol and psirz sign */
  for(i = 0; i < nhw; i++) psirz[i] = absv(psirz[i]);
  for(i = 0; i < nw; i++) fpol[i] = absv(fpol[i]);
  simag = absv(simag);
  sibry = absv(sibry);


  fscanf(fp, "%d %d", &nbbbs, &limitr);
  /* allocate 1d array[nbbbs] */
  rbbbs = (double*)malloc(sizeof(double)*nbbbs);
  zbbbs = (double*)malloc(sizeof(double)*nbbbs);
  /* allocate 1d array[limitr] */
  rlim = (double*)malloc(sizeof(double)*limitr);
  zlim = (double*)malloc(sizeof(double)*limitr);

  for(i = 0; i < nbbbs; i++)  fscanf(fp, "%lf %lf", rbbbs+i, zbbbs+i);
  for(i = 0; i < limitr; i++) fscanf(fp, "%lf %lf", rlim+i, zlim+i);
  /* read from file ---------------------------------------------------- */

  /* number of data for interpolation is num = eqd_num - 1. 
     So we have num + 3 = eqd_num + 2 */
  g_eqd.NR   = g_eqd.nw-1;
  g_eqd.NZ   = g_eqd.nh-1;
  g_eqd.Npsi = g_eqd.nw-1;

  /* set rectangular eqd-domain information */
  g_eqd.Rmin = ucf_len*rleft;
  g_eqd.Rmax = ucf_len*(rleft + rdim);
  g_eqd.Zmin = ucf_len*(zmid - 0.5*zdim);
  g_eqd.Zmax = ucf_len*(zmid + 0.5*zdim);
  g_eqd.psi0 = ucf_flux*simag;
  g_eqd.psi1 = ucf_flux*sibry;
  g_eqd.dR = (g_eqd.Rmax - g_eqd.Rmin)/g_eqd.NR;
  g_eqd.dZ = (g_eqd.Zmax - g_eqd.Zmin)/g_eqd.NZ;
  g_eqd.dpsi = (g_eqd.psi1 - g_eqd.psi0)/g_eqd.Npsi;

  /* allocate variables */
  g_eqd.I_CS = (double*)malloc(sizeof(double)*(g_eqd.Npsi+3));
  g_eqd.psi_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));
  g_eqd.I_2d_CS = (double*)malloc(sizeof(double)*(g_eqd.NR+3)*(g_eqd.NZ+3));
  /* initialize interpolation */
  for(i = 0; i <= g_eqd.Npsi; i++) g_eqd.I_CS[idx_1d_CS(i)] = ucf_I*fpol[i];
  init_intp_1d_CS_umfpack(g_eqd.psi0, g_eqd.psi1, g_eqd.Npsi, g_eqd.I_CS);
#ifdef PR_DEBUG
  my_mpi_barrier("eqdsk: pass-1");
#endif
  sprintf(fname, "test-geqdsk2.txt", g_sys.work_dir);
  fp_test = fopen(fname, "w"); 
  assert(fp_test != NULL);
  for(i = 0; i <= g_eqd.NR; i++) 
  {
    for(j = 0; j <= g_eqd.NZ; j++)
    {
      psi = ucf_flux*psirz[j+nh*i];

      g_eqd.psi_2d_CS[idx_2d_CS(i,j,g_eqd.NR,g_eqd.NZ)] = psi;

      fprintf(fp_test, "%e %e %e\n", 
          1.0e2*(g_eqd.Rmin + i*g_eqd.dR), 
          1.0e2*(g_eqd.Zmin + j*g_eqd.dZ), psi*1.0e8);
    }
    fprintf(fp_test, "\n");
  }
  fclose(fp_test);

  init_intp_2d_CS_umfpack(g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ, 
      g_eqd.psi_2d_CS);

#ifdef PR_DEBUG
  my_mpi_barrier("eqdsk: pass-2");
#endif

  /* set plasma center */
  g_eqd.Rc = ucf_len*rmaxis;
  g_eqd.Zc = ucf_len*zmaxis;
  intp_dfn_2d_CS(g_eqd.Rc, g_eqd.Zc, 
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);
  g_eqd.psic = pv[0];
  g_eqd.Ic = intp_1d_CS(pv[0],g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  g_eqd.Bc = double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+g_eqd.Ic*g_eqd.Ic)/g_eqd.Rc;

#ifdef PR_DEBUG
  my_mpi_barrier("eqdsk: pass-3");
#endif

  /* initialize global variables describing edge geometry */
  g_eqd.nbbbs = nbbbs;
  g_eqd.rbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  g_eqd.zbbbs = (double*)malloc(sizeof(double)*g_eqd.nbbbs);
  g_eqd.rbbbs_max = g_eqd.rbbbs_min = g_eqd.Rc;
  g_eqd.zbbbs_max = g_eqd.zbbbs_min = g_eqd.Zc;
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    g_eqd.rbbbs[j] = ucf_len*rbbbs[j];
    g_eqd.zbbbs[j] = ucf_len*zbbbs[j];

    g_eqd.rbbbs_max = max(g_eqd.rbbbs_max, g_eqd.rbbbs[j]);
    g_eqd.rbbbs_min = min(g_eqd.rbbbs_min, g_eqd.rbbbs[j]);

    g_eqd.zbbbs_max = max(g_eqd.zbbbs_max, g_eqd.rbbbs[j]);
    g_eqd.zbbbs_min = min(g_eqd.zbbbs_min, g_eqd.rbbbs[j]);
  }

  /* set minor radius */
  g_eqd.a0 = g_eqd.rbbbs_max - g_eqd.Rc;

  g_eqd.limitr = limitr;
  g_eqd.rlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  g_eqd.zlim = (double*)malloc(sizeof(double)*g_eqd.limitr);
  for(j = 0; j < g_eqd.limitr; j++)
  {
    g_eqd.rlim[j] = ucf_len*rlim[j];
    g_eqd.zlim[j] = ucf_len*zlim[j];
  }

  /* set x-point */
  for(i = 0; i < N_NULLX; i++)
  {
#ifdef PR_DEBUG
    my_mpi_barrier("eqdsk: find x-point");
#endif

    j = 0;
    g_eqd.Rx[i] = g_eqd.rbbbs[j];
    g_eqd.Zx[i] = g_eqd.zbbbs[j];
    g_eqd.ix[i] = j;

    if(i == 0)
    {
      /* lower null */
      for(j = 1; j < g_eqd.nbbbs; j++)
      {
        if(g_eqd.Zx[i] > g_eqd.zbbbs[j])
        {
          g_eqd.Rx[i] = g_eqd.rbbbs[j];
          g_eqd.Zx[i] = g_eqd.zbbbs[j];
          g_eqd.ix[i] = j;
        }
      }
    }
    else
    {
      /* upper null */
      for(j = 1; j < g_eqd.nbbbs; j++)
      {
        if(g_eqd.Zx[i] < g_eqd.zbbbs[j])
        {
          g_eqd.Rx[i] = g_eqd.rbbbs[j];
          g_eqd.Zx[i] = g_eqd.zbbbs[j];
          g_eqd.ix[i] = j;
        }
      }
    }

    if(g_pall.ipe == 0)
    {
      printf("===> X-point(%d) = (%e, %e)\n", i, g_eqd.Rx[i], g_eqd.Zx[i]);
    }
  }

  /* for double null case, sort x-point according to Z */
  if(N_NULLX == 2)
  {
    if(g_eqd.Zx[0] > g_eqd.Zx[1]) 
    {
      pv[0] = g_eqd.Rx[0];
      pv[1] = g_eqd.Zx[0];
      pv[2] = g_eqd.psix[0];
      ix = g_eqd.ix[0];

      g_eqd.Rx[0] = g_eqd.Rx[1];
      g_eqd.Zx[0] = g_eqd.Zx[1];
      g_eqd.psix[0] = g_eqd.psix[1];
      g_eqd.ix[0] = g_eqd.ix[1];

      g_eqd.Rx[1] = pv[0];
      g_eqd.Zx[1] = pv[1];
      g_eqd.psix[1] = pv[2];
      g_eqd.ix[1] = ix;
    }
  }

  /* set circular approximate dimension : just for comparison */
  g_cir.a0 = g_eqd.a0*1.0e2;
  g_cir.R0 = g_eqd.Rc*1.0e2;
  g_cir.B0 = g_eqd.Bc*1.0e4;

  refine_my_geqdsk();

  /* initialize I_2d_CS interpolation.
     In geqdsk, I is defined only inside of plasma boundary (last closed
     surface). For the simulation, the definition is extended to cover
     entire R, Z box. Outside of plasma value is set equal to the value
     at the last closed surface */
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

  /* free allocations */
  free(fpol); free(pres); free(ffprim); free(pprime); free(qpsi);
  free(psirz);
  free(rbbbs); free(zbbbs); free(rlim); free(zlim);
}


/* convert geqdsk variables from MKS to cgs */
void my_geqdsk_to_cgs(void)
{
  int i, j;
  double ucf_len, ucf_flux, ucf_I, ucf_B;

  ucf_len = 1.0e2;
  ucf_flux = 1.0e4*1.0e2*1.0e2;
  ucf_I = 1.0e4*1.0e2;
  ucf_B = 1.0e4;

  g_eqd.Rmin *= ucf_len;
  g_eqd.Rmax *= ucf_len;
  g_eqd.Zmin *= ucf_len;
  g_eqd.Zmax *= ucf_len;
  g_eqd.dR *= ucf_len;
  g_eqd.dZ *= ucf_len;

  g_eqd.psi0 *= ucf_flux;
  g_eqd.psi1 *= ucf_flux;
  g_eqd.dpsi *= ucf_flux;

  g_eqd.a0 *= ucf_len;
  g_eqd.Rc *= ucf_len;
  g_eqd.Zc *= ucf_len;
  g_eqd.psic *= ucf_flux;
  g_eqd.Bc *= ucf_B;
  g_eqd.Ic *= ucf_I;

  for(j = 0; j < N_NULLX; j++)
  {
    g_eqd.Rx[j] *= ucf_len;
    g_eqd.Zx[j] *= ucf_len;
    g_eqd.psix[j] *= ucf_flux;
  }

  for(i = 0; i < g_eqd.Npsi+3; i++)
    g_eqd.I_CS[i] *= ucf_I;

  for(i = 0; i < (g_eqd.NR+3)*(g_eqd.NZ+3); i++)
  {
    g_eqd.I_2d_CS[i] *= ucf_I;
    g_eqd.psi_2d_CS[i] *= ucf_flux;
  }

  g_eqd.rbbbs_max *= ucf_len;
  g_eqd.rbbbs_min *= ucf_len;
  g_eqd.zbbbs_max *= ucf_len;
  g_eqd.zbbbs_min *= ucf_len;
  for(j = 0; j < g_eqd.nbbbs; j++)
  {
    g_eqd.rbbbs[j] *= ucf_len;
    g_eqd.zbbbs[j] *= ucf_len;
  }

  for(j = 0; j < g_eqd.limitr; j++)
  {
    g_eqd.rlim[j] *= ucf_len;
    g_eqd.zlim[j] *= ucf_len;
  }
}


void refine_my_geqdsk(void)
{
  double x[2], fval, pv[6];
  int niter;

  /* set and refine plasma center ---------------------------------------- */
  if(g_pall.ipe == 0) printf("refine plasma center....\n");

  /* guess from original geqdsk file */
  x[0] = g_eqd.Rc;
  x[1] = g_eqd.Zc;

  if(g_pall.ipe == 0)
  {
    printf("before refinement : psi(%.10e,%.10e) = %.10e\n", 
        x[0], x[1], nrc_wrap_eqd_psi_fn(x));
    printf("Bc, Ic = %.10e %.10e\n\n", g_eqd.Bc, g_eqd.Ic);
  }

  /* minimize psi to refine plasma center */
  frprmn(x, 2, 1.0e-30, &niter, &fval, 
      nrc_wrap_eqd_psi_fn, nrc_wrap_eqd_dpsi_fn);

  /* refined plasma center */
  g_eqd.Rc = x[0];
  g_eqd.Zc = x[1];

  intp_dfn_2d_CS(g_eqd.Rc, g_eqd.Zc, 
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);

  /* readjust uniform psi grid */
  g_eqd.psi0 = pv[0];

  /* refined psi & B at plasma center */
  g_eqd.psic = pv[0];
  g_eqd.Ic = intp_1d_CS(pv[0],g_eqd.psi0,g_eqd.psi1,g_eqd.Npsi,g_eqd.I_CS);
  g_eqd.Bc = double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+g_eqd.Ic*g_eqd.Ic)/g_eqd.Rc;

  if(g_pall.ipe == 0)
  {
    printf("after refinement(%d) : psi(%.10e,%.10e) = %.10e(%.10e)\n", 
        niter, g_eqd.Rc, g_eqd.Zc, g_eqd.psic, fval);
    printf("Bc, Ic = %.10e %.10e\n\n", g_eqd.Bc, g_eqd.Ic);
#ifdef PR_DEBUG
    fflush(stdout);
#endif
  }
  /* set and refine plasma center ---------------------------------------- */



  /* set and refine x-points --------------------------------------------- 
     for(i = 0; i < N_NULLX; i++)
     {
     if(g_pall.ipe == 0)
     printf("refine plasma x-point(%d)....\n", i);

  // guess of x-point by naive search along g_eqd(r,z)bbbx 
  x[0] = g_eqd.Rx[i];
  x[1] = g_eqd.Zx[i];

  if(g_pall.ipe == 0)
  printf("before refinement : psi(%.10e,%.10e) = %.10e\n", 
  x[0], x[1], nrc_wrap_eqd_psi_fn(x));

  // minimiza gradpsi^2 to refine x-point 
  frprmn(x, 2, 1.0e-30, &niter, &fval, 
  nrc_wrap_eqd_gdp2_fn, nrc_wrap_eqd_dgdp2_fn);
  g_eqd.rbbbs[g_eqd.ix[i]] = g_eqd.Rx[i] = x[0];
  g_eqd.zbbbs[g_eqd.ix[i]] = g_eqd.Zx[i] = x[1];
  g_eqd.psix[i] = nrc_wrap_eqd_psi_fn(x);
  g_eqd.psi1 = max(g_eqd.psi1, g_eqd.psix[i]);

  if(g_pall.ipe == 0)
  {
  printf("after refinement(%d) : psi(%.10e,%.10e) = %.10e (%.10e)\n", 
  niter, g_eqd.Rx[i], g_eqd.Zx[i], g_eqd.psix[i], fval);
#ifdef PR_DEBUG
fflush(stdout);
#endif
}
}*/

/* for double null case */
if(N_NULLX == 2) assert(g_eqd.Zx[0] < g_eqd.Zx[1]);

/* rescale uniform psi grid */
g_eqd.dpsi = (g_eqd.psi1 - g_eqd.psi0)/g_eqd.Npsi;
/* set and refine x-points --------------------------------------------- */
}


void write_my_geqdsk(void)
{
  FILE *fp_eqd;
  int i, j, nw, nh, npsi, idum;
  double pv[6];
  double xdum, rdim, zdim, rcenter, rleft, zmid, 
       rmaxis, zmaxis, simag, sibry, bcenter, current;
  double ucf_len, ucf_flux, ucf_I, ucf_B;

  ucf_len = 1.0e2;
  ucf_flux = 1.0e4*1.0e2*1.0e2;
  ucf_I = 1.0e4*1.0e2;
  ucf_B = 1.0e4;

  idum = 0;
  nw = g_eqd.NR+1;
  nh = g_eqd.NZ+1;
  assert(g_eqd.NR == g_eqd.Npsi);

  rdim = g_eqd.Rmax - g_eqd.Rmin;
  rleft = g_eqd.Rmin;

  zmid = 0.5*(g_eqd.Zmax + g_eqd.Zmin);
  zdim = g_eqd.Zmax - g_eqd.Zmin;

  rmaxis = rcenter = g_eqd.Rc;
  zmaxis = 0.0;

  bcenter = g_eqd.Bc;
  simag = g_eqd.psi0;
  sibry = g_eqd.psi1;

  current = 0.0;
  xdum = 0.0;

  rdim /= ucf_len;
  rleft /= ucf_len;
  zdim /= ucf_len;
  zmid /= ucf_len;
  rmaxis /= ucf_len;
  zmaxis /= ucf_len;
  rcenter /= ucf_len;
  bcenter /= ucf_B;
  simag /= ucf_flux;
  sibry /= ucf_flux;

  assert(g_pall.ipe == 0);

  fp_eqd = fopen("my_geqdsk.txt", "w");
  assert(fp_eqd != NULL);

  fprintf(fp_eqd, "EFITD    01/23/2002    # 13728  3300ms   ");
  fprintf(fp_eqd, "%d %d %d\n", idum, nw, nh);
  fprintf(fp_eqd, "%.10e %.10e %.10e %.10e %.10e\n",
      rdim, zdim, rcenter, rleft, zmid);
  fprintf(fp_eqd, "%.10e %.10e %.10e %.10e %.10e\n",
      rmaxis, zmaxis, simag, sibry, bcenter);
  fprintf(fp_eqd, "%.10e %.10e %.10e %.10e %.10e\n",
      current, simag, xdum, rmaxis, xdum);
  fprintf(fp_eqd, "%.10e %.10e %.10e %.10e %.10e\n",
      zmaxis, xdum, sibry, xdum, xdum);

  /* fpol */
  for(i = 0; i <= g_eqd.Npsi; i++)
  {
    fprintf(fp_eqd, "%.10e ",
        intp_1d_CS(g_eqd.psi0+i*g_eqd.dpsi, g_eqd.psi0, g_eqd.psi1,
          g_eqd.Npsi, g_eqd.I_CS)/ucf_I);
  } 
  fprintf(fp_eqd, "\n");

  /* pres */
  for(i = 0; i <= g_eqd.Npsi; i++)
  {
    fprintf(fp_eqd, "%.10e ", 0.0);
  } 
  fprintf(fp_eqd, "\n");

  /* ffprim */
  for(i = 0; i <= g_eqd.Npsi; i++)
  {
    fprintf(fp_eqd, "%.10e ", 0.0);
  } 
  fprintf(fp_eqd, "\n");

  /* pprime */
  for(i = 0; i <= g_eqd.Npsi; i++)
  {
    fprintf(fp_eqd, "%.10e ", 0.0);
  } 
  fprintf(fp_eqd, "\n");

  /* psirz */
  for(i = 0; i <= g_eqd.NR; i++) 
  {
    for(j = 0; j <= g_eqd.NZ; j++)
    {
      pv[0] = intp_2d_CS(g_eqd.Rmin+i*g_eqd.dR,g_eqd.Zmin+j*g_eqd.dZ,
          g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
          g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
          g_eqd.psi_2d_CS);
      assert(pv[0] >= g_eqd.psi0);

      fprintf(fp_eqd, "%.10e ", pv[0]/ucf_flux);
    }
  }
  fprintf(fp_eqd, "\n");

  /* qpsi */
  for(i = 0; i <= g_eqd.Npsi; i++)
  {
    fprintf(fp_eqd, "%.10e ", 0.0);
  } 
  fprintf(fp_eqd, "\n");


  fprintf(fp_eqd, "%d %d\n", g_eqd.nbbbs, g_eqd.limitr);
  /* rbbbs, zbbbs */
  for(i = 0; i < g_eqd.nbbbs; i++)
  {
    fprintf(fp_eqd, "%.10e %.10e", 
        g_eqd.rbbbs[i]/ucf_len, g_eqd.zbbbs[i]/ucf_len);
  } 
  fprintf(fp_eqd, "\n");
  /* rlim, zlim */
  for(i = 0; i < g_eqd.limitr; i++)
  {
    fprintf(fp_eqd, "%.10e %.10e", 
        g_eqd.rlim[i]/ucf_len, g_eqd.zlim[i]/ucf_len);
  } 
  fprintf(fp_eqd, "\n");

  fclose(fp_eqd);
}


/* compare theta of given grid points */
int comp_grid_pt(const void *x, const void *y)
{
  double R1, Z1, R2, Z2, th1, th2;

  R1 = ((grid_pt_t*)x)->R;
  Z1 = ((grid_pt_t*)x)->Z;

  R2 = ((grid_pt_t*)y)->R;
  Z2 = ((grid_pt_t*)y)->Z;

  th1 = th_fn(R1, Z1);
  if(th1 < 0.0) th1 += g_c.twopi;
  else if(th1 >= g_c.twopi) th1 -= g_c.twopi;

  th2 = th_fn(R2, Z2);
  if(th2 < 0.0) th2 += g_c.twopi;
  else if(th2 >= g_c.twopi) th2 -= g_c.twopi;

  if(th1 < th2) return -1;
  else if(th1 == th2) return 0;
  else return 1;
}

/* sort grid points accroding to theta value */
void sort_grid_pt(grid_pt_t *curve, int num)
{
  qsort(curve, num, sizeof(grid_pt_t), comp_grid_pt);
}


/* nrc wrappers for various minization ------------------------------------ */
double nrc_wrap_eqd_psi_fn(double *x)
{
  return intp_2d_CS(x[0], x[1],
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS);
}
void nrc_wrap_eqd_dpsi_fn(double *x, double *dy)
{
  double pv[6];

  intp_dfn_2d_CS(x[0], x[1],
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);

  dy[0] = pv[1];
  dy[1] = pv[2];
}
double nrc_wrap_eqd_gdp2_fn(double *x)
{
  double pv[6];
  //double err=1.0e-10;

  //printf("gdp2 : (%e,%e)\n", x[0], x[1]);
  intp_dfn_2d_CS(x[0], x[1],
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);

  //return double_log(pv[1]*pv[1] + pv[2]*pv[2] + err);
  return pv[1]*pv[1] + pv[2]*pv[2];
}
void nrc_wrap_eqd_dgdp2_fn(double *x, double *dy)
{
  double pv[6];
  //double err=1.0e-10, dem;

  //printf("dgdp2 : (%e,%e)\n", x[0], x[1]);
  intp_dfn_2d_CS(x[0], x[1],
      g_eqd.Rmin, g_eqd.Rmax, g_eqd.NR,
      g_eqd.Zmin, g_eqd.Zmax, g_eqd.NZ,
      g_eqd.psi_2d_CS, pv);

  /*dem = pv[1]*pv[1] + pv[2]*pv[2] + err;
    dy[0] = 2.0*(pv[1]*pv[3] + pv[2]*pv[5])/dem;
    dy[1] = 2.0*(pv[1]*pv[5] + pv[2]*pv[4])/dem;*/
  dy[0] = 2.0*(pv[1]*pv[3] + pv[2]*pv[5]);
  dy[1] = 2.0*(pv[1]*pv[5] + pv[2]*pv[4]);
}
/* nrc wrappers for various minization ------------------------------------ */



/* find minimum grad psi^2 point along g_eqd(r,z)bbbs */
void find_min_gdp2_along_bbbs(int istart, int num, 
    int *nmin, double *psi, double *gdp2)
{
  int i, ip;
  double val, x[2];

  while(istart >= g_eqd.nbbbs) istart -= g_eqd.nbbbs;
  while(istart < 0) istart += g_eqd.nbbbs;

  x[0] = g_eqd.rbbbs[istart];
  x[1] = g_eqd.zbbbs[istart];
  *nmin = istart;
  *gdp2 = nrc_wrap_eqd_gdp2_fn(x);

  for(i = istart; i < istart + num; i++)
  {
    ip = i;
    while(ip >= g_eqd.nbbbs) ip -= g_eqd.nbbbs;
    while(ip < 0) ip += g_eqd.nbbbs;

    x[0] = g_eqd.rbbbs[ip];
    x[1] = g_eqd.zbbbs[ip];
    val = nrc_wrap_eqd_gdp2_fn(x);

    if(*gdp2 >= val) 
    {
      *gdp2 = val;
      *nmin = ip;
    }
  }

  assert(*nmin >= 0 && *nmin < g_eqd.nbbbs);
  x[0] = g_eqd.rbbbs[*nmin];
  x[1] = g_eqd.zbbbs[*nmin];
  *psi = nrc_wrap_eqd_psi_fn(x);
}


/* nrc-miz.c */

void frprmn(double *p, int n, double ftol, int *iter, double *fret,
    double (*func)(double []), void (*dfunc)(double [], double []));
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));
double f1dim(double x);
double brent(double ax, double bx, double cx, double (*f)(double), double tol,
    double *xmin);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
    double (*func)(double));

/* frprmn.c ----------------------------------------------------------------- */
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free(xi);free(h);free(g);

void frprmn(double *p, int n, double ftol, int *iter, double *fret,
    double (*func)(double []), void (*dfunc)(double [], double []))
{
  int j,its;
  double gg,gam,fp,dgg;
  double *g,*h,*xi;

  g = (double*)malloc(sizeof(double)*n);
  h = (double*)malloc(sizeof(double)*n);
  xi = (double*)malloc(sizeof(double)*n);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  for (j=0;j<n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      FREEALL
        return;
    }
    fp=(*func)(p);
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=0;j<n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      FREEALL
        return;
    }
    gam=dgg/gg;
    for (j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }

  printf("Too many iterations in frprmn");
  exit(1);
}
#undef ITMAX
#undef EPS
#undef FREEALL
/* frprmn.c ----------------------------------------------------------------- */


/* linmin.c ----------------------------------------------------------------- */
#define TOL 2.0e-4

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;

  ncom=n;
  pcom=(double*)malloc(sizeof(double)*n);
  xicom=(double*)malloc(sizeof(double)*n);
  nrfunc=func;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free(xicom);
  free(pcom);
}
#undef TOL
/* linmin.c ----------------------------------------------------------------- */


/* f1dim.c ------------------------------------------------------------------ */
double f1dim(double x)
{
  int j;
  double f,*xt;

  xt=(double*)malloc(sizeof(double)*ncom);
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free(xt);
  return f;
}
/* f1dim.c ------------------------------------------------------------------ */


/* brent.c ------------------------------------------------------------------ */
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double brent(double ax, double bx, double cx, double (*f)(double), double tol,
    double *xmin)
{
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
        d=p/q;
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
        SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  printf("Too many iterations in brent");
  *xmin=x;
  return fx;
}
#undef SIGN
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
/* brent.c ------------------------------------------------------------------ */


/* mnbrak.c ----------------------------------------------------------------- */
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
    double (*func)(double))
{
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(max(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
          SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
  }
}
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
/* mnbrak.c ----------------------------------------------------------------- */



/* (C) Copr. 1986-92 Numerical Recipes Software *!%-V. */













