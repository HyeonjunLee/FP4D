#include "cnst-srt.h"
#ifdef FLAG_KAIROS
  #include <umfpack.h>
#else // hnuc
  #include <suitesparse/umfpack.h>
#endif

#define idx_1d_CS(i) ((i)+1)
#define idx_2d_CS(i,j,nx,ny) ((j)+1 + (ny+3)*((i)+1))
#define idx_xth_CS(i,j,nx,nt) ((j) + (nt)*((i)+1))
/* the 2nd dimension is periodic 
 *    (i,j) = [-1, nx+1] x [0, ny-1] ==> total (nx+3)*nt */

void constr_umfpackM_1d_CS(int nn, int **Ap, int **Ai, double **Ax, void **Num)
{
  int i, ip, ii;
  size_t dim;
  void *sym;
  double  C[3] = {1.0/6.0, 4.0/6.0, 1.0/6.0},
        //dC[3] = {1.0, 0.0, -1.0},
        ddC[3] = {1.0, -2.0, 1.0};

  dim = nn + 3;

  *Ap = (int*)malloc(sizeof(int)*(dim+1));
  *Ai = (int*)malloc(sizeof(int)*(dim*9));
  *Ax = (double*)malloc(sizeof(double)*(dim*9));

  (*Ap)[0] = 0;
  ip = 0;
  for(i = -1; i <= nn+1; i++) 
  {
    (*Ap)[idx_1d_CS(i)+1] = (*Ap)[idx_1d_CS(i)]+3;

    if(i == -1)
    {
      for(ii = -1; ii <= 1; ii++) 
      {
        (*Ai)[ip] = idx_1d_CS(ii);
        (*Ax)[ip] = ddC[ii+1];
        ip++;
      }
    }
    else if(i == nn+1)
    {
      for(ii = nn-1; ii <= nn+1; ii++) 
      {
        (*Ai)[ip] = idx_1d_CS(ii);
        (*Ax)[ip] = ddC[ii-nn+1];
        ip++;
      }
    }
    else
    {
      for(ii = i-1; ii <= i+1; ii++) 
      {
        (*Ai)[ip] = idx_1d_CS(ii);
        (*Ax)[ip] = C[ii-i+1];
        ip++;
      }
    }
  }

  assert(ip == 3*dim);

  umfpack_di_symbolic(dim,dim,*Ap,*Ai,*Ax,&sym,NULL,NULL);
  umfpack_di_numeric(*Ap,*Ai,*Ax,sym,Num,NULL,NULL);
  umfpack_di_free_symbolic(&sym);
}




void init_intp_1d_CS_umfpack(double x1, double x2, int nn, double *data)
{
  size_t dim;
  int i, ip, *Ap, *Ai;
  double *Ax, *b;
  void *Num;

  dim = nn + 3;

  constr_umfpackM_1d_CS(nn, &Ap, &Ai, &Ax, &Num);

  b = (double*)malloc(sizeof(double)*dim);

  // set components of source vector ------------------------------------
  for(i = -1; i <= nn+1; i++) 
  {
    ip = idx_1d_CS(i);

    if(i == -1 || i == nn+1) b[ip] = 0.0;
    else b[ip] = data[ip];
  }
  // set components of source vector ------------------------------------

  umfpack_di_solve(UMFPACK_At, Ap, Ai, Ax, data, b, Num, NULL, NULL);

  umfpack_di_free_numeric(&Num);
  free(Ap);
  free(Ai);
  free(Ax);
  free(b);
}




void constr_umfpackM_2d_CS(int nx, int ny, 
    int **Ap, int **Ai, double **Ax, void **Num)
  /*void constr_umfpackM_2d_CS(int nx, int ny, 
    int *Ap, int *Ai, double *Ax, void *Num)*/
{
  int i, j, ij, ijp, ii, jj;
  size_t dim;
  void *sym;
  double  C[3] = {1.0/6.0, 4.0/6.0, 1.0/6.0},
        dC[3] = {1.0, 0.0, -1.0},
        ddC[3] = {1.0, -2.0, 1.0};

  dim = (nx+3)*(ny+3);

  *Ap = (int*)malloc(sizeof(int)*(dim+1));
  *Ai = (int*)malloc(sizeof(int)*(dim*9));
  *Ax = (double*)malloc(sizeof(double)*(dim*9));

  (*Ap)[0] = 0;
  ijp = 0;
  for(i = -1; i <= nx+1; i++) for(j = -1; j <= ny+1; j++)
  {
    ij = idx_2d_CS(i, j, nx, ny);
    (*Ap)[ij+1] = (*Ap)[ij]+9;
    //Ap[ij+1] = Ap[ij]+9;

    if(i == -1 || i == nx+1 || j == -1 || j == ny+1)
    {
      if(i == -1 && j == -1)
      {
        for(ii = -1; ii <= 1; ii++) for(jj = -1; jj <= 1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = dC[ii+1]*dC[jj+1];
          ijp++;
        }
      }
      else if(i == -1 && j == ny+1)
      {
        for(ii = -1; ii <= 1; ii++) for(jj = ny-1; jj <= ny+1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = dC[ii+1]*dC[jj-ny+1];
          ijp++;
        }
      }
      else if(i == nx+1 && j == ny+1)
      {
        for(ii = nx-1; ii <= nx+1; ii++) for(jj = ny-1; jj <= ny+1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = dC[ii-nx+1]*dC[jj-ny+1];
          ijp++;
        }
      }
      else if(i == nx+1 && j == -1)
      {
        for(ii = nx-1; ii <= nx+1; ii++) for(jj = -1; jj <= +1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = dC[ii-nx+1]*dC[jj+1];
          ijp++;
        }
      }
      else if(i == -1)
      {
        for(ii = -1; ii <= 1; ii++) for(jj = j-1; jj <= j+1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = ddC[ii+1]*C[jj-j+1];
          ijp++;
        }
      }
      else if(i == nx+1)
      {
        for(ii = nx-1; ii <= nx+1; ii++) for(jj = j-1; jj <= j+1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = ddC[ii-nx+1]*C[jj-j+1];
          ijp++;
        }
      }
      else if(j == -1)
      {
        for(ii = i-1; ii <= i+1; ii++) for(jj = -1; jj <= +1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = C[ii-i+1]*ddC[jj+1];
          ijp++;
        }
      }
      else if(j == ny+1)
      {
        for(ii = i-1; ii <= i+1; ii++) for(jj = ny-1; jj <= ny+1; jj++)
        {
          (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
          (*Ax)[ijp] = C[ii-i+1]*ddC[jj-ny+1];
          ijp++;
        }
      }
      else
      {
        printf("you should not reach to this point...\n");
        exit(1);
      }
    }
    else
    {
      for(ii = i-1; ii <= i+1; ii++) for(jj = j-1; jj <= j+1; jj++)
      {
        (*Ai)[ijp] = idx_2d_CS(ii,jj,nx,ny);
        (*Ax)[ijp] = C[ii-i+1]*C[jj-j+1];
        ijp++;
      }
    }
  }

  assert(ijp == dim*9);

  umfpack_di_symbolic(dim,dim,*Ap,*Ai,*Ax,&sym,NULL,NULL);
  umfpack_di_numeric(*Ap,*Ai,*Ax,sym,Num,NULL,NULL);
  umfpack_di_free_symbolic(&sym);
}


void init_intp_2d_CS_umfpack(double x1, double x2, int nx, 
    double y1, double y2, int ny,
    double *data)
{
  double *b, *Ax;
  size_t dim;
  int i, j, ij, *Ap, *Ai;
  void *Num;

  dim = (nx+3)*(ny+3);

  b = (double*)malloc(sizeof(double)*dim);

  constr_umfpackM_2d_CS(nx, ny, &Ap, &Ai, &Ax, &Num);

  // set components of source vector ------------------------------------
  for(i = -1; i <= nx+1; i++) for(j = -1; j <= ny+1; j++)
  {
    ij = idx_2d_CS(i, j, nx, ny);

    if(i == -1 || i == nx+1 || j == -1 || j == ny+1) b[ij] = 0.0;
    else b[ij] = data[ij];
  }
  // set components of source vector ------------------------------------

  umfpack_di_solve(UMFPACK_At, Ap, Ai, Ax, data, b, Num, NULL, NULL);

  umfpack_di_free_numeric(&Num);
  free(Ap);
  free(Ai);
  free(Ax);
  free(b);
}

/* ========================================================================= */

void constr_umfpackM_xth_CS(int nx, int nt, 
    int **Ap, int **Ai, double **Ax, void **Num)
{
  int i, j, ij, ip, jp, jp_in, cnt;
  size_t dim;
  void *sym;
  double  C[3] = {1.0/6.0, 4.0/6.0, 1.0/6.0},
        ddC[3] = {1.0, -2.0, 1.0}, Cj;

  dim = (nx+3)*nt;

  *Ap = (int*)malloc(sizeof(int)*(dim+1));
  *Ai = (int*)malloc(sizeof(int)*(dim*9));
  *Ax = (double*)malloc(sizeof(double)*(dim*9));

  (*Ap)[0] = 0;
  cnt = 0;
  for(i = -1; i <= nx+1; i++) for(j = 0; j < nt; j++)
  {
    ij = idx_xth_CS(i,j,nx,nt);
    (*Ap)[ij+1] = (*Ap)[ij]+9;

    if(i >= 0 && i <= nx)
    {
      for(ip = i-1; ip <= i+1; ip++) for(jp = 0; jp < nt; jp++)
      {
        jp_in = 0;

        if(j == 0)         { if(jp <= 1 || jp == nt-1)  jp_in = 1; }
        else if(j == nt-1) { if(jp >= nt-2 || jp == 0)  jp_in = 1; }
        else               { if(jp >= j-1 && jp <= j+1) jp_in = 1; }

        if(jp_in == 1)
        {
          if(j == jp) Cj = 4.0/6.0;
          else Cj = 1.0/6.0;

          (*Ai)[cnt] = idx_xth_CS(ip,jp,nx,nt);
          (*Ax)[cnt] = C[ip-i+1]*Cj;
          cnt++;
        }
      }
    }
    else if(i == -1)
    {
      for(ip = -1; ip <= 1; ip++) for(jp = 0; jp < nt; jp++)
      {
        jp_in = 0;

        if(j == 0)         { if(jp <= 1 || jp == nt-1)  jp_in = 1; }
        else if(j == nt-1) { if(jp >= nt-2 || jp == 0)  jp_in = 1; }
        else               { if(jp >= j-1 && jp <= j+1) jp_in = 1; }

        if(jp_in == 1)
        {
          if(j == jp) Cj = 4.0/6.0;
          else Cj = 1.0/6.0;

          (*Ai)[cnt] = idx_xth_CS(ip,jp,nx,nt);
          (*Ax)[cnt] = ddC[ip+1]*Cj;
          cnt++;
        }
      }
    }
    else /* i = nx+1 */
    {
      for(ip = nx-1; ip <= nx+1; ip++) for(jp = 0; jp < nt; jp++)
      {
        jp_in = 0;

        if(j == 0)         { if(jp <= 1 || jp == nt-1)  jp_in = 1; }
        else if(j == nt-1) { if(jp >= nt-2 || jp == 0)  jp_in = 1; }
        else               { if(jp >= j-1 && jp <= j+1) jp_in = 1; }

        if(jp_in == 1)
        {
          if(j == jp) Cj = 4.0/6.0;
          else Cj = 1.0/6.0;

          (*Ai)[cnt] = idx_xth_CS(ip,jp,nx,nt);
          (*Ax)[cnt] = ddC[ip-nx+1]*Cj;
          cnt++;
        }
      }
    }
  }

  assert(cnt == dim*9);

  umfpack_di_symbolic(dim,dim,*Ap,*Ai,*Ax,&sym,NULL,NULL);
  umfpack_di_numeric(*Ap,*Ai,*Ax,sym,Num,NULL,NULL);
  umfpack_di_free_symbolic(&sym);
}


void init_intp_xth_CS_umfpack(double x1, double x2, int nx, 
    double y1, double y2, int nt,
    double *data)
{
  double *b, *Ax;
  size_t dim;
  int i, j, ij, *Ap, *Ai;
  void *Num;

  dim = (nx+3)*nt;

  b = (double*)malloc(sizeof(double)*dim);

  constr_umfpackM_xth_CS(nx, nt, &Ap, &Ai, &Ax, &Num);

  // set components of source vector ------------------------------------
  for(i = -1; i <= nx+1; i++) for(j = 0; j < nt; j++)
  {
    ij = idx_xth_CS(i, j, nx, nt);

    if(i == -1 || i == nx+1) b[ij] = 0.0;
    else b[ij] = data[ij];
  }
  // set components of source vector ------------------------------------

  umfpack_di_solve(UMFPACK_At, Ap, Ai, Ax, data, b, Num, NULL, NULL);

  umfpack_di_free_numeric(&Num);
  free(Ap);
  free(Ai);
  free(Ax);
  free(b);
}
/* ========================================================================= */





double intp_xth_CS(double x, double y, 
    double x1, double x2, int nx, 
    double y1, double y2, int nt,
    double *array)
{
  double tfi, tfj, h, t, wi[4], wj[4], val, dx, dy;
  int tmpi, tmpj, i[4], j[4], ip, jp, ij;

  dx = (x2 - x1)/nx;
  dy = (y2 - y1)/nt;

  tfi = (x - x1)/dx;
  tmpi = (int)tfi; h = tfi - tmpi; t = 1.0 - h;
  i[0] = tmpi - 1;
  i[1] = tmpi;
  i[2] = tmpi + 1;
  i[3] = tmpi + 2;
  wi[0] = t*t*t;
  wi[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wi[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wi[3] = h*h*h;

  tfj = (y - y1)/dy;
  tmpj = (int)tfj; h = tfj - tmpj; t = 1.0 - h;
  j[0] = tmpj - 1;
  j[1] = tmpj;
  j[2] = tmpj + 1;
  j[3] = tmpj + 2;
  wj[0] = t*t*t;
  wj[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wj[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wj[3] = h*h*h;

  /* ensure the radial array boundary */
  for(ip = 0; ip <= 3; ip++) i[ip] = min(max(i[ip], -1), nx+1);

  /* make theta boundary periodic */
  for(jp = 0; jp <= 3; jp++)
  {
    if(j[jp] >= nt) j[jp] = j[jp]%nt;
    else if(j[jp] < 0) j[jp] = j[jp]%nt + nt;
    assert(j[jp] >= 0 && j[jp] < nt);
  }

  val = 0.0;
  for(jp = 0; jp <= 3; jp++) for(ip = 0; ip <= 3; ip++)
  {
    ij = idx_xth_CS(i[ip],j[jp],nx,nt);
    val += wi[ip]*wj[jp]*array[ij];
  }

  return val/36.0;
}


void intp_dfn_xth_CS(double x, double y, 
    double x1, double x2, int nx,
    double y1, double y2, int nt,
    double *array, double *val)
{
  double tfi, tfj, CR[4], CZ[4], dCR[4], dCZ[4], ddCR[4], ddCZ[4],
  dx, dy, h, t, tt, tt1, tt2;
  int tmpi, tmpj, i[4], j[4], ip, jp, ij;


  dx = (x2 - x1)/nx;
  dy = (y2 - y1)/nt;

  /* ------------------------------------------------------------- */
  tfi = (x - x1)/dx;
  tmpi = (int)tfi; i[0] = tmpi - 1;
  h = tfi - tmpi; t = 1.0 - h;

  CR[0] = t*t*t;
  CR[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  CR[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  CR[3] = h*h*h;

  for(ip = 0; ip <= 3; ip++)
  {
    i[ip] = i[0] + ip;
    tt = tfi - i[ip];

    tt1 = absv(tt + 0.5); tt2 = absv(tt - 0.5);
    dCR[ip] = (tsc(tt1) - tsc(tt2))/dx;

    tt1 = absv(tt + 1.0); tt2 = absv(tt - 1.0); tt = absv(tt);
    ddCR[ip] = (cic(tt1) + cic(tt2) - 2.0*cic(tt))/(dx*dx);
  }
  /* ------------------------------------------------------------- */


  /* ------------------------------------------------------------- */
  tfj = (y - y1)/dy;
  tmpj = (int)tfj; j[0] = tmpj - 1;
  h = tfj - tmpj; t = 1.0 - h;

  CZ[0] = t*t*t;
  CZ[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  CZ[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  CZ[3] = h*h*h;

  for(jp = 0; jp <= 3; jp++)
  {
    j[jp] = j[0] + jp;;
    tt = tfj - j[jp];

    tt1 = absv(tt + 0.5); tt2 = absv(tt - 0.5);
    dCZ[jp] = (tsc(tt1) - tsc(tt2))/dy;

    tt1 = absv(tt + 1.0); tt2 = absv(tt - 1.0); tt = absv(tt);
    ddCZ[jp] = (cic(tt1) + cic(tt2) - 2.0*cic(tt))/(dy*dy);
  }
  /* ------------------------------------------------------------- */

  /* ensure the radial array boundary */
  for(ip = 0; ip <= 3; ip++) 
  {
    i[ip] = min(max(i[ip], -1), nx+1);
    CR[ip] /= 6.0;
  }

  /* make theta boundary periodic */
  for(jp = 0; jp <= 3; jp++)
  {
    if(j[jp] >= nt) j[jp] = j[jp]%nt;
    else if(j[jp] < 0) j[jp] = j[jp]%nt + nt;
    assert(j[jp] >= 0 && j[jp] < nt);

    CZ[jp] /= 6.0;
  }

  for(ip = 0; ip <= 5; ip++) val[ip] = 0.0;

  for(jp = 0; jp <= 3; jp++) for(ip = 0; ip <= 3; ip++)
  {
    ij = idx_xth_CS(i[ip], j[jp], nx, nt);

    /* interpolated psi */
    val[0] += CR[ip]*CZ[jp]*array[ij];

    /* 1st derivatives */
    val[1] += dCR[ip]*CZ[jp]*array[ij];
    val[2] += CR[ip]*dCZ[jp]*array[ij];

    /* 2nd derivatives */
    val[3] += ddCR[ip]*CZ[jp]*array[ij];
    val[4] += CR[ip]*ddCZ[jp]*array[ij];
    val[5] += dCR[ip]*dCZ[jp]*array[ij];
  }
}



double intp_1d_CS(double x, 
		double x1, double x2, int nx, 
		double *array)
{
  double tfi, h, t, wi[4], val=0.0;
  int tmpi, i[4], ip;

  tfi = (x - x1)*nx/(x2 - x1);
  tmpi = (int)tfi; h = tfi - tmpi; t = 1.0 - h;
  i[0] = tmpi - 1;
  i[1] = tmpi;
  i[2] = tmpi + 1;
  i[3] = tmpi + 2;
  wi[0] = t*t*t;
  wi[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wi[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wi[3] = h*h*h;

  /* ensure the array boundary */
  for(ip = 0; ip <= 3; ip++) i[ip] = min(max(i[ip], -1), nx+1);

  for(ip = 0; ip <= 3; ip++) val += wi[ip]*array[idx_1d_CS(i[ip])];

  return val/6.0;
}

void intp_dfn_1d_CS(double x,
		    double x1, double x2, int nx,
		    double *array, double *val)
{
  double dx, tfi, CR[4], dCR[4], ddCR[4], h, t, tt, tt1, tt2;
  int tmpi, i[4], ip, ic;

  //assert(x >= x1 && x <= x2);

  dx = (x2 - x1)/nx;

  /* ------------------------------------------------------------- */
  tfi = (x - x1)/dx;
  tmpi = (int)tfi; i[0] = tmpi - 1;
  h = tfi - tmpi; t = 1.0 - h;

  CR[0] = t*t*t;
  CR[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  CR[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  CR[3] = h*h*h;

  for(ip = 0; ip <= 3; ip++)
  {
    i[ip] = i[0] + ip;
    tt = tfi - i[ip];

    tt1 = absv(tt + 0.5); tt2 = absv(tt - 0.5);
    dCR[ip] = (tsc(tt1) - tsc(tt2))/dx;

    tt1 = absv(tt + 1.0); tt2 = absv(tt - 1.0); tt = absv(tt);
    ddCR[ip] = (cic(tt1) + cic(tt2) - 2.0*cic(tt))/(dx*dx);
  }
  /* ------------------------------------------------------------- */


  /* ensure the array boundary */
  for(ip = 0; ip <= 3; ip++) i[ip] = min(max(i[ip], -1), nx+1);

  val[0] = 0.0;
  val[1] = 0.0;
  val[2] = 0.0;

  for(ip = 0; ip <= 3; ip++)
  {
    ic = idx_1d_CS(i[ip]);
	
    /* interpolated radial */
    val[0] += CR[ip]*array[ic];

    /* 1st derivatives */
    val[1] += dCR[ip]*array[ic];

    /* 2nd derivatives */
    val[2] += ddCR[ip]*array[ic];
  }

  val[0] /= 6.0;
}




double intp_2d_CS(double x, double y, 
    double x1, double x2, int nx, 
    double y1, double y2, int ny, 
    double *array)
{
  double tfi, tfj, h, t, wi[4], wj[4], val, dx, dy;
  int tmpi, tmpj, i[4], j[4], ip, jp, ij;

  dx = (x2 - x1)/nx;
  dy = (y2 - y1)/ny;

  tfi = (x - x1)/dx;
  tmpi = (int)tfi; h = tfi - tmpi; t = 1.0 - h;
  i[0] = tmpi - 1;
  i[1] = tmpi;
  i[2] = tmpi + 1;
  i[3] = tmpi + 2;
  wi[0] = t*t*t;
  wi[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wi[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wi[3] = h*h*h;

  tfj = (y - y1)/dy;
  tmpj = (int)tfj; h = tfj - tmpj; t = 1.0 - h;
  j[0] = tmpj - 1;
  j[1] = tmpj;
  j[2] = tmpj + 1;
  j[3] = tmpj + 2;
  wj[0] = t*t*t;
  wj[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wj[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wj[3] = h*h*h;

  /* ensure the array boundary */
  for(ip = 0; ip <= 3; ip++) 
  {
    i[ip] = min(max(i[ip], -1), nx+1);
    j[ip] = min(max(j[ip], -1), ny+1);
  }

  val = 0.0;
  for(jp = 0; jp <= 3; jp++) for(ip = 0; ip <= 3; ip++)
  {
    //ij = j[jp]+1 + (ny+3)*(i[ip]+1);
    ij = idx_2d_CS(i[ip],j[jp],nx,ny);
    val += wi[ip]*wj[jp]*array[ij];
  }

  return val/36.0;
}

void intp2_2d_CS(double x, double y, 
    double x1, double x2, int nx, 
    double y1, double y2, int ny, 
    double *array1, double *array2, double *val1, double *val2)
{
  double tfi, tfj, h, t, wi[4], wj[4], dx, dy;
  int tmpi, tmpj, i[4], j[4], ip, jp, ij;

  dx = (x2 - x1)/nx;
  dy = (y2 - y1)/ny;

  tfi = (x - x1)/dx;
  tmpi = (int)tfi; h = tfi - tmpi; t = 1.0 - h;
  i[0] = tmpi - 1;
  i[1] = tmpi;
  i[2] = tmpi + 1;
  i[3] = tmpi + 2;
  wi[0] = t*t*t;
  wi[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wi[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wi[3] = h*h*h;

  tfj = (y - y1)/dy;
  tmpj = (int)tfj; h = tfj - tmpj; t = 1.0 - h;
  j[0] = tmpj - 1;
  j[1] = tmpj;
  j[2] = tmpj + 1;
  j[3] = tmpj + 2;
  wj[0] = t*t*t;
  wj[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  wj[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  wj[3] = h*h*h;

  /* ensure the array boundary */
  for(ip = 0; ip <= 3; ip++) 
  {
    i[ip] = min(max(i[ip], -1), nx+1);
    j[ip] = min(max(j[ip], -1), ny+1);
  }

  *val1 = 0.0;
  *val2 = 0.0;
  for(jp = 0; jp <= 3; jp++) for(ip = 0; ip <= 3; ip++)
  {
    //ij = j[jp]+1 + (ny+3)*(i[ip]+1);
    ij = idx_2d_CS(i[ip],j[jp],nx,ny);
    *val1 += wi[ip]*wj[jp]*array1[ij];
    *val2 += wi[ip]*wj[jp]*array2[ij];
  }

  *val1 /= 36.0;
  *val2 /= 36.0;
}

void intp_dfn_2d_CS(double x, double y, 
    double x1, double x2, int nx,
    double y1, double y2, int ny,
    double *array, double *val)
{
  double tfi, tfj, CR[4], CZ[4], dCR[4], dCZ[4], ddCR[4], ddCZ[4],
  dx, dy, h, t, tt, tt1, tt2;
  int tmpi, tmpj, i[4], j[4], ip, jp, ij;


  dx = (x2 - x1)/nx;
  dy = (y2 - y1)/ny;

  /* ------------------------------------------------------------- */
  tfi = (x - x1)/dx;
  tmpi = (int)tfi; i[0] = tmpi - 1;
  h = tfi - tmpi; t = 1.0 - h;

  CR[0] = t*t*t;
  CR[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  CR[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  CR[3] = h*h*h;

  for(ip = 0; ip <= 3; ip++)
  {
    i[ip] = i[0] + ip;
    tt = tfi - i[ip];

    tt1 = absv(tt + 0.5); tt2 = absv(tt - 0.5);
    dCR[ip] = (tsc(tt1) - tsc(tt2))/dx;

    tt1 = absv(tt + 1.0); tt2 = absv(tt - 1.0); tt = absv(tt);
    ddCR[ip] = (cic(tt1) + cic(tt2) - 2.0*cic(tt))/(dx*dx);
  }
  /* ------------------------------------------------------------- */


  /* ------------------------------------------------------------- */
  tfj = (y - y1)/dy;
  tmpj = (int)tfj; j[0] = tmpj - 1;
  h = tfj - tmpj; t = 1.0 - h;

  CZ[0] = t*t*t;
  CZ[1] = 1.0 + 3.0*(t + t*t - t*t*t);
  CZ[2] = 1.0 + 3.0*(h + h*h - h*h*h);
  CZ[3] = h*h*h;

  for(jp = 0; jp <= 3; jp++)
  {
    j[jp] = j[0] + jp;;
    tt = tfj - j[jp];

    tt1 = absv(tt + 0.5); tt2 = absv(tt - 0.5);
    dCZ[jp] = (tsc(tt1) - tsc(tt2))/dy;

    tt1 = absv(tt + 1.0); tt2 = absv(tt - 1.0); tt = absv(tt);
    ddCZ[jp] = (cic(tt1) + cic(tt2) - 2.0*cic(tt))/(dy*dy);
  }
  /* ------------------------------------------------------------- */


  /* ensure the array boundary */
  for(ip = 0; ip <= 3; ip++)
  {
    i[ip] = min(max(i[ip], -1), nx+1);
    j[ip] = min(max(j[ip], -1), ny+1);

    CR[ip] /= 6.0; CZ[ip] /= 6.0;
  }

  for(ip = 0; ip <= 5; ip++) val[ip] = 0.0;

  for(jp = 0; jp <= 3; jp++) for(ip = 0; ip <= 3; ip++)
  {
    //ij = j[jp]+1 + (ny+3)*(i[ip]+1);
    ij = idx_2d_CS(i[ip], j[jp], nx, ny);

    /* interpolated psi */
    val[0] += CR[ip]*CZ[jp]*array[ij];

    /* 1st derivatives */
    val[1] += dCR[ip]*CZ[jp]*array[ij];
    val[2] += CR[ip]*dCZ[jp]*array[ij];

    /* 2nd derivatives */
    val[3] += ddCR[ip]*CZ[jp]*array[ij];
    val[4] += CR[ip]*ddCZ[jp]*array[ij];
    val[5] += dCR[ip]*dCZ[jp]*array[ij];
  }
}

/* (nr, nz) is the actual dimension of the matrix data */
void init_intp_2d_CS(double x1, double x2, int nx,
    double y1, double y2, int ny,
    double *data, int bc_x, int bc_y)
{
  int i, j, l, nxy, in, ip, jn, jp;
  double *s3, q0_inv, q1_inv, q2_inv,
       e_val, w_val, n_val, s_val, ne_val, nw_val, se_val, sw_val;

  q0_inv = 6.0*6.0/(4.0*4.0);
  q1_inv = 0.25;
  q2_inv = 0.25*0.25;

  nxy = (nx+3)*(ny+3);
  s3 = (double*)malloc(sizeof(double)*nxy);
  assert(s3 != NULL);

  /*  +++++++++++++++++++
      +#################+
      +#################+
      +#################+
      +#################+
      +#################+
      +#################+
      +++++++++++++++++++ */

  /*  +#################+
      +#################+
      +#################+
      +#################+
      +#################+
      +#################+
      +#################+
      +#################+ */

  /*  ###################
###################
###################
###################
###################
###################
###################
################### */

  for(i = 0; i <= nx; i++)
    for(j = 0; j <= ny; j++) 
      s3[idx_2d_CS(i,j,nx,ny)] = data[idx_2d_CS(i,j,nx,ny)]; // initial guess

  if(bc_y == 0) // periodic
  {
    for(i = 0; i <= nx; i++)
    {
      s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
      s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
    }
  }
  else if(bc_y == 1) // vanishing 1st derivative
  {
    for(i = 0; i <= nx; i++)
    {
      s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
      s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
    }
  }
  else if(bc_y == 2) // vanishing 2nd derivative
  {
    for(i = 0; i <= nx; i++)
    {
      s3[idx_2d_CS(i,-1, nx,ny)] = 
        2.0*s3[idx_2d_CS(i,0, nx,ny)] - s3[idx_2d_CS(i,1, nx,ny)];
      s3[idx_2d_CS(i,ny+1, nx,ny)] = 
        2.0*s3[idx_2d_CS(i,ny, nx,ny)] - s3[idx_2d_CS(i,ny-1, nx,ny)];
    }
  }

  if(bc_x == 0) // periodic
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
    }
  }
  else if(bc_x == 1) // vanishing 1st derivative
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
  }
  else if(bc_x == 2) // vanishing 2nd derivative
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = 
        2.0*s3[idx_2d_CS(0,j, nx,ny)] - s3[idx_2d_CS(1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = 
        2.0*s3[idx_2d_CS(nx,j, nx,ny)] - s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
  }

  /* iterate until desired accuracy is met */
  for(l = 0; l < N_CS_ITER; l++)
  {
    for(i = 0; i <= nx; i++)
    {
      ip = i+1; in = i-1;
      for(j = 0; j <= ny; j++)
      {
        jp = j+1; jn = j-1;

        e_val = s3[idx_2d_CS(ip,j, nx,ny)];
        w_val = s3[idx_2d_CS(in,j, nx,ny)];

        n_val = s3[idx_2d_CS(i,jp, nx,ny)];
        ne_val = s3[idx_2d_CS(ip,jp, nx,ny)];
        nw_val = s3[idx_2d_CS(in,jp, nx,ny)];

        s_val = s3[idx_2d_CS(i,jn, nx,ny)];
        se_val = s3[idx_2d_CS(ip,jn, nx,ny)];
        sw_val = s3[idx_2d_CS(in,jn, nx,ny)];

        s3[idx_2d_CS(i,j, nx,ny)] = data[idx_2d_CS(i,j, nx,ny)]*q0_inv
          - (n_val + s_val + e_val + w_val)*q1_inv
          - (ne_val + nw_val + se_val + sw_val)*q2_inv;
      }
      //}

      if(bc_y == 0) // periodic
      {
        for(i = 0; i <= nx; i++)
        {
          s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
          s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
        }
      }
      else if(bc_y == 1) // vanishing 1st derivative
      {
        for(i = 0; i <= nx; i++)
        {
          s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
          s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
        }
      }
      else if(bc_y == 2) // vanishing 2nd derivative
      {
        for(i = 0; i <= nx; i++)
        {
          s3[idx_2d_CS(i,-1, nx,ny)] = 
            2.0*s3[idx_2d_CS(i,0, nx,ny)] - s3[idx_2d_CS(i,1, nx,ny)];
          s3[idx_2d_CS(i,ny+1, nx,ny)] = 
            2.0*s3[idx_2d_CS(i,ny, nx,ny)] - s3[idx_2d_CS(i,ny-1, nx,ny)];
        }
      }
  }

  if(bc_x == 0) // periodic
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
    }
  }
  else if(bc_x == 1) // vanishing 1st derivative
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
  }
  else if(bc_x == 2) // vanishing 2nd derivative
  {
    for(j = -1; j <= ny+1; j++)
    {
      s3[idx_2d_CS(-1,j, nx,ny)] = 
        2.0*s3[idx_2d_CS(0,j, nx,ny)] - s3[idx_2d_CS(1,j, nx,ny)];
      s3[idx_2d_CS(nx+1,j, nx,ny)] = 
        2.0*s3[idx_2d_CS(nx,j, nx,ny)] - s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
  }

  /*if(bc == 0) // periodic
    {
    for(i = 0; i <= nx; i++)
    {
    s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
    s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
    }
    for(j = -1; j <= ny+1; j++)
    {
    s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
    s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
    }
    }
    else if(bc == 1) // vanishing 1st derivative
    {
    for(i = 0; i <= nx; i++)
    {
    s3[idx_2d_CS(i,-1, nx,ny)] = s3[idx_2d_CS(i,1, nx,ny)];
    s3[idx_2d_CS(i,ny+1, nx,ny)] = s3[idx_2d_CS(i,ny-1, nx,ny)];
    }
    for(j = -1; j <= ny+1; j++)
    {
    s3[idx_2d_CS(-1,j, nx,ny)] = s3[idx_2d_CS(1,j, nx,ny)];
    s3[idx_2d_CS(nx+1,j, nx,ny)] = s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
    }
    else if(bc == 2) // vanishing 2nd derivative
    {
    for(i = 0; i <= nx; i++)
    {
    s3[idx_2d_CS(i,-1, nx,ny)] = 
    2.0*s3[idx_2d_CS(i,0, nx,ny)] - s3[idx_2d_CS(i,1, nx,ny)];
    s3[idx_2d_CS(i,ny+1, nx,ny)] = 
    2.0*s3[idx_2d_CS(i,ny, nx,ny)] - s3[idx_2d_CS(i,ny-1, nx,ny)];
    }
    for(j = -1; j <= ny+1; j++)
    {
    s3[idx_2d_CS(-1,j, nx,ny)] = 
    2.0*s3[idx_2d_CS(0,j, nx,ny)] - s3[idx_2d_CS(1,j, nx,ny)];
    s3[idx_2d_CS(nx+1,j, nx,ny)] = 
    2.0*s3[idx_2d_CS(nx,j, nx,ny)] - s3[idx_2d_CS(nx-1,j, nx,ny)];
    }
    }*/

  /*for(i = 0; i <= nx; i++)
    {
    ip = i+1; in = i-1;
    for(j = 0; j <= ny; j++)
    {
    jp = j+1; jn = j-1;

    e_val = s3[idx_2d_CS(ip,j, nx,ny)];
    w_val = s3[idx_2d_CS(in,j, nx,ny)];

    n_val = s3[idx_2d_CS(i,jp, nx,ny)];
    ne_val = s3[idx_2d_CS(ip,jp, nx,ny)];
    nw_val = s3[idx_2d_CS(in,jp, nx,ny)];

    s_val = s3[idx_2d_CS(i,jn, nx,ny)];
    se_val = s3[idx_2d_CS(ip,jn, nx,ny)];
    sw_val = s3[idx_2d_CS(in,jn, nx,ny)];

    s3[idx_2d_CS(i,j, nx,ny)] = data[idx_2d_CS(i,j, nx,ny)]*q0_inv
    - (n_val + s_val + e_val + w_val)*q1_inv
    - (ne_val + nw_val + se_val + sw_val)*q2_inv;
    }
    }*/
}

for(l = 0; l < nxy; l++) data[l] = s3[l];
free(s3);
s3 = NULL;
}
