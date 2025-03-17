#include "cnst-srt.h"
#include <math.h>

#define idx_xfth_CS(i,j) ((j)+1 + (N_FTH+3)*((i)+1))
#define idx_xth_CS(i,j,nx,nt) ((j) + (nt)*((i)+1))

void farr_to_carr(int *intpar, double *dblpar)//,int *nxx,int *nth)
{

  g_inpv.nxx = intpar[0];
  g_inpv.nth = intpar[1];
  g_inpv.ifd0 = intpar[2];
  g_inpv.ifd1 = intpar[3];
  g_inpv.op_geo_info = intpar[4]; 
  g_pall.ipe = intpar[5]; // youknow

  g_inpv.xx0 = dblpar[0] ;
  g_inpv.xx1 = dblpar[1];
  g_inpv.cir_q0 = dblpar[2];
  g_inpv.cir_q1 = dblpar[3];
  g_inpv.cir_q2 = dblpar[4];
  g_inpv.cir_q3 = dblpar[5];
  g_inpv.cir_q4 = dblpar[6];
  g_inpv.cir_a0 = dblpar[7];
  g_inpv.cir_R0 = dblpar[8];
  g_inpv.cir_B0 = dblpar[9];
  g_inpv.ase_q0 = dblpar[10];
  g_inpv.ase_a0 = dblpar[11];
  g_inpv.ase_R0 = dblpar[12];
  g_inpv.ase_B0 = dblpar[13];
  g_inpv.ase_kappa = dblpar[14];
  g_inpv.ase_del = dblpar[15];
  g_inpv.ase_shift = dblpar[16];



  g_inpv.rrl_cir=1.1;

  /* for ase or circle data,
 * this parameter indicates the position of circular limiter
 * in the unit of the minor radius.
 *
  g_eqd.Rmin = mR0*(1.0 - g_ase.shift) - rrm;
  g_eqd.Rmax = mR0*(1.0 - g_ase.shift) + rrm;
  g_eqd.Zmin = -g_ase.kappa*rrm;
  g_eqd.Zmax =  g_ase.kappa*rrm;
 */

  g_inpv.op_RZ_xfth_file=0;
  // 1: read RZ_xfth.dat file for g_xth.epm_CS, R_CS, Z_CS etc..
  g_inpv.op_test_output=0; 
  // 1: write "test-lim1.txt" file for g_RZ.lim_R, lim_Z
  
    if (g_pall.ipe == 0) {
        printf("Generating output of eqmodule\n");
    } else if (g_pall.ipe == 1) {
        printf("Not Generating output of eqmodule\n");
    } else {
        printf("Invalid value\n");
    }

  //  g_pall.ipe=1; 
  //  g_pall.ipe=0; 

  // 0: write test-ase1.txt, test-lim1.txt, xth-test3 4.txt, 
  // RZ_xfth.dat, RZ_xfth_test.txt" and some explanation for geqdsk


//  printf("inpv.:%d)\n",g_inpv.nxx);
//  printf("inpv.:%d)\n",g_inpv.nth);
//  printf("inpv.:%f)\n",g_inpv.ase_del);
//  printf("inpv.:%f)\n",g_inpv.ase_shift);
//  printf("pall.ipe:%d)\n",g_pall.ipe);
}


void init_diag(double *ptr, int length)
{
  int ij, i, j,k, nwr;
  double ep, th, dtdR, dtdZ, pv[6], pv2[6], pv3[6], R, Z;
  double ucf_len, ucf_flux, ucf_I, ucf_B;
  double *R_arr,*Z_arr,*th_arr,*B_arr,*Bp_arr,*I_arr,*J_arr,*gradpsi2,*gradth2,*gradpsith;
  double *psi_arr,*dpdr_arr, *dpdz_arr, *dpdrr_arr, *dpdzz_arr, *dpdrz_arr, *temp_arr, *temp2_arr;
  double *dBdp_arr, *dBdt_arr;
  char fname[CHARLEN];
  FILE *fp;
  
  R_arr      = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  Z_arr      = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  th_arr     = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  psi_arr    = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dpdr_arr   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dpdz_arr   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dpdrr_arr  = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dpdzz_arr  = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1)); 
  dpdrz_arr  = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  B_arr      = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  Bp_arr     = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  I_arr      = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  J_arr      = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  gradpsi2   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  gradth2    = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  gradpsith  = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dBdp_arr   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  dBdt_arr   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  temp_arr   = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  temp2_arr  = (double*)malloc(sizeof(double)*(g_xth.nxx+1)*(g_xth.nth+1));
  

  ucf_len = 1.0e2;
  ucf_flux = 1.0e4*1.0e2*1.0e2;
  ucf_I = 1.0e4*1.0e2;
  ucf_B = 1.0e4;

  // g_eqd.Rc 

  for(i=0; i<=g_xth.nxx; i++)
  {
    for(j=0; j<=g_xth.nth-1; j++)
    {
      ij = j+(g_xth.nth)*(i);
      
      th=j*g_xth.dth;
      th_arr[ij]=mod_th(th);
      intp_RZ_xfth(g_xth.xx[i],th_arr[ij],&R,&Z);
//     R_arr[ij]=g_xth.R_CS[(i+1)*(N_FTH+3)+j];
//     Z_arr[ij]=g_xth.Z_CS[(i+1)*(N_FTH+3)+j];
//     R_arr[ij]=g_xth.eqR_CS[ij];
//     Z_arr[ij]=g_xth.eqZ_CS[ij];
      R_arr[ij]=R;
      Z_arr[ij]=Z;
      intp_dpsi_RZ(R_arr[ij],Z_arr[ij],pv);
     // printf("R,Z,xx,th:%f, %f, %f, %f\n",R_arr[ij],Z_arr[ij],g_xth.xx[i],th_arr[ij]);
      psi_arr[ij]=pv[0];
      dpdr_arr[ij]=pv[1];
      dpdz_arr[ij]=pv[2];
      dpdrr_arr[ij]=pv[3];
      dpdzz_arr[ij]=pv[4];
      dpdrz_arr[ij]=pv[5];
      ep=double_sqrt((R_arr[ij]-1.0)*(R_arr[ij]-1.0)+Z_arr[ij]*Z_arr[ij]);
      dtdR=-double_sin(th)/(ep+1.0e-20);
      dtdZ=double_cos(th)/(ep+1.0e-20);
      Bp_arr[ij]=double_sqrt(pv[1]*pv[1]+pv[2]*pv[2])/R_arr[ij];
      temp_arr[ij]=(pv[1]*pv[1]+pv[2]*pv[2])/R_arr[ij];
      J_arr[ij]=R_arr[ij]/(pv[2]*dtdR-pv[1]*dtdZ); // 1 / B dot grad theta
      gradpsi2[ij]=pv[1]*pv[1]+pv[2]*pv[2];
      gradth2[ij]=dtdR*dtdR+dtdZ*dtdZ;
      gradpsith[ij]=pv[1]*dtdR+pv[2]*dtdZ;
      intp_dI_RZ(R_arr[ij],Z_arr[ij],pv2);
      I_arr[ij]=pv2[0];
      B_arr[ij]=double_sqrt(pv[1]*pv[1]+pv[2]*pv[2]+pv2[0]*pv2[0])/R_arr[ij];
      intp_dfn_2d_CS(g_xth.xx[i], th_arr[ij], g_xth.xx0, g_xth.xx1, g_xth.nxx,
           0.0, g_c.twopi, N_FTH, g_xth.Bmod_CS, pv);
      //printf("xx,xx0,xx1,th,0,2pi: %g, %g, %g, %g, %g, %g\n",g_xth.xx[i],g_xth.xx0,g_xth.xx1,th_arr[ij],0.0,g_c.twopi);
      temp_arr[ij]=pv[0];
      temp2_arr[ij]=0;
      dBdp_arr[ij]=pv[1]/2.0/g_RZ.psi_last/g_xth.xx[i];
      dBdt_arr[ij]=pv[2];


      // eqdsk value
      R_arr[ij]=R_arr[ij]*g_eqd.Rc/ucf_len;
      Z_arr[ij]=Z_arr[ij]*g_eqd.Rc/ucf_len;
      psi_arr[ij]=(psi_arr[ij]*g_eqd.Rc*g_eqd.Rc*g_eqd.Bc+g_eqd.psic)/ucf_flux+g_eqd.simag;
      dpdr_arr[ij]=dpdr_arr[ij]*g_eqd.Rc*g_eqd.Bc/ucf_len/ucf_B;
      dpdz_arr[ij]=dpdz_arr[ij]*g_eqd.Rc*g_eqd.Bc/ucf_len/ucf_B;
      dpdrr_arr[ij]=dpdrr_arr[ij]*g_eqd.Bc/ucf_B;
      dpdzz_arr[ij]=dpdzz_arr[ij]*g_eqd.Bc/ucf_B;
      dpdrz_arr[ij]=dpdrz_arr[ij]*g_eqd.Bc/ucf_B;
      B_arr[ij]=B_arr[ij]*g_eqd.Bc/ucf_B;
      Bp_arr[ij]=Bp_arr[ij]*g_eqd.Bc/ucf_B;
      I_arr[ij]=I_arr[ij]*g_eqd.Ic/ucf_I;
      J_arr[ij]=J_arr[ij]/g_eqd.Bc*g_eqd.Rc*ucf_B/ucf_len;
      gradpsi2[ij]=gradpsi2[ij]*g_eqd.Rc*g_eqd.Rc*g_eqd.Bc*g_eqd.Bc/ucf_I/ucf_I;
      gradpsith[ij]=gradpsith[ij]*g_eqd.Bc/ucf_B;
      gradth2[ij]=gradth2[ij]/g_eqd.Rc/g_eqd.Rc*ucf_len*ucf_len;
      dBdp_arr[ij]=dBdp_arr[ij]/g_eqd.Rc/g_eqd.Rc*ucf_len*ucf_len;
      dBdt_arr[ij]=dBdt_arr[ij]*g_eqd.Bc/ucf_B;
      temp_arr[ij]=temp_arr[ij]*g_eqd.Bc/ucf_B;
      temp2_arr[ij]=temp2_arr[ij]*g_eqd.Bc/ucf_B;

    }  
  }


    nwr=(g_xth.nxx+1)*(g_xth.nth);

    for(j=0;j<nwr;j++){
     ptr[j]=R_arr[j];
     ptr[j+1*nwr]=Z_arr[j];
     ptr[j+2*nwr]= th_arr[j];
     ptr[j+3*nwr]= psi_arr[j];
     ptr[j+4*nwr]= dpdr_arr[j];
     ptr[j+5*nwr]= dpdz_arr[j];
     ptr[j+6*nwr]= dpdrr_arr[j];
     ptr[j+7*nwr]= dpdzz_arr[j];
     ptr[j+8*nwr]= dpdrz_arr[j];
     ptr[j+9*nwr]= B_arr[j];
     ptr[j+10*nwr]= Bp_arr[j];
     ptr[j+11*nwr]= I_arr[j];
     ptr[j+12*nwr]= J_arr[j];
     ptr[j+13*nwr]= gradpsi2[j];
     ptr[j+14*nwr]= gradth2[j];
     ptr[j+15*nwr]= gradpsith[j];
     ptr[j+16*nwr]= dBdp_arr[j];
     ptr[j+17*nwr]= dBdt_arr[j];
     ptr[j+18*nwr]= temp_arr[j];
     ptr[j+19*nwr]= temp2_arr[j];
    } 



  if (g_pall.ipe==0) {
//  printf("\nI_arr\n");  
//  for(i = 0 ;i < nwr ; i++) {
//    printf("%lf ",I_arr[i]);  
//  }

// "/home/sjhsong/Codes/eq_util/2/arr_I.txt

  printf("print geo_arr_RZ.txt\n");

  sprintf(fname, "geo_arr_RZ.txt", g_sys.work_dir);
  fp = fopen(fname, "w"); assert(fp != NULL);
  for(j = 0; j < nwr ; j++){
    k=j/g_xth.nth;
//    fprintf(fp, "%e %e %e %e %g %i %e %e %e %e %e \n",R_arr[j],Z_arr[j],th_arr[j],psi_arr[j],g_xth.xx[k],k,dpdr_arr[j],dpdz_arr[j],dpdrr_arr[j],dpdzz_arr[j],dpdrz_arr[j]);
    fprintf(fp, "%e %e %e %e %e %e %e %e %e \n",R_arr[j],Z_arr[j],th_arr[j],psi_arr[j],dpdr_arr[j],dpdz_arr[j],dpdrr_arr[j],dpdzz_arr[j],dpdrz_arr[j]);
  }
  fflush(fp);

  printf("print geo_arr_B.txt\n");
  sprintf(fname, "geo_arr_B.txt", g_sys.work_dir);
  fp = fopen(fname, "w");
  for(j = 0; j < nwr ; j++) {
    fprintf(fp, "%g %g %g %g %g %g %g %g %g %g %g\n",B_arr[j],Bp_arr[j],I_arr[j],J_arr[j],gradpsi2[j],gradth2[j],gradpsith[j],dBdp_arr[j],dBdt_arr[j],temp_arr[j],temp2_arr[j]);
  }
  fflush(fp);

  printf("print geo_arr_rho.txt\n");
  sprintf(fname, "geo_arr_rho.txt", g_sys.work_dir);
  fp = fopen(fname, "w");
  for(j = 0; j < sizeof(g_xth.xx)+1 ; j++) {
    fprintf(fp, "%g %i \n",g_xth.xx[j],j);
  }
  fflush(fp);

  printf("g_eqd.Rc: %g, g_eqd.Bc: %g, g_RZ.psi_last: %g\n",g_eqd.Rc,g_eqd.Bc,g_RZ.psi_last);
  printf("g_eqd.psic: %g, g_eqd.psi0: %g, g_eqd.psi1: %g\n",g_eqd.psic,g_eqd.psi0,g_eqd.psi1);
  printf("psi_arr[0]: %g, psi_arr[end]: %g\n",psi_arr[0],psi_arr[nwr-1]);
  }
//
//  free(R_arr);
//  free(Z_arr);
//  free(th_arr);
//  free(psi_arr);
//  free(dpdr_arr);
//  free(dpdz_arr);
//  free(dpdrr_arr);
//  free(dpdzz_arr);
//  free(dpdrz_arr);
//  free(B_arr);
//  free(Bp_arr);
//  free(I_arr);
//  free(J_arr);
//  free(gradpsi2_arr);
//  free(gradth2_arr);
//  free(gradpsith_arr)

}
