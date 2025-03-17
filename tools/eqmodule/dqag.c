/* dqag.c -- modified version of QUADPACK routine DQAG.
 * (C)1999, C. Bond. All right reserved.
 *
 * There are no changes to the basic computational method. Only
 * the temporary storage strategy is changed to utilize the
 * local stack at the appropriate level. This reduces the
 * need for memory allocation of arrays at higher levels and
 * the resulting passing of memory pointers down the line.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "cnst-srt.h"

/* DQAG - Approximation to definite integral. (From QUADPACK)
 *
 *  Calls DQAGE with appropriate parameters assigned.
 *
 * PARAMETERS:
 *
 *	f() - double precision function to be integrated.
 *
 *	a - lower limit of integration.
 *
 *	b - upper limit of integration.
 *
 *	epsabs - absolute accuracy requested.
 *
 *	epsrel - relative accuracy requested.
 *
 *	irule - integration rule to be used as follows:
 *		irule = 1 -- G_K 7-15
 *		irule = 2 -- G_K 10-21
 *		irule = 3 -- G_K 15-31
 *		irule = 4 -- G_K 20-41
 *		irule = 5 -- G_K 25-51
 *		irule = 6 -- G_K 30-61
 */ 	
double dqag(double f(),double a,double b,double epsabs,
    double epsrel,int irule,double *abserr,int *neval,int *ier)
{
  double result;
  int last;

  result = dqage(f,a,b,epsabs,epsrel,irule,abserr,neval,ier,&last); 

  return result;
}


	
/* DQAGE - Approximation to definite integral. (From QUADPACK)
 *
 *	Allows user's choice of Gauss-Kronrod integration rule.
 *
 * PARAMETERS:
 *
 *	f() - double precision function to be integrated.
 *
 *	a - lower limit of integration.
 *
 *	b - upper limit of integration.
 *
 *	epsabs - absolute accuracy requested.
 *
 *	epsrel - relative accuracy requested.
 *
 *	irule - integration rule to be used as follows:
 *		irule = 1 -- G_K 7-15
 *		irule = 2 -- G_K 10-21
 *		irule = 3 -- G_K 15-31
 *		irule = 4 -- G_K 20-41
 *		irule = 5 -- G_K 25-51
 *		irule = 6 -- G_K 30-61
 *
 *	limit - maximum number of subintervals.
 */ 	
double dqage(double f(),double a,double b,double epsabs,double epsrel,
    int irule,double *abserr,int *neval,int *ier,int *last)
{
  double area,area1=0.0,area2=0.0,area12,a1,a2,b1,b2,c,defabs;
  double defab1,defab2,errbnd,errmax,error1,error2;
  double erro12,errsum,resabs,result;
  double alist[LIMIT],blist[LIMIT],rlist[LIMIT],elist[LIMIT];
  int iroff1,iroff2,k,keyf,maxerr,nrmax,iord[LIMIT],limit;

  limit = LIMIT - 1;
  *ier = 0;
  *neval = 0;
  *last = 0;
  result = 0.0;
  *abserr = 0.0;
  alist[0] = a;
  blist[0] = b;
  rlist[0] = 0.0;
  elist[0] = 0.0;
  iord[0] = 0;
  defabs = 0.0;
  resabs = 0.0;
  if ((epsabs < 0.0) && (epsrel < 0.0))
    *ier = 6;
  if (*ier == 6) return result; 

  /* First approximation to the integral. */
  keyf = irule;
  if (irule <= 0) keyf = 1;
  if (irule >= 7) keyf = 6;
  c = keyf;
  *neval = 0;
  switch (keyf) {
    case 1: 
      result = G_K15(f,a,b,abserr,&defabs,&resabs);
      break;
    case 2:
      result = G_K21(f,a,b,abserr,&defabs,&resabs);
      break;
    case 3:
      result = G_K31(f,a,b,abserr,&defabs,&resabs);
      break;
    case 4:
      result = G_K41(f,a,b,abserr,&defabs,&resabs);
      break;
    case 5:
      result = G_K51(f,a,b,abserr,&defabs,&resabs);
      break;
    case 6:
      result = G_K61(f,a,b,abserr,&defabs,&resabs);
      break;
  }
  *last = 0;
  rlist[0] = result;
  elist[0] = *abserr;
  iord[0] = 0;

  /* Test on accuracy. */
  errbnd = max(epsabs,epsrel * fabs(result));
  if ((*abserr <= 50.0 * epmach * defabs) && (*abserr > errbnd))
    *ier = 2;
  if (limit == 0) *ier = 1;
  if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
      (*abserr == 0.0)) goto _60;

  /* Initialization. */
  errmax = *abserr;
  maxerr = 0;
  area = result;
  errsum = *abserr;
  nrmax = 0;
  iroff1 = 0;
  iroff2 = 0;

  /* Main Loop. */
  for (*last = 1; *last <= limit; (*last)++) {
    /* Bisect the subinterval with the largest error estimate. */
    a1 = alist[maxerr];
    b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
    a2 = b1;
    b2 = blist[maxerr];
    switch (keyf) {
      case 1:
        area1 = G_K15(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K15(f,a2,b2,&error2,&resabs,&defab2);
        break;
      case 2:
        area1 = G_K21(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K21(f,a2,b2,&error2,&resabs,&defab2);
        break;
      case 3:
        area1 = G_K31(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K31(f,a2,b2,&error2,&resabs,&defab2);
        break;
      case 4:
        area1 = G_K41(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K41(f,a2,b2,&error2,&resabs,&defab2);
        break;
      case 5:
        area1 = G_K51(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K51(f,a2,b2,&error2,&resabs,&defab2);
        break;
      case 6:
        area1 = G_K61(f,a1,b1,&error1,&resabs,&defab1);
        area2 = G_K61(f,a2,b2,&error2,&resabs,&defab2);
        break;
    }

    /* Improve previous approximations to integral and error,
       and test for accuracy. */
    (*neval) += 1;
    area12 = area1 + area2;
    erro12 = error1 + error2;
    errsum = errsum + erro12 - errmax;
    area = area + area12 - rlist[maxerr];
    if ((defab1 != error1) && (defab2 != error2)) {
      if ((fabs(rlist[maxerr]-area12) <= 1.0e-5 * fabs(area12)) &&
          (erro12 >= .99 * errmax)) 
        iroff1++;
      if ((*last > 9) && (erro12 > errmax)) 
        iroff2++;
    }
    rlist[maxerr] = area1;
    rlist[*last] = area2;
    errbnd = max(epsabs,epsrel * fabs(area));
    if (errsum > errbnd)  {

      /* Test for roundoff error and eventually set error flag. */
      if ((iroff1 > 6) || (iroff2 > 20))  
        *ier = 2;

      /* Set error flag in the case that the number of subintervals
         equals the limit. */
      if (*last == limit)
        *ier = 1;

      /* Set error flag in the case of bad integrand behavior at a
         point of the integration range. */
      if (max(fabs(a1),fabs(b2)) <= (1.0 + c * 1000.0 * epmach) *
          (fabs(a2)+1.0e4 * uflow))
        *ier = 3;
    }
    /* Append the newly-created intervals to the list. */

    if (error2 <= error1) {
      alist[*last] = a2;
      blist[maxerr] = b1;
      blist[*last] = b2;
      elist[maxerr] = error1;
      elist[*last] = error2;
    }
    else {
      alist[maxerr] = a2;
      alist[*last] = a1;
      blist[*last] = b1;
      rlist[maxerr] = area2;
      rlist[*last] = area1;
      elist[maxerr] = error2;
      elist[*last] = error1;
    }

    /* Call DQSORT to maintain the descending ordering in the list of
       error estimates and select the subinterval with the
       largest error estimate (to be bisected next). */

    dqsort(limit,*last,&maxerr,&errmax,elist,iord,&nrmax);
    if ((*ier != 0) || (errsum <= errbnd)) break;
  }

  /* Compute final result. */

  result = 0.0;
  for (k = 0; k <= *last; k++) {
    result += rlist[k];
  }
  *abserr = errsum;
_60:
  if (keyf != 1)
    *neval = (10 * keyf + 1) * (2 * (*neval) + 1);
  else
    *neval = 30 * (*neval) + 15;

  return result;
}	



void dqsort(int limit,int last,int *maxerr,double *ermax,double elist[],
    int iord[],int *nrmax)
{
  double errmax,errmin;
  int i,ibeg,ido,isucc,j,jbnd,jupbn,k;

  if (last > 1) goto _10;
  iord[0] = 0;
  iord[1] = 1;
  goto _90;
_10:
  errmax = elist[*maxerr];
  if (*nrmax == 0) goto _30;
  ido = (*nrmax) - 1;
  for (i = 0;i <= ido; i++) {
    isucc = iord[*nrmax-1];
    if (errmax <= elist[isucc]) goto _30;
    iord[*nrmax] = isucc;
    (*nrmax)--;
  }
_30:
  jupbn = last;
  if (last > (limit/2 + 2))
    jupbn = limit + 3 - last;
  errmin = elist[last];
  jbnd = jupbn - 1;
  ibeg = *nrmax + 1;
  if (ibeg > jbnd) goto _50;
  for (i = ibeg; i <= jbnd; i++) {
    isucc = iord[i];
    if (errmax >= elist[isucc]) goto _60;
    iord[i-1] = isucc;
  }
_50: 
  iord[jbnd] = *maxerr;
  iord[jupbn] = last;
  goto _90;
_60:
  iord[i-1] = *maxerr;
  k = jbnd;
  for (j = i;j <= jbnd; j++) {
    isucc = iord[k];
    if (errmin < elist[isucc]) goto _80;
    iord[k+1] = isucc;
    k--;
  }
  iord[i] = last;
  goto _90;
_80:
  iord[k+1] = last;
_90:
  *maxerr = iord[*nrmax];
  *ermax = elist[*maxerr];
  return;  
}	

double G_K15(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  static long double XGK15[8] = {
    0.99145537112081263921,
    0.94910791234275852453,
    0.86486442335976907279,
    0.74153118559939443986,
    0.58608723546769113029,
    0.40584515137739716691,
    0.20778495500789846760,
    0.00000000000000000000};
  static long double WGK15[8] = {
    0.02293532201052922496,
    0.06309209262997855329,
    0.10479001032225018384,
    0.14065325971552591875,
    0.16900472663926790283,
    0.19035057806478540991,
    0.20443294007529889241,
    0.20948214108472782801};
  static long double WG7[4] = {
    0.12948496616886969327,
    0.27970539148927666790,
    0.38183005050511894495,
    0.41795918367346938776};
  double fv1[7],fv2[7];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  fc=(*f)(centr);
  resg = fc * WG7[3];
  resk = fc * WGK15[7];
  *resabs = fabs(resk);
  for (j = 0; j < 3; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK15[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG7[j] * fsum;
    resk += WGK15[jtw] * fsum;
    *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 4; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK15[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk = resk + WGK15[jtwm1] * fsum;
    *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK15[7] * fabs(fc - reskh);
  for (j = 0; j < 7; j++ )
    *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}

double G_K21(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  static long double XGK21[11] = {
    0.99565716302580808074,
    0.97390652851717172008,
    0.93015749135570822600,
    0.86506336668898451073,
    0.78081772658641689706,
    0.67940956829902440623,
    0.56275713466860468334,
    0.43339539412924719080,
    0.29439286270146019813,
    0.14887433898163121088,
    0.00000000000000000000};
  static long double WGK21[11] = {
    0.01169463886737187428,
    0.03255816230796472748,
    0.05475589657435199603,
    0.07503967481091995277,
    0.09312545458369760554,
    0.10938715880229764190,
    0.12349197626206585108,
    0.13470921731147332593,
    0.14277593857706008080,
    0.14773910490133849137,
    0.14944555400291690566};
  static long double WG10[5] = {
    0.06667134430868813759,
    0.14945134915058059315,
    0.21908636251598204400,
    0.26926671930999635509,
    0.29552422471475287017};
  double fv1[10],fv2[10];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  resg = 0.0;
  fc=(*f)(centr);
  resk = fc * WGK21[10];
  *resabs = fabs(resk);
  for (j = 0; j < 5; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK21[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG10[j] * fsum;
    resk += WGK21[jtw] * fsum;
    *resabs = *resabs + WGK21[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 5; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK21[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk = resk + WGK21[jtwm1] * fsum;
    *resabs = (*resabs) + WGK21[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK21[10] * fabs(fc - reskh);
  for (j = 0; j < 10; j++ )
    *resasc = (*resasc) + WGK21[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}

double G_K31(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  static long double XGK31[16] = {
    0.99800229869339706029,
    0.98799251802048542849,
    0.96773907567913913426,
    0.93727339240070590431,
    0.89726453234408190088,
    0.84820658341042721620,
    0.79041850144246593297,
    0.72441773136017004742,
    0.65099674129741697053,
    0.57097217260853884754,
    0.48508186364023968069,
    0.39415134707756336990,
    0.29918000715316881217,
    0.20119409399743452230,
    0.10114206691871749903,
    0.00000000000000000000};
  static long double WGK31[16] = {
    0.00537747987292334899,
    0.01500794732931612254,
    0.02546084732671532019,
    0.03534636079137584622,
    0.04458975132476487661,
    0.05348152469092808727,
    0.06200956780067064029,
    0.06985412131872825871,
    0.07684968075772037889,
    0.08308050282313302104,
    0.08856444305621177065,
    0.09312659817082532123,
    0.09664272698362367851,
    0.09917359872179195933,
    0.10076984552387559504,
    0.10133000701479154902};
  static long double WG15[8] = {
    0.03075324199611726835,
    0.07036604748810812471,
    0.10715922046717193501,
    0.13957067792615431445,
    0.16626920581699393355,
    0.18616100001556221103,
    0.19843148532711157646,
    0.20257824192556127288};

  double fv1[15],fv2[15];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  fc=(*f)(centr);
  resg = fc * WG15[7];
  resk = fc * WGK31[15];
  *resabs = fabs(resk);
  for (j = 0; j < 7; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK31[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG15[j] * fsum;
    resk += WGK31[jtw] * fsum;
    *resabs = *resabs + WGK31[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 8; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK31[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk = resk + WGK31[jtwm1] * fsum;
    *resabs = (*resabs) + WGK31[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK31[15] * fabs(fc - reskh);
  for (j = 0; j < 15; j++ )
    *resasc = (*resasc) + WGK31[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}

double G_K41(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  /* Gauss-Kronrod abscissae and weights for 41 - 20 rules */
  static long double XGK41[21] = {
    0.99885903158827766384,
    0.99312859918509492479,
    0.98150787745025025919,
    0.96397192727791379127,
    0.94082263383175475352,
    0.91223442825132590587,
    0.87827681125228197608,
    0.83911697182221882339,
    0.79504142883755119835,
    0.74633190646015079261,
    0.69323765633475138481,
    0.63605368072651502545,
    0.57514044681971031534,
    0.51086700195082709800,
    0.44359317523872510320,
    0.37370608871541956067,
    0.30162786811491300432,
    0.22778585114164507808,
    0.15260546524092267551,
    0.07652652113349733375,
    0.00000000000000000000};
  static long double WGK41[21] = {
    0.00307358371852053150,
    0.00860026985564294220,
    0.01462616925697125298,
    0.02038837346126652360,
    0.02588213360495115883,
    0.03128730677703279896,
    0.03660016975820079803,
    0.04166887332797368626,
    0.04643482186749767472,
    0.05094457392372869193,
    0.05519510534828599474,
    0.05911140088063957237,
    0.06265323755478116803,
    0.06583459713361842211,
    0.06864867292852161935,
    0.07105442355344406831,
    0.07303069033278666750,
    0.07458287540049918899,
    0.07570449768455667466,
    0.07637786767208073671,
    0.07660071191799965645};
  static long double WG20[10] = {
    0.01761400713915211831,
    0.04060142980038694133,
    0.06267204833410906357,
    0.08327674157670474872,
    0.10193011981724043504,
    0.11819453196151841731,
    0.13168863844917662690,
    0.14209610931838205133,
    0.14917298647260374679,
    0.15275338713072585070};	
  double fv1[20],fv2[20];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  resg = 0.0;
  fc=(*f)(centr);
  resk = fc * WGK41[20];
  *resabs = fabs(resk);
  for (j = 0; j < 10; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK41[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG20[j] * fsum;
    resk += WGK41[jtw] * fsum;
    *resabs = *resabs + WGK41[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 10; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK41[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk = resk + WGK41[jtwm1] * fsum;
    *resabs = (*resabs) + WGK41[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK41[20] * fabs(fc - reskh);
  for (j = 0; j < 20; j++ )
    *resasc = (*resasc) + WGK41[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}

double G_K51(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  /* Gauss-Kronrod abscissae and weights for 51 - 25 rules */
  static long double XGK51[26] = {
    0.99926210499260983419,
    0.99555696979049809791,
    0.98803579453407724764,
    0.97666392145951751150,
    0.96161498642584251242,
    0.94297457122897433941,
    0.92074711528170156175,
    0.89499199787827536885,
    0.86584706529327559545,
    0.83344262876083400142,
    0.79787379799850005941,
    0.75925926303735763058,
    0.71776640681308438819,
    0.67356636847346836449,
    0.62681009901031741279,
    0.57766293024122296772,
    0.52632528433471918260,
    0.47300273144571496052,
    0.41788538219303774885,
    0.36117230580938783774,
    0.30308953893110783017,
    0.24386688372098843205,
    0.18371893942104889202,
    0.12286469261071039639,
    0.06154448300568507889,
    0.00000000000000000000};
  static long double WGK51[26] = {
    0.00198738389233031593,
    0.00556193213535671376,
    0.00947397338617415161,
    0.01323622919557167481,
    0.01684781770912829823,
    0.02043537114588283546,
    0.02400994560695321622,
    0.02747531758785173780,
    0.03079230016738748889,
    0.03400213027432933784,
    0.03711627148341554356,
    0.04008382550403238207,
    0.04287284502017004948,
    0.04550291304992178891,
    0.04798253713883671391,
    0.05027767908071567196,
    0.05236288580640747586,
    0.05425112988854549014,
    0.05595081122041231731,
    0.05743711636156783285,
    0.05868968002239420796,
    0.05972034032417405998,
    0.06053945537604586295,
    0.06112850971705304831,
    0.06147118987142531666,
    0.06158081806783293508};
  static long double WG25[13] = {
    0.01139379850102628795,
    0.02635498661503213726,
    0.04093915670130631266,
    0.05490469597583519193,
    0.06803833381235691721,
    0.08014070033500101801,
    0.09102826198296364981,
    0.10053594906705064420,
    0.10851962447426365312,
    0.11485825914571164834,
    0.11945576353578477223,
    0.12224244299031004169,
    0.12317605372671545120};

  double fv1[25],fv2[25];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  fc=(*f)(centr);
  resg = fc * WG25[12];
  resk = fc * WGK51[25];
  *resabs = fabs(resk);
  for (j = 0; j < 12; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK51[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG25[j] * fsum;
    resk += WGK51[jtw] * fsum;
    *resabs = *resabs + WGK51[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 13; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK51[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk = resk + WGK51[jtwm1] * fsum;
    *resabs = (*resabs) + WGK51[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK51[25] * fabs(fc - reskh);
  for (j = 0; j < 25; j++ )
    *resasc = (*resasc) + WGK51[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}

double G_K61(double f(),double a,double b,double *abserr,
    double *resabs,double *resasc)
{
  /* Gauss-Kronrod abscissae and weights for 61 - 30 rules */

  static long double XGK61[31] = {
    0.99948441005049063757,
    0.99689348407464954027,
    0.99163099687040459486,
    0.98366812327974720997,
    0.97311632250112626837,
    0.96002186496830751222,
    0.94437444474855997942,
    0.92620004742927432588,
    0.90557330769990779855,
    0.88256053579205268154,
    0.85720523354606109896,
    0.82956576238276839744,
    0.79972783582183908301,
    0.76777743210482619492,
    0.73379006245322680473,
    0.69785049479331579693,
    0.66006106412662696137,
    0.62052618298924286114,
    0.57934523582636169176,
    0.53662414814201989926,
    0.49248046786177857499,
    0.44703376953808917678,
    0.40040125483039439254,
    0.35270472553087811347,
    0.30407320227362507737,
    0.25463692616788984644,
    0.20452511668230989144,
    0.15386991360858354696,
    0.10280693796673703015,
    0.05147184255531769583,
    0.00000000000000000000};
  static long double WGK61[31] = {
    0.00138901369867700762,
    0.00389046112709988405,
    0.00663070391593129217,
    0.00927327965951776343,
    0.01182301525349634174,
    0.01436972950704580481,
    0.01692088918905327263,
    0.01941414119394238117,
    0.02182803582160919230,
    0.02419116207808060137,
    0.02650995488233310161,
    0.02875404876504129284,
    0.03090725756238776247,
    0.03298144705748372603,
    0.03497933802806002414,
    0.03688236465182122922,
    0.03867894562472759295,
    0.04037453895153595911,
    0.04196981021516424615,
    0.04345253970135606932,
    0.04481480013316266319,
    0.04605923827100698812,
    0.04718554656929915395,
    0.04818586175708712914,
    0.04905543455502977889,
    0.04979568342707420636,
    0.05040592140278234684,
    0.05088179589874960649,
    0.05122154784925877217,
    0.05142612853745902593,
    0.05149472942945156756};
  static long double WG30[15] = {
    0.00796819249616660562,
    0.01846646831109095914,
    0.02878470788332336935,
    0.03879919256962704960,
    0.04840267283059405290,
    0.05749315621761906648,
    0.06597422988218049513,
    0.07375597473770520627,
    0.08075589522942021535,
    0.08689978720108297980,
    0.09212252223778612872,
    0.09636873717464425964,
    0.09959342058679526706,
    0.10176238974840550460,
    0.10285265289355884034};

  double fv1[30],fv2[30];
  double absc,centr,dhlgth;
  double fc,fsum,fval1,fval2,hlgth;
  double resg,resk,reskh,result;
  int j,jtw,jtwm1;

  centr = 0.5 * (a + b);
  hlgth = 0.5 * (b - a);
  dhlgth = fabs(hlgth);

  resg = 0.0;
  fc=(*f)(centr);
  resk = fc * WGK61[30];
  *resabs = fabs(resk);
  for (j = 0; j < 15; j++) {
    jtw = 2 * j + 1;
    absc = hlgth * XGK61[jtw];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    fsum = fval1 + fval2;
    resg += WG30[j] * fsum;
    resk += WGK61[jtw] * fsum;
    *resabs = *resabs + WGK61[jtw] * (fabs(fval1) + fabs(fval2));
  }
  for (j = 0; j < 15; j++) {
    jtwm1 = j * 2;
    absc = hlgth * XGK61[jtwm1];
    fval1 = (*f)(centr-absc);
    fval2 = (*f)(centr+absc);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    fsum = fval1 + fval2;
    resk += WGK61[jtwm1] * fsum;
    *resabs = *resabs + WGK61[jtwm1] * (fabs(fval1) + fabs(fval2));
  }
  reskh = resk * 0.5;
  *resasc = WGK61[30] * fabs(fc - reskh);
  for (j = 0; j < 30; j++ )
    *resasc = (*resasc) + WGK61[j] * (fabs(fv1[j] - reskh) +
        fabs(fv2[j] - reskh));
  result = resk * hlgth;
  *resabs = (*resabs) * dhlgth;
  *resasc = (*resasc) * dhlgth;
  *abserr = fabs((resk - resg) * hlgth);
  if ((*resasc != 0.0) && (*abserr != 0.0))
    *abserr = (*resasc) * min(1.0,pow((200.0 * (*abserr)/(*resasc)),1.5));
  if (*resabs > uflow/(50.0 * epmach))
    *abserr = max(epmach * 50.0 * (*resabs),(*abserr)); 	
  return result;
}
