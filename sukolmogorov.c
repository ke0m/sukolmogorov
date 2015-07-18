/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKOLMOGOROV: $Revision: 1.0 $ ; $Date: 2015/07/17 19:45:15 $		*/
 
#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUKOLMOGOROV - The kolmogorov method for finding the time domain minimum phase equivalent of a spectrum",
"									     ",
"  sukolmogorov <stdin >stdout [optional parameters]				     ",
"									     ",
" Optional Parameters:							     ",
" log=0              =1 to take in account that input is already the log of spectrum  ",
" specsq=0           =1 to take in account that input is already the spectrum sqaured ",
" Notes: Forthcoming...								     ",
"									     ",
NULL};

/* Credits:
 *	SEP: 
 *	CWP: 
 *	TCCS: 
 * Technical Reference:
 *
 * Trace header fields accessed: ns, dt
 */
/**************** end self doc *******************************************/

#define LOOKFAC 2 /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft           */

segy tr;

int
main(int argc, char **argv)
{
	int nf;			/* number of frequency samples per trace */
  int nfft;   /* For computing the FFT */
	int i;			/* time sample index */
  int logar;   /* Check if the trace has had a log applied (for mie scattering) */
  int specsq; /* Check if the input is the spectrum or spectrum squared */
  const double eps = 1.e-32;
  complex *ct      = NULL;
  float   *rt      = NULL;

	/* hook up getpar */
	initargs(argc, argv);
	requestdoc(1);
  
  /* Read in the trace */
  if (!gettr(&tr)) err("can't get first trace");
  nf = 2*tr.ns;
  nfft = npfaro(nf,LOOKFAC*nf);

  if (!getparint("log",&logar)) logar = 0;
  if (!getparint("specsq",&specsq)) specsq = 0;

  /* Allocating memory */
  ct = alloc1complex(nf);

  fprintf(stderr, "input spec: \n");
  for(i = 0; i < nf; i++) {
    fprintf(stderr, "i=%d input=%f\n",i,tr.data[i]);
  }
  fprintf(stderr, "\n");

  /* Squaring the spectrum */
  if(!specsq) {
    for(i = 0; i < nf; i++) {
      tr.data[i] *= tr.data[i];
    }
  }

  /* Take the log and create a complex type */
  if(!logar) {
    for(i = 0; i < nf; i++) { //This needs to be the length of the data
      ct[i].r = log(tr.data[i]+eps)/nfft; // Therefore, the nfft needs to be twice that length
      ct[i].i = 0.f;
    }
  }
  else {
    for(i = 0; i < nf; i++) {
      ct[i].r = tr.data[i];
      ct[i].i = 0.f;
    }
  }

  fprintf(stderr, "log of spec(nfft=%d): \n",nfft);
  for(i = 0; i < nf; i++) {
    fprintf(stderr, "i=%d real=%f imag=%f\n",i,ct[i].r,ct[i].i);
  }
  fprintf(stderr, "\n");

  /* Find the inverse FFT */
  fprintf(stderr, "the ifft: \n");
  pfacr(-1, nfft, ct, tr.data);
  for(i = 0; i < nfft; i++) {
    fprintf(stderr, "i=%d ifft=%f\n",i,tr.data[i]);
  }
  fprintf(stderr, "\n");

  tr.data[0]      *= 0.5;
  tr.data[nfft/2] *= 0.5;
  for(i=1+nfft/2; i < nfft; i++) {
    tr.data[i] = 0.;
  }

  pfarc(1,nfft,tr.data,ct);

  for(int i=0; i < nf; i++) {
    ct[i] = crmul(cwp_cexp(ct[i]),1./nfft);
  }

  // Put in cosine taper window

  fprintf(stderr, "output: \n");
  pfacr(-1,nfft,ct,tr.data);
  for(i = 0; i < nfft; i++) {
    fprintf(stderr, "i=%d output=%f\n",i,tr.data[i]);
  }
 
  tr.ns = nfft;
  tr.f1 = 0.f;

  puttr(&tr);
	return(CWP_Exit());
}

