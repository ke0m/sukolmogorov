/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKOLMOGOROV: $Revision: 0.1 $ ; $Date: 2016/07/12 19:45:15 $		*/
 
#include "su.h"
#include "segy.h"
#include "kiss_fft.h"
#include "kiss_fftr.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUKOLMOGOROV - The kolmogorov method for finding the time domain minimum phase equivalent of a spectrum",
"									     ",
"  sukolmogorov <stdin >stdout [optional parameters]				     ",
"									     ",
" Optional Parameters:							     ",
" log=0           =1 to take in account that input is already log of spectrum  ",
" specsq=0           =1 to take in account that input is already log of spectrum  ",
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

segy tr;
void create_spectrum(int n, float o1, float d1, float* spectrum);

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
  kiss_fftr_cfg forw, invs;

	/* hook up getpar */
	initargs(argc, argv);
	requestdoc(1);
  //requestdoc(0); // For now (testing), stdin is not used
  
  /* Read in the trace */
  if (!gettr(&tr)) err("can't get first trace");
  nf = tr.ns;
  //nfft = 2*(nf - 1); // I am not sure if this should be -2 or -1
  nfft = npfar(nf);

  if (!getparint("log",&logar)) logar = 0;
  if (!getparint("specsq",&specsq)) specsq = 0;

  /* Allocating memory */
  ct       = alloc1complex(nf);
  rt       = ealloc1float(nfft);
  
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
  // Note that I am no longer dividing by nfft. Therefore, the outputs will be different
  if(!logar) {
    for(i = 0; i < nf; i++) { //This needs to be the length of the data
      ct[i].r = log(tr.data[i]+eps)/nfft; // Therefore, the nfft needs to be twice that length
      ct[i].i = 0.f;
      fprintf(stderr, "i=%d fftr=%f ffti=%f\n", i, ct[i].r, ct[i].i);
    }
  }
  else {
    for(i = 0; i < nf; i++) {
      ct[i].r = tr.data[i];
      ct[i].i = 0.f;
    }
  }

  fprintf(stderr, "\n");
  fprintf(stderr, "The ifft: nfft=%d \n", nfft);
  /* Find the inverse FFT */
  //kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data); //change to pfarc
  pfacr(-1, nfft, ct, tr.data);
  for(i = 0; i < nfft; i++) {
    fprintf(stderr, "i=%d rt=%f\n", i, tr.data[i]);
  }
  fprintf(stderr, "\n");

  tr.data[0]      *= 0.5;
  tr.data[nfft/2] *= 0.5;
  for(i=1+nfft/2; i < nfft; i++) {
    tr.data[i] = 0.;
  }

  fprintf(stderr, "The ffft: nfft=%d \n", nfft);
  //kiss_fftr(forw, tr.data, (kiss_fft_cpx *) fft);
  pfarc(1,nfft,tr.data,ct);
  for(i = 0; i < nf; i++) { //This needs to be the length of the data
    fprintf(stderr, "i=%d fftr=%f ffti=%f\n", i, ct[i].r, ct[i].i);
  }
  fprintf(stderr, "\n");

  for(int i=0; i < nf; i++) {
    ct[i] = crmul(cwp_cexp(ct[i]),1./nfft);
  }

  // Put in cosine taper window

  //kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data);
  pfacr(-1,nfft,ct,tr.data);
  //for(i = 0; i < nfft; i++) {
  //  fprintf(stderr, "i=%d rt=%f\n", i, tr.data[i]);
  //}
 
  tr.ns = nfft;
  tr.f1 = 0.f;


  puttr(&tr);
	return(CWP_Exit());
}

/* Functions from Jons GEE example */
void create_spectrum(int n, float o1, float d1, float* spectrum) {

  float x = 0.f;
  for(int i = 0; i < n+1; i++) {
    x = o1 + d1*i;
    spectrum[i] = exp(-fabs(6*(x)));
    //fprintf(stderr, "i=%d spectrum=%f\n", i,spectrum[i]);
  }
}
