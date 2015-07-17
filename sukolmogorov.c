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
	int nw;			/* number of frequency samples per trace */
  int nfft;   /* For computing the FFT */
	int i;			/* time sample index */
  int logar;   /* Check if the trace has had a log applied (for mie scattering) */
  int specsq; /* Check if the input is the spectrum or spectrum squared */
  const double eps = 1.e-32;
  complex *fft     = NULL;
  kiss_fftr_cfg forw, invs;

	/* hook up getpar */
	initargs(argc, argv);
	requestdoc(1);
  //requestdoc(0); // For now (testing), stdin is not used
  
  if (!gettr(&tr)) err("can't get first trace");
  nw = tr.ns;
  nfft = 2*(nw-1);

  if (!getparint("log",&logar)) logar = 0;
  if (!getparint("specsq",&specsq)) specsq = 0;

  /* Allocating memory */
  fft       = alloc1complex(nfft/2+1);
  forw      = kiss_fftr_alloc(nfft,0,NULL,NULL);
  invs      = kiss_fftr_alloc(nfft,1,NULL,NULL);
  if(NULL == forw || NULL == invs)
    err("KISS FFT allocation error");

  /* Squaring the spectrum */
  if(!specsq) {
    for(i = 0; i < nw; i++) {
      tr.data[i] *= tr.data[i];
    }
  }

  /* Take the log and create a complex type */
  if(!logar) {
    for(i = 0; i < nw; i++) {
      fft[i] = cmplx(log(tr.data[i]+eps)/nfft,0.);  
    }
  }
  else {
    for(i = 0; i < nw; i++) {
      fft[i] = cmplx(tr.data[i]/nfft,0.); 
    }
  }

  /* Find the inverse FFT */
  kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data); //change to pfarc

  tr.data[0]      *= 0.5;
  tr.data[nfft/2] *= 0.5;
  for(i=1+nfft/2; i < nfft; i++) {
    tr.data[i] = 0.;
  }

  kiss_fftr(forw, tr.data, (kiss_fft_cpx *) fft);

  for(int i=0; i < nw; i++) {
    fft[i] = crmul(cwp_cexp(fft[i]),1./nfft);
  }

  // Put in cosine taper window

  kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data);

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
