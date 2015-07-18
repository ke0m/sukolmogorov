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
	//int nt;			/* number of frequency samples per trace */
  int nfft;   /* For computing the FFT */
	//float dw;		/* time sampling interval */
	//int i1;			/* time sample index */
  //int log;   /* Check if the trace has had a log applied (for mie scattering) */
  int j;     /* Other indices */
  const double eps = 1.e-32;
  kiss_fftr_cfg forw, invs;

  /* Variables for Jons example */
  float *spectrum  = NULL;
  float *rspectrum = NULL;
  complex *fft     = NULL;
  int n            = 64;
  int nw           = n/2+1;
  float o1         = -M_PI/2;
  float d1         = M_PI/n;

  nfft = n;

	/* hook up getpar */
	initargs(argc, argv);
	//requestdoc(1);
  requestdoc(0); // For now (testing), stdin is not used

  /* Allocating memory */
  spectrum  = ealloc1float((n+1));
  rspectrum = ealloc1float((n/2+1));
  fft       = alloc1complex(nfft/2+1);
  forw      = kiss_fftr_alloc(nfft,0,NULL,NULL);
  invs      = kiss_fftr_alloc(nfft,1,NULL,NULL);
  if(NULL == forw || NULL == invs)
    err("KISS FFT allocation error");
  memset( (void *) tr.data, 0, (n+1) * FSIZE);

  tr.dt = d1;
  tr.f1 = o1;
  tr.ns = n;
  create_spectrum(n, o1, d1, spectrum);

  /* Squaring the spectrum */
  j = n/2;
  for(int i = 0; i < nw; i++, j++) {
    rspectrum[i] = spectrum[j]*spectrum[j];
    tr.data[i] = rspectrum[i];
  }

 
  fprintf(stderr, "Created input: \n"); 
  for(int i = 0; i < nw; i++) {
    fprintf(stderr, "i=%d input=%f\n",i,rspectrum[i]);
  }
  fprintf(stderr, "\n");

  // Take the log and create a complex type
  fprintf(stderr, "Log of spectrum: \n"); 
  for(int i = 0; i < nw; i++) {
    fft[i] = cmplx(log(rspectrum[i]+eps)/nfft,0.);
  }
  for(int i = 0; i < nw; i++) {
    fprintf(stderr, "i=%d real=%f imag=%f\n",i,fft[i].r,fft[i].i);
  }
  fprintf(stderr, "\n");

  // Find the inverse FFT
  kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data);

  tr.data[0]      *= 0.5;
  tr.data[nfft/2] *= 0.5;
  for(int i=1+nfft/2; i < nfft; i++) {
    tr.data[i] = 0;
  }

  kiss_fftr(forw, tr.data, (kiss_fft_cpx *) fft);

  for(int i=0; i < nw; i++) {
    fft[i] = crmul(cwp_cexp(fft[i]),1./nfft);
  }

  kiss_fftri(invs,(const kiss_fft_cpx *) fft, tr.data);

  for(int i = 0; i < nfft; i++) {
    fprintf(stderr, "i=%d output=%f\n", i, tr.data[i]);
  }

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
