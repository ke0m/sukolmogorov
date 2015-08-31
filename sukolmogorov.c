/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUKOLMOGOROV: $Revision: 1.0 $ ; $Date: 2015/07/17 19:45:15 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
  "									     ",
  " SUKOLMOGOROFF - The kolmogorff method for spectral factorization",
  "									     ",
  "  sukolmogoroff <stdin >stdout [optional parameters]				     ",
  "									     ",
  " Optional Parameters:							     ",
  " log=0              =1 to take in account that input is already the log of spectrum  ",
  " specsq=0           =1 to take in account that input is already the spectrum sqaured ",
  " verbose=0          =1 for advisory messages",
  " lag=0                 lag for cosine window applied to the trace",
  " Notes: Forthcoming...								     ",
  "									     ",
  NULL};

/* Credits: Based on Sergey Fomel's sfkolmog
 *	SEP: 
 *	CWP: 
 *	TCCS: 
 * Technical Reference: TODO: a geophysics style biblio for GIEE.
 *  
 *
 * Trace header fields accessed: ns, trid, d1
 */

/*
 * TODO:
 * 1. Apply a shift to line up the traces
 *  Look for the first non-zero sample
 *  This sample will be the tmin for sushift
 *  Then shift the trace to this point.
 */

/**************** end self doc *******************************************/

#define LOOKFAC 2 /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  /* Largest allowed nfft           */
#define AMP_SPEC 118 /* trid for the amplitude spectrum of a complex trace */

segy tr;

  int
main(int argc, char **argv)
{
  int nt;           /* number of time samples */
  size_t ntsize;    /* Size of nf in bytes */
  int nf;			      /* number of frequency samples. Right now is 2*nt*/
  int nfft;         /* For computing the FFT */
  size_t nzeros;    /* Size of nf in bytes */
  int i;			      /* time sample index */
  int db;           /* Check if the trace is in dB (for mie scattering) */
  int specsq;       /* Check if the input is the spectrum or spectrum squared */
  int verbose;      /* flag to get advisory messages */
  int lag;          /* lag for asymmetric part */
  float asym;       /*(?) For applying the cosine window to the trace */
  float newd1;      /* New sampling interval for amplitude spectra input */
  cwp_Bool seismic; /* Check if the input is seismic data */

  const double eps = 1.e-32; /* Small value for log in case of log(0) */
  complex *ct      = NULL;   /* Complex array for fft */
  float   *rt      = NULL;   /* Real array for fft */

  /* hook up getpar */
  initargs(argc, argv);
  requestdoc(1);

  /* Getting and setting parameters */
  if (!getparint("lag"    ,  &lag    ))    lag     = 0;
  if (!getparint("db"     ,  &db     ))    db      = 0;
  if (!getparint("specsq" ,  &specsq ))    specsq  = 0;
  if (!getparint("verbose",  &verbose))    verbose = 0;

  /* Read in the trace */
  if (!gettr(&tr)) err("can't get first trace");
  seismic = ISSEISMIC(tr.trid);

  if (seismic) {
    if (verbose)  warn("input is seismic data, trid=%d",tr.trid);
    nt   = tr.ns;
    nfft = 2*npfaro((nt+1)/2,LOOKFAC*nt);
    if (nfft >= SU_NFLTS || nfft >= PFA_MAX) 
      err("Padded nt=%d -- too big", nfft);
    nf   = nfft/2+1;

    ntsize = nt*FSIZE;
    nzeros = (nfft-nt)*FSIZE;

    /* Allocating memory */
    rt = ealloc1float(nfft);
    ct = alloc1complex(nf);
    memcpy( (void *) rt, (const void *) tr.data, ntsize);
    memset((void *) (rt + nt), 0, nzeros);

    /* Forward FFT */
    pfarc(1, nfft, rt, ct);

    /* Computing square of amplitude spectrum */
    for(i = 0; i < nf; i++) {
      rt[i] = (ct[i].r * ct[i].r) + (ct[i].i * ct[i].i);
    }
  }
  else if(tr.trid == AMP_SPEC) {
    //TODO: How do I know what is the correct nfft?
    if (verbose)  warn("input is amplitude spectrum data, trid=%d",tr.trid);
    nt = tr.ns;
    nf = nt; //Should I scale by 2?
    ntsize = nt*FSIZE;
    nfft = npfaro(nf,LOOKFAC*nf);
    if (nfft >= SU_NFLTS || nfft >= PFA_MAX) 
      err("Padded nt=%d -- too big", nfft);
    nzeros = (nfft-nt)*FSIZE;

    /* Computing sampling interval */
    if(tr.d1) {
      newd1 = 1/(nfft*tr.d1); //TODO:I am not sure if this is correct...
    }
    else {
      if(tr.dt) newd1 = (float) (((double) tr.dt)/1000000.0);
      else {
        warn("You did not specify a d1 or dt. Setting d1 to 1...");
        newd1 = 1.f;
      }
    }

    /* Allocating memory */
    rt = ealloc1float(nfft);
    ct = alloc1complex(nf);
    memcpy( (void *) rt, (const void *) tr.data, ntsize);
    memset((void *) (rt + nt), 0, nzeros);

  }
  else {
    err("I do not know how to work with the data you provided. \
        I can only work on time series or amplitude spectra");
  }
  
  /* Squaring the spectrum */
  if(!specsq) {
    for(i = 0; i < nf; i++) {
      rt[i] *= rt[i];
    }
  }

  /* Take the log and create a complex type */
  if(!db) {
    for(i = 0; i < nf; i++) {
      ct[i].r = log(rt[i]+eps)/nfft;
      ct[i].i = 0.f;
    }
  }
  else {
    for(i = 0; i < nf; i++) {
      // Undoing the dB
      rt[i] /= 20.f;
      rt[i] = pow(rt[i],10.f);
      // Squaring the spectrum
      rt[i] *= rt[i];
      // Finding the log of the spectrum
      ct[i].r = log(rt[i]+eps)/nfft;
      ct[i].i = 0.f;
    }
  }

  /* Find the inverse FFT */
  pfacr(-1, nfft, ct, rt);
  
  rt[0]      *= 0.5;
  rt[nfft/2] *= 0.5;
  for(i=1+nfft/2; i < nfft; i++) {
    rt[i] = 0.;
  }

  //TODO: I am really not sure how to use this...
  for (i=1; i < lag; i++) {
    asym = cosf(0.5*PI*i/(lag-1.0));  /* tapering weight */
    asym *= 0.5*rt[i]*asym; 
    rt[i]      -= asym;
    rt[nfft-i] += asym;
  }

  pfarc(1,nfft,rt,ct);

  for(int i=0; i < nf; i++) {
    ct[i] = crmul(cwp_cexp(ct[i]),1./nfft);
  }

  pfacr(-1,nfft,ct,rt);

  /* Copy data back to tr */
  if(seismic) {
    for(i = 0; i < nt; i++) {
      tr.data[i] = rt[i];
    }
  }
  else {
    if(tr.trid == AMP_SPEC) {
      tr.trid = TREAL;
      tr.ns = nfft;
      tr.d1 = newd1;
      if(!tr.f1) {
        tr.f1 = 0.f;
      }
      /* Copy back to tr */
      for(i = 0; i < nfft; i++) {
        tr.data[i] = rt[i];
      }
    }
    else {
      err("I do not know how to work with the data you provided. \
        I can only work on time series or amplitude spectra");
    }
  }
  
  puttr(&tr);

  /* Deallocate memory */
  free1float(rt);
  free1complex(ct);

  return(CWP_Exit());
}

