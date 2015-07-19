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
  " verbose=0          =1 for advisory messages",
  " Notes: Forthcoming...								     ",
  "									     ",
  NULL};

/* Credits:
 *	SEP: 
 *	CWP: 
 *	TCCS: 
 * Technical Reference:
 *
 * Trace header fields accessed: ns,
 */

/*
 * TODO:
 * 1. Put in the cosine taper window
 * 2. Take in a time domain (or other type of signal)
 *    a. With this, make sure that the sample rates and lengths are the same
 * 3. Test code on 1,-0.9 and -0.9,1 wavelets
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
  int logar;        /* Check if the trace has had a log applied (for mie scattering) */
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
  if (!getparint("log"    ,  &logar  ))    logar   = 0;
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
    fprintf(stderr, "nt=%d nfft=%d nf=%d \n", nt,nfft,nf);
  }
  else if(tr.trid == AMP_SPEC) {
    //TODO: How do I know what is the correct nfft?
    // Attempt to mimic what suifft does.
    if (verbose)  warn("input is amplitude spectrum data, trid=%d",tr.trid);
    nt   = tr.ns;
    nf = 2*nt; //Should I scale by 2?
    ntsize = nt*FSIZE;
    nfft = npfaro(nf,LOOKFAC*nf);
    if (nfft >= SU_NFLTS || nfft >= PFA_MAX) 
      err("Padded nt=%d -- too big", nfft);
    nzeros = (nfft-nt)*FSIZE;

    /* Computing sampling interval */
    if(tr.d1) {
      newd1 = 1/(nfft*tr.d1);
    }
    else {
      if(tr.dt) newd1 = (float) (((double) tr.dt)/1000000.0);
      else {
        warn("You did not specify a d1 or dt. Setting to 1...");
        newd1 = 1.f;
      }
    }

    /* Allocating memory */
    rt = ealloc1float(nfft);
    ct = alloc1complex(nf);
    memcpy( (void *) rt, (const void *) tr.data, ntsize);
    memset((void *) (rt + nt), 0, nzeros);

    /* Squaring the spectrum */
    if(!specsq) {
      for(i = 0; i < nf; i++) {
        rt[i] *= rt[i];
      }
    }
    fprintf(stderr, "nt=%d nfft=%d nf=%d \n", nt,nfft,nf);
  }
  else {
    err("I do not know how to work with the data you provided. \
        I can only work on time series or amplitude spectra");
  }

  fprintf(stderr, "input spec: \n");
  for(i = 0; i < nfft; i++) {
    fprintf(stderr, "i=%d input=%f\n",i,rt[i]);
  }
  fprintf(stderr, "\n");

  /* Take the log and create a complex type */
  if(!logar) {
    for(i = 0; i < nf; i++) { //This needs to be the length of the data
      ct[i].r = log(rt[i]+eps)/nfft; // Therefore, the nfft needs to be twice that length
      ct[i].i = 0.f;
    }
  }
  else {
    for(i = 0; i < nf; i++) {
      ct[i].r = rt[i];
      ct[i].i = 0.f;
    }
  }

  //fprintf(stderr, "log of spec(nfft=%d): \n",nfft);
  //for(i = 0; i < nf; i++) {
  //  fprintf(stderr, "i=%d real=%f imag=%f\n",i,ct[i].r,ct[i].i);
  //}
  //fprintf(stderr, "\n");

  /* Find the inverse FFT */
  //fprintf(stderr, "the ifft: \n");
  pfacr(-1, nfft, ct, rt);
  /*
     for(i = 0; i < nfft; i++) {
     fprintf(stderr, "i=%d ifft=%f\n",i,tr.data[i]);
     }
     fprintf(stderr, "\n");
     */

  rt[0]      *= 0.5;
  rt[nfft/2] *= 0.5;
  for(i=1+nfft/2; i < nfft; i++) {
    rt[i] = 0.;
  }

  // Put in cosine taper window
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

  fprintf(stderr, "output: \n");
  pfacr(-1,nfft,ct,rt);

  float o1 = -M_PI/2;
  float dt = 0.004;
  //for(i = 0; i < nfft; i++) {
  //  fprintf(stderr, "i=%d t=%f output=%f\n",i,o1+i*dt,rt[i]);
  //}

  /* Copy data back to tr */
  if(seismic) {
    for(i = 0; i < nt; i++) { // I am not sure but this might cause problems for the time domain case???
      tr.data[i] = rt[i];
    }
  }
  else {
    //TODO: What other fields should I set here:
    //      d1, dt, ns?
    if(tr.trid == AMP_SPEC) {
      tr.trid = TREAL;
      tr.ns = nfft;
      //tr.d1 = newd1;
      tr.dt = 0.004; // for testing...
      if(!tr.f1) {
        tr.f1 = 0.f;
      }
      /* Copy back to tr */
      for(i = 0; i < nfft; i++) {
        tr.data[i] = rt[i];
        fprintf(stderr, "i=%d t=%f output=%f\n",i,o1+i*dt,rt[i]);
      }
    }
  }
  
  puttr(&tr);

  /* Deallocate memory */
  free1float(rt);
  free1complex(ct);

  return(CWP_Exit());
}

