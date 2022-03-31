#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <complex.h> 
#include <fstream>
#include <stdlib.h>
#include<time.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "ldc_code.h"


int main( int argc, const char* argv[] ){
	srand(time(0));
	int n_waves; n_waves=(int)atoi(argv[1]);
	printf("\nn_waves= = %d", n_waves);
	int Tobs, dt, *timeInfo; timeInfo = (int *)malloc(3*sizeof(int));
	double *parameters, *XFFT, *YFFT, *ZFFT;
	parameters = (double *)malloc(NUM_PARAMETERS_GALBIN*sizeof(double));
	readInData("data/sim.txt", "data/parameters.txt", timeInfo, parameters, n_waves); Tobs=timeInfo[1]; dt=timeInfo[2];

	int N =1024;
	XFFT = (double *)malloc(2*N*sizeof(double)); YFFT = (double *)malloc(2*N*sizeof(double)); ZFFT = (double *)malloc(2*N*sizeof(double));
	
	struct Waveform *wfm = (struct Waveform*) malloc(sizeof(struct Waveform));
	Fast_GB(parameters, wfm, N, Tobs, dt, XFFT, YFFT, ZFFT, NUM_PARAMETERS_GALBIN);
	for (int ii=54; ii<64; ii++){
		printf("\n XFFT[%d]=%f", ii, XFFT[ii]*pow(10,18));
	}
	
    return 0;
}

int buffersize(double f0, double T){
	int mult=1;
	if(T <= 8.0*YEAR){mult = 8;}
	if(T <= 4.0*YEAR){mult = 4;}
	if(T <= 2.0*YEAR){mult = 2;}
	if(T <= 1.0*YEAR){mult = 1;}
	int N=0;
	if(f0>0.001){N=64*mult;}
	else if(f0>0.01){N=256*mult;}
	else if(f0>0.03){N=512*mult;}
	else if(f0>0.1){N=1024*mult;}
	else{N=32*mult;}
	return N;
}

void readInData(const char file_sim[19], const char file_GB_parameter[19], int *timeInfo, double *parameters, int n_waves){
    FILE *f;
    f = fopen(file_sim,"r");
	int t0, T, dt;
    if (f != NULL){ 
		fscanf(f, "%d %d %d", &t0,&T,&dt);
		timeInfo[0]=t0;timeInfo[1]=T;timeInfo[2]=dt;
		fclose(f);
		}

	FILE *g;
    g = fopen(file_GB_parameter,"r");
	int i = 0;
    long double datatmp[NUM_PARAMETERS_GALBIN];
	
    while (i <= n_waves){
        fscanf(g, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf", &datatmp[0],&datatmp[1],&datatmp[2],&datatmp[3],&datatmp[4],&datatmp[5],&datatmp[6],&datatmp[7]);

		printf("\n parameters[0]=%Lf", datatmp[0]);
		printf("\n parameters[1]=%Lf", datatmp[1]);
		printf("\n parameters[2]=%Lf", datatmp[2]);
		printf("\n parameters[3]=%Lf", datatmp[3]);
		printf("\n parameters[4]=%Lf", datatmp[4]);
		printf("\n parameters[5]=%Lf", datatmp[5]);
		printf("\n parameters[6]=%Lf", datatmp[6]);
		printf("\n parameters[7]=%Lf", datatmp[7]);

		// datatmp[0]=M_PI/2;
		// datatmp[1]=M_PI/2;
		// datatmp[2]=M_PI/2;
		// datatmp[3]=1;
		// datatmp[4]=0;
		// datatmp[5]=1e-3+1e-6;
		// datatmp[6]=0;
		// datatmp[7]=0;
		if (i == n_waves) {
			// datatmp[5]=(double)((round(datatmp[5+n_waves*7]*pow(10,9)))*pow(10,-9));
			parameters[0]=datatmp[5]*T; // frequency
			parameters[1]=cos(PI/2-datatmp[1]); // cos(theta) i.e. cos(pi/2-beta)
			parameters[2]=datatmp[0]; // phi i.e. lambda
			parameters[3]=log(datatmp[3]); // amplitude
			parameters[4]=cos(datatmp[4]); // cos(iota)
			parameters[5]=datatmp[2]; // psi
			parameters[6]=-1*datatmp[6]; // phi_0
		}
		i++;
    }
    fclose(g);
	printf("\n parameters[0]=%f", parameters[0]);
	printf("\n parameters[1]=%f", parameters[1]);
	printf("\n parameters[2]=%f", parameters[2]);
	printf("\n parameters[3]=%f", parameters[3]);
	printf("\n parameters[4]=%f", parameters[4]);
	printf("\n parameters[5]=%f", parameters[5]);
	printf("\n parameters[6]=%f", parameters[6]);
}

void Fast_GB(double *params, struct Waveform* wfm, long N, double Tobs, double dt, double *XSL, double *YSL, double *ZSL, int NP)
{		
	long n; // iterator
	double t; // time
	double* orbit_params = (double*) malloc(3*sizeof(double));
  	orbit_params[0] = LARM; orbit_params[1] = LAMBDA; orbit_params[2] = KAPPA;

	struct LISA *lisa = (struct LISA*) malloc(sizeof(struct LISA));
	init_lisa(lisa, orbit_params[0], orbit_params[1], orbit_params[2]);
	double fstar = (CLIGHT/lisa->spacecraft_separation)/(2*PI); //0.01908538063694777;

	wfm->N  = N; // set number of samples
	wfm->T  = Tobs; // set observation period
	wfm->NP = NP; // inform model of number of parameters being used
	alloc_waveform(wfm); // allocate memory to hold pieces of waveform
	copy_params(wfm, params); // copy parameters to waveform structure

	get_basis_tensors(wfm); //  Tensor construction for building slowly evolving LISA response
	for(n=0; n<N; n++)
	{
	  	t = wfm->T*(double)(n)/(double)N; // First time sample must be at t=0 for phasing
	  	calc_xi_f(wfm, lisa, fstar, t);  // calc frequency and time variables
	  	calc_sep_vecs(wfm, lisa->spacecraft_separation); // calculate the S/C separation vectors
	  	calc_d_matrices(wfm); // calculate pieces of waveform
	  	calc_kdotr(wfm); // calculate dot product
	  	get_transfer(wfm, t); // Calculating Transfer function
	  	fill_time_series(wfm, n); // Fill  time series data arrays with slowly evolving signal.
	}
	fft_data(wfm); // Numerical Fourier transform of slowly evolving signal
	unpack_data(wfm); // Unpack arrays from FFT and normalize

	XYZ(wfm->d, wfm->params[0]/wfm->T, wfm->q, N, dt, Tobs, lisa->spacecraft_separation, fstar, XSL, YSL, ZSL);

	free_waveform(wfm); // Deallocate memory
	free(wfm);

	return;
}

void calc_xi_f(struct Waveform *wfm, struct LISA* lisa, double fstar, double t)
{
	long i;

	double f0;

	f0       = wfm->params[0]/wfm->T;

	// Calculate position of each spacecraft at time t
	for(i=0; i<3; i++){
	  wfm->x[i] = AnalyticOrbits_get_position_x(lisa, i+1, t);
	  wfm->y[i] = AnalyticOrbits_get_position_y(lisa, i+1, t);
	  wfm->z[i] = AnalyticOrbits_get_position_z(lisa, i+1, t);
	}

	for(i=0; i<3; i++)
	{
		wfm->kdotx[i] = (wfm->x[i]*wfm->k[0] + wfm->y[i]*wfm->k[1] + wfm->z[i]*wfm->k[2])/CLIGHT;
		//Wave arrival time at spacecraft i
		wfm->xi[i]    = t - wfm->kdotx[i];
		//First order approximation to frequency at spacecraft i
		wfm->f[i]     = f0;
		//Ratio of true frequency to transfer frequency
		wfm->fonfs[i] = wfm->f[i]/fstar;

		double tmp=wfm->xi[i];

		if(i==0){
			FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/xi1.bin","ab");
			fwrite(&tmp, sizeof(double), 1, ptrPS);
			fclose(ptrPS);
		}
		if(i==1){
			FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/xi2.bin","ab");
			fwrite(&tmp, sizeof(double), 1, ptrPS);
			fclose(ptrPS);
		}
		if(i==2){
			FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/xi3.bin","ab");
			fwrite(&tmp, sizeof(double), 1, ptrPS);
			fclose(ptrPS);
		}
	}

	return;
}

void copy_params(struct Waveform *wfm, double *params)
{
//	long i;
    int i;

	int NP = wfm->NP;

	wfm->params = (double*) malloc(NP*sizeof(double));

	for (i=0; i<NP; i++){
	    wfm->params[i] = params[i];
	    }

	wfm->q  = (long)(params[0]); //Calculate carrier frequency bin

	return;
}

void fill_time_series(struct Waveform *wfm, int n)
{
	wfm->data12[2*n]   = wfm->TR[0][1];
	wfm->data21[2*n]   = wfm->TR[1][0];
	wfm->data31[2*n]   = wfm->TR[2][0];
	wfm->data12[2*n+1] = wfm->TI[0][1];
	wfm->data21[2*n+1] = wfm->TI[1][0];
	wfm->data31[2*n+1] = wfm->TI[2][0];
	wfm->data13[2*n]   = wfm->TR[0][2];
	wfm->data23[2*n]   = wfm->TR[1][2];
	wfm->data32[2*n]   = wfm->TR[2][1];
	wfm->data13[2*n+1] = wfm->TI[0][2];
	wfm->data23[2*n+1] = wfm->TI[1][2];
	wfm->data32[2*n+1] = wfm->TI[2][1];


   FILE *ptr12=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y12.bin","ab");
   double tmpR=wfm->data12[2*n];
   double tmpI=wfm->data12[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr12);
   fwrite(&tmpI, sizeof(double), 1, ptr12);
   fclose(ptr12);

   FILE *ptr13=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y13.bin","ab");
   tmpR=wfm->data13[2*n];
   tmpI=wfm->data13[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr13);
   fwrite(&tmpI, sizeof(double), 1, ptr13);
   fclose(ptr13);

   FILE *ptr23=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y23.bin","ab");
   tmpR=wfm->data23[2*n];
   tmpI=wfm->data23[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr23);
   fwrite(&tmpI, sizeof(double), 1, ptr23);
   fclose(ptr23);

   FILE *ptr21=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y21.bin","ab");
   tmpR=wfm->data21[2*n];
   tmpI=wfm->data21[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr21);
   fwrite(&tmpI, sizeof(double), 1, ptr21);
   fclose(ptr21);

   FILE *ptr31=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y31.bin","ab");
   tmpR=wfm->data31[2*n];
   tmpI=wfm->data31[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr31);
   fwrite(&tmpI, sizeof(double), 1, ptr31);
   fclose(ptr31);

   FILE *ptr32=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/y32.bin","ab");
   tmpR=wfm->data32[2*n];
   tmpI=wfm->data32[2*n+1];
   fwrite(&tmpR, sizeof(double), 1, ptr32);
   fwrite(&tmpI, sizeof(double), 1, ptr32);
   fclose(ptr32);

	return;
}

void unpack_data(struct Waveform *wfm)
{
	long i;
	long N = wfm->N;

	for(i=0; i<N; i++)
	{	// populate from most negative (Nyquist) to most positive (Nyquist-1)
		wfm->a12[i]   = wfm->data12[N+i]/(double)N;
		wfm->a21[i]   = wfm->data21[N+i]/(double)N;
		wfm->a31[i]   = wfm->data31[N+i]/(double)N;
		wfm->a12[i+N] = wfm->data12[i]/(double)N;
		wfm->a21[i+N] = wfm->data21[i]/(double)N;
		wfm->a31[i+N] = wfm->data31[i]/(double)N;
		wfm->a13[i]   = wfm->data13[N+i]/(double)N;
		wfm->a23[i]   = wfm->data23[N+i]/(double)N;
		wfm->a32[i]   = wfm->data32[N+i]/(double)N;
		wfm->a13[i+N] = wfm->data13[i]/(double)N;
		wfm->a23[i+N] = wfm->data23[i]/(double)N;
		wfm->a32[i+N] = wfm->data32[i]/(double)N;
	}

	//   Renormalize so that the resulting time series is real
	for(i=0; i<2*N; i++)
	{
		wfm->d[0][1][i] = 0.5*wfm->a12[i];
		wfm->d[1][0][i] = 0.5*wfm->a21[i];
		wfm->d[2][0][i] = 0.5*wfm->a31[i];
		wfm->d[0][2][i] = 0.5*wfm->a13[i];
		wfm->d[1][2][i] = 0.5*wfm->a23[i];
		wfm->d[2][1][i] = 0.5*wfm->a32[i];
	}

	return;
}

void fft_data(struct Waveform *wfm)
{
	long N = wfm->N;
	gsl_fft_complex_radix2_forward(wfm->data12, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data21, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data31, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data13, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data23, 1, N);
	gsl_fft_complex_radix2_forward(wfm->data32, 1, N);


	FILE *ptr12=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT12.bin","ab");
	double* tmpR=wfm->data12;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr12);
	fclose(ptr12);

	FILE *ptr13=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT13.bin","ab");
	tmpR=wfm->data13;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr13);
	fclose(ptr13);

	FILE *ptr23=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT23.bin","ab");
	tmpR=wfm->data23;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr23);
	fclose(ptr23);

	FILE *ptr21=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT21.bin","ab");
	tmpR=wfm->data21;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr21);
	fclose(ptr21);

	FILE *ptr31=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT31.bin","ab");
	tmpR=wfm->data31;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr31);
	fclose(ptr31);

	FILE *ptr32=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/Doppler/yFFT32.bin","ab");
	tmpR=wfm->data32;
	fwrite(tmpR, 2*N*sizeof(double), 1, ptr32);
	fclose(ptr32);

	return;
}

void alloc_waveform(struct Waveform *wfm)
{
	long i, j, n;
	long N;

	N = wfm->N;

	wfm->k = (double*) malloc(3*sizeof(double));

	wfm->kdotx = (double*) malloc(3*sizeof(double));
	wfm->kdotr = (double**) malloc(3*sizeof(double *));
	for (i=0; i<3; i++) wfm->kdotr[i] = (double*) malloc(3*sizeof(double));

	for (i=0; i<3; i++)
	{
	  for (j=0; j<3; j++) wfm->kdotr[i][j] = 0.;
	  wfm->kdotx[i] = 0.;
	}

	wfm->xi    = (double*) malloc(3*sizeof(double));
	wfm->f     = (double*) malloc(3*sizeof(double));
	wfm->fonfs = (double*) malloc(3*sizeof(double));
	for (i=0; i<3; i++)
	{
		wfm->xi[i]    = 0.;
		wfm->f[i]     = 0.;
		wfm->fonfs[i] = 0.;
	}

	// Polarization basis tensors
	wfm->eplus = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->eplus[i] = (double*) malloc(3*sizeof(double));
	wfm->ecross = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->ecross[i] = (double*) malloc(3*sizeof(double));

	wfm->dplus = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->dplus[i] = (double*) malloc(3*sizeof(double));
	wfm->dcross = (double**) malloc(3*sizeof(double*));
	for (i=0; i<3; i++) wfm->dcross[i] = (double*) malloc(3*sizeof(double));

	wfm->r12 = (double*) malloc(3*sizeof(double));
	wfm->r21 = (double*) malloc(3*sizeof(double));
	wfm->r31 = (double*) malloc(3*sizeof(double));
	wfm->r13 = (double*) malloc(3*sizeof(double));
	wfm->r23 = (double*) malloc(3*sizeof(double));
	wfm->r32 = (double*) malloc(3*sizeof(double));


	wfm->data12 = (double*) malloc(2*N*sizeof(double));
	wfm->data21 = (double*) malloc(2*N*sizeof(double));
	wfm->data31 = (double*) malloc(2*N*sizeof(double));
	wfm->data13 = (double*) malloc(2*N*sizeof(double));
	wfm->data23 = (double*) malloc(2*N*sizeof(double));
	wfm->data32 = (double*) malloc(2*N*sizeof(double));
	for (i=0; i<2*N; i++)
	{
		wfm->data12[i] = 0.;
		wfm->data21[i] = 0.;
		wfm->data31[i] = 0.;
		wfm->data13[i] = 0.;
		wfm->data23[i] = 0.;
		wfm->data32[i] = 0.;
	}

	wfm->a12 = (double*) malloc(2*N*sizeof(double));
	wfm->a21 = (double*) malloc(2*N*sizeof(double));
	wfm->a31 = (double*) malloc(2*N*sizeof(double));
	wfm->a13 = (double*) malloc(2*N*sizeof(double));
	wfm->a23 = (double*) malloc(2*N*sizeof(double));
	wfm->a32 = (double*) malloc(2*N*sizeof(double));
	for (i=0; i<2*N; i++)
	{
		wfm->a12[i] = 0.;
		wfm->a21[i] = 0.;
		wfm->a31[i] = 0.;
		wfm->a13[i] = 0.;
		wfm->a23[i] = 0.;
		wfm->a32[i] = 0.;
	}

	wfm->TR = (double**) malloc(3*sizeof(double*));
	for (i=0;i<3; i++) wfm->TR[i] = (double*) malloc(3*sizeof(double));
	wfm->TI = (double**) malloc(3*sizeof(double*));
	for (i=0;i<3; i++) wfm->TI[i] = (double*) malloc(3*sizeof(double));

	wfm->x = (double*) malloc(3*sizeof(double));
	wfm->y = (double*) malloc(3*sizeof(double));
	wfm->z = (double*) malloc(3*sizeof(double));

	for (i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			wfm->eplus[i][j]  = 0.;
			wfm->ecross[i][j] = 0.;
			wfm->dplus[i][j]  = 0.;
			wfm->dcross[i][j] = 0.;
			wfm->TR[i][j]     = 0.;
			wfm->TI[i][j]     = 0.;
		}
		wfm->x[i]   = 0.;
		wfm->y[i]   = 0.;
		wfm->z[i]   = 0.;
		wfm->r12[i] = 0.;
		wfm->r21[i] = 0.;
		wfm->r31[i] = 0.;
		wfm->r13[i] = 0.;
		wfm->r23[i] = 0.;
		wfm->r32[i] = 0.;
	}

	wfm->d = (double***) malloc(3*sizeof(double**));
	for (i=0; i<3; i++) wfm->d[i] = (double**) malloc(3*sizeof(double*));

	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			wfm->d[i][j] =(double*) malloc(2*N*sizeof(double));
		}
	}

	for (i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			for (n=0; n<2*N; n++)
			{
				wfm->d[i][j][n] = 0.;
			}
		}
	}

	return;
}

void free_waveform(struct Waveform *wfm)
{
	long i, j;

	free(wfm->kdotx);
	free(wfm->k);
	for (i=0; i<3; i++) free(wfm->kdotr[i]);
	free(wfm->kdotr);

	free(wfm->xi);
	free(wfm->f);
	free(wfm->fonfs);

	for (i=0; i<3; i++){
	  free(wfm->eplus[i]);
	  free(wfm->ecross[i]);
	  free(wfm->dplus[i]);
	  free(wfm->dcross[i]);
	}
	free(wfm->eplus);
	free(wfm->ecross);
	free(wfm->dplus);
	free(wfm->dcross);

	free(wfm->r12);
	free(wfm->r21);
	free(wfm->r31);
	free(wfm->r13);
	free(wfm->r23);
	free(wfm->r32);

	free(wfm->data12);
	free(wfm->data21);
	free(wfm->data31);
	free(wfm->data13);
	free(wfm->data23);
	free(wfm->data32);

	free(wfm->a12);
	free(wfm->a21);
	free(wfm->a31);
	free(wfm->a13);
	free(wfm->a23);
	free(wfm->a32);

	free(wfm->x);
	free(wfm->y);
	free(wfm->z);

	for (i=0; i<3; i++){
	  free(wfm->TR[i]);
	  free(wfm->TI[i]);
	}
	free(wfm->TR);
	free(wfm->TI);

	for (i=0; i<3; i++)
	  for(j=0; j<3; j++)
	    free(wfm->d[i][j]);

	for (i=0; i<3; i++)
	  free(wfm->d[i]);
	free(wfm->d);

	free(wfm->params);

	return;
}

void calc_d_matrices(struct Waveform *wfm)
{
	long i, j;

	//Zero arrays to be summed
	wfm->dplus [0][1] = wfm->dplus [0][2] = wfm->dplus [1][0] = 0.;
	wfm->dplus [1][2] = wfm->dplus [2][0] = wfm->dplus [2][1] = 0.;
	wfm->dcross[0][1] = wfm->dcross[0][2] = wfm->dcross[1][0] = 0.;
	wfm->dcross[1][2] = wfm->dcross[2][0] = wfm->dcross[2][1] = 0.;

	//Convenient quantities d+ & dx
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			wfm->dplus [0][1] += wfm->r12[i]*wfm->r12[j]*wfm->eplus[i][j];
			wfm->dcross[0][1] += wfm->r12[i]*wfm->r12[j]*wfm->ecross[i][j];
			wfm->dplus [1][2] += wfm->r23[i]*wfm->r23[j]*wfm->eplus[i][j];
			wfm->dcross[1][2] += wfm->r23[i]*wfm->r23[j]*wfm->ecross[i][j];
			wfm->dplus [0][2] += wfm->r13[i]*wfm->r13[j]*wfm->eplus[i][j];
			wfm->dcross[0][2] += wfm->r13[i]*wfm->r13[j]*wfm->ecross[i][j];
		}
	}
	//Makng use of symmetry
	wfm->dplus[1][0] = wfm->dplus[0][1];  wfm->dcross[1][0] = wfm->dcross[0][1];
	wfm->dplus[2][1] = wfm->dplus[1][2];  wfm->dcross[2][1] = wfm->dcross[1][2];
	wfm->dplus[2][0] = wfm->dplus[0][2];  wfm->dcross[2][0] = wfm->dcross[0][2];

	return;
}

void calc_sep_vecs(struct Waveform *wfm, double Larm)
{
	long i;

	//Unit separation vector from spacecrafts i to j
	wfm->r12[0] = (wfm->x[1] - wfm->x[0])/Larm;
	wfm->r13[0] = (wfm->x[2] - wfm->x[0])/Larm;
	wfm->r23[0] = (wfm->x[2] - wfm->x[1])/Larm;
	wfm->r12[1] = (wfm->y[1] - wfm->y[0])/Larm;
	wfm->r13[1] = (wfm->y[2] - wfm->y[0])/Larm;
	wfm->r23[1] = (wfm->y[2] - wfm->y[1])/Larm;
	wfm->r12[2] = (wfm->z[1] - wfm->z[0])/Larm;
	wfm->r13[2] = (wfm->z[2] - wfm->z[0])/Larm;
	wfm->r23[2] = (wfm->z[2] - wfm->z[1])/Larm;

	//Make use of symmetry
	for(i=0; i<3; i++)
	{
		wfm->r21[i] = -wfm->r12[i];
		wfm->r31[i] = -wfm->r13[i];
		wfm->r32[i] = -wfm->r23[i];
	}

	double tmpx=wfm->r12[0];
	double tmpy=wfm->r12[1];
	double tmpz=wfm->r12[2];
	FILE *ptrPS1=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/n12.bin","ab");
	fwrite(&tmpx, sizeof(double), 1, ptrPS1);
	fwrite(&tmpy, sizeof(double), 1, ptrPS1);
	fwrite(&tmpz, sizeof(double), 1, ptrPS1);
	fclose(ptrPS1);
	tmpx=wfm->r23[0];
	tmpy=wfm->r23[1];
	tmpz=wfm->r23[2];
	FILE *ptrPS2=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/n23.bin","ab");
	fwrite(&tmpx, sizeof(double), 1, ptrPS2);
	fwrite(&tmpy, sizeof(double), 1, ptrPS2);
	fwrite(&tmpz, sizeof(double), 1, ptrPS2);
	fclose(ptrPS2);
	tmpx=wfm->r31[0];
	tmpy=wfm->r31[1];
	tmpz=wfm->r31[2];
	FILE *ptrPS3=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/n31.bin","ab");
	fwrite(&tmpx, sizeof(double), 1, ptrPS3);
	fwrite(&tmpy, sizeof(double), 1, ptrPS3);
	fwrite(&tmpz, sizeof(double), 1, ptrPS3);
	fclose(ptrPS3);


	return;
}

void calc_kdotr(struct Waveform *wfm)
{
	long i;

	//Zero arrays to be summed
	wfm->kdotr[0][1] = wfm->kdotr[0][2] = wfm->kdotr[1][0] = 0.;
	wfm->kdotr[1][2] = wfm->kdotr[2][0] = wfm->kdotr[2][1] = 0.;

	for(i=0; i<3; i++)
	{
		wfm->kdotr[0][1] += wfm->k[i]*wfm->r12[i];
		wfm->kdotr[0][2] += wfm->k[i]*wfm->r13[i];
		wfm->kdotr[1][2] += wfm->k[i]*wfm->r23[i];
	}

	//Making use of antisymmetry
	wfm->kdotr[1][0] = -wfm->kdotr[0][1];
	wfm->kdotr[2][0] = -wfm->kdotr[0][2];
	wfm->kdotr[2][1] = -wfm->kdotr[1][2];

	return;
}

void set_const_trans(struct Waveform *wfm)
{
	double amp, cosiota;
	double Aplus, Across;
	double psi;
	double sinps, cosps;

	amp      = exp(wfm->params[3]);
	cosiota  = wfm->params[4];
	psi      = wfm->params[5];

	//Calculate GW polarization amplitudes
	Aplus  = amp*(1. + cosiota*cosiota);
	Across = -2.0*amp*cosiota;

	//Calculate cos and sin of polarization
	cosps = cos(2.*psi);
	sinps = sin(2.*psi);

	//Calculate constant pieces of transfer functions
	wfm->DPr    =  Aplus*cosps;
	wfm->DPi    = -Across*sinps;
	wfm->DCr    = -Aplus*sinps;
	wfm->DCi    = -Across*cosps;

	return;
}

void get_basis_tensors(struct Waveform *wfm)
{
	long i, j;

	double *u, *v;		  // GW basis vectors

	u = (double*) malloc(3*sizeof(double));
	v = (double*) malloc(3*sizeof(double));

	set_const_trans(wfm);  // set the constant pieces of transfer function

	get_basis_vecs(wfm->params, u, v, wfm->k); //Gravitational Wave source basis vectors
	//GW polarization basis tensors
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			wfm->eplus[i][j]  = v[i]*v[j] - u[i]*u[j];
			wfm->ecross[i][j] = u[i]*v[j] + v[i]*u[j];
			double tmpp=wfm->eplus[i][j];
			double tmpc=wfm->ecross[i][j];
			FILE *ptrp=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/eplus.bin","ab");
			fwrite(&tmpp, sizeof(double), 1, ptrp);
			fclose(ptrp);

			FILE *ptrc=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/ecross.bin","ab");
			fwrite(&tmpc, sizeof(double), 1, ptrc);
			fclose(ptrc);
		}
	}

	free(u);
	free(v);

	return;
}

void get_basis_vecs(double *params, double *u, double *v, double *k)
{
	long i;

	double phi;
	double costh, sinth, cosph, sinph;

	for (i=0; i<3; i++)
	{
		u[i] = 0.;
		v[i] = 0.;
		k[i] = 0.;
	}

	phi	  = params[2];
	costh = params[1];

	sinth = sqrt(1.0-costh*costh);

	cosph = cos(phi);
	sinph = sin(phi);

	u[0] =  costh*cosph;  u[1] =  costh*sinph;  u[2] = -sinth;
	v[0] =  sinph;        v[1] = -cosph;        v[2] =  0.;
	k[0] = -sinth*cosph;  k[1] = -sinth*sinph;  k[2] = -costh;

	printf("k[0]=%f, k[1]=%f, k[2]=%f", k[0],k[1],k[2]);
	return;
}

void get_transfer(struct Waveform *wfm, double t)
{
	long i, j;
	long q;

	double tran1r, tran1i;
	double tran2r, tran2i;
	double aevol;			// amplitude evolution factor
	double arg1, arg2, sinc;
	double f0;
	double df, phi0;

	f0       = wfm->params[0]/wfm->T;
	// printf("\n params[0]=%f, f0=%f , q=%ld", wfm->params[0], wfm->params[0]/wfm->T, wfm->q);
	phi0     = wfm->params[6];

	q  = wfm->q;
	df = PI*2*(((double)q)/wfm->T);

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if(i!=j)
			{
				//Argument of transfer function
				arg1 = 0.5*wfm->fonfs[i]*(1. + wfm->kdotr[i][j]);
				//Argument of complex exponentials
				arg2 =  PI*2*f0*wfm->xi[i] + phi0 - df*t;
				//Transfer function
				sinc = 0.25*sin(arg1)/arg1;
				//Evolution of amplitude
				aevol = 1.0;
				///Real and imaginary pieces of time series (no complex exponential)
				tran1r = aevol*(wfm->dplus[i][j]*wfm->DPr + wfm->dcross[i][j]*wfm->DCr);
				tran1i = aevol*(wfm->dplus[i][j]*wfm->DPi + wfm->dcross[i][j]*wfm->DCi);
				//Real and imaginry components of complex exponential
				tran2r = cos(arg1 + arg2);
				tran2i = sin(arg1 + arg2);
				//Real & Imaginary part of the slowly evolving signal
				wfm->TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i);
				wfm->TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);

				double phaseShiftR = sinc*tran2r;
                double phaseShiftI = sinc*tran2i;

				if(i==0 && j==1){
					FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/phaseShift12.bin","ab");
					fwrite(&phaseShiftR, sizeof(double), 1, ptrPS);
					fwrite(&phaseShiftI, sizeof(double), 1, ptrPS);
					fclose(ptrPS);
					FILE *ptr2=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/dApsi12.bin","ab");
					fwrite(&tran1r, sizeof(double), 1, ptr2);
                    fwrite(&tran1i, sizeof(double), 1, ptr2);
                    fclose(ptr2);
				}
				if(i==1 && j==2){
					FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/phaseShift23.bin","ab");
					fwrite(&phaseShiftR, sizeof(double), 1, ptrPS);
					fwrite(&phaseShiftI, sizeof(double), 1, ptrPS);
					fclose(ptrPS);
					FILE *ptr2=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/dApsi23.bin","ab");
					fwrite(&tran1r, sizeof(double), 1, ptr2);
                    fwrite(&tran1i, sizeof(double), 1, ptr2);
                    fclose(ptr2);
				}
				if(i==2 && j==0){
					FILE *ptrPS=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/phaseShift31.bin","ab");
					fwrite(&phaseShiftR, sizeof(double), 1, ptrPS);
					fwrite(&phaseShiftI, sizeof(double), 1, ptrPS);
					fclose(ptrPS);
					FILE *ptr2=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/zwischenergebnisse/dApsi31.bin","ab");
					fwrite(&tran1r, sizeof(double), 1, ptr2);
                    fwrite(&tran1i, sizeof(double), 1, ptr2);
                    fclose(ptr2);
				}
			}
		}
	}

	return;
}

void XYZ(double ***d, double f0, long q, long M, double dt, double Tobs, double Larm, double fstar,double *XSL, double *YSL, double *ZSL)
{
	int i;
	double fonfs;
	double c3, s3, c2, s2, c1, s1;
	double f;
	double *X, *Y, *Z;
	double *XX, *YY, *ZZ;
	double phiSL, cSL, sSL;

	X   = (double*) malloc(2*M*sizeof(double)); XX   = (double*) malloc(2*M*sizeof(double));
	Y   = (double*) malloc(2*M*sizeof(double)); YY   = (double*) malloc(2*M*sizeof(double));
	Z   = (double*) malloc(2*M*sizeof(double)); ZZ   = (double*) malloc(2*M*sizeof(double));

	phiSL = PI/2.0-2.0*PI*f0*(Larm/CLIGHT);
	cSL = cos(phiSL);
	sSL = sin(phiSL);

	for(i=0; i<M; i++)
	{	
		f = ((double)(q + i - M/2))/Tobs;
		fonfs = f/fstar;
		c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(1.*fonfs);
		s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(1.*fonfs);

		X[2*i]   = (d[0][1][2*i]-d[0][2][2*i])*c3 + (d[0][1][2*i+1]-d[0][2][2*i+1])*s3 +
		           (d[1][0][2*i]-d[2][0][2*i])*c2 + (d[1][0][2*i+1]-d[2][0][2*i+1])*s2 +
		           (d[0][2][2*i]-d[0][1][2*i])*c1 + (d[0][2][2*i+1]-d[0][1][2*i+1])*s1 +
		           (d[2][0][2*i]-d[1][0][2*i]);

		X[2*i+1] = (d[0][1][2*i+1]-d[0][2][2*i+1])*c3 - (d[0][1][2*i]-d[0][2][2*i])*s3 +
		           (d[1][0][2*i+1]-d[2][0][2*i+1])*c2 - (d[1][0][2*i]-d[2][0][2*i])*s2 +
		           (d[0][2][2*i+1]-d[0][1][2*i+1])*c1 - (d[0][2][2*i]-d[0][1][2*i])*s1 +
		           (d[2][0][2*i+1]-d[1][0][2*i+1]);

		Y[2*i]   = (d[1][2][2*i]-d[1][0][2*i])*c3 + (d[1][2][2*i+1]-d[1][0][2*i+1])*s3 +
		           (d[2][1][2*i]-d[0][1][2*i])*c2 + (d[2][1][2*i+1]-d[0][1][2*i+1])*s2+
		           (d[1][0][2*i]-d[1][2][2*i])*c1 + (d[1][0][2*i+1]-d[1][2][2*i+1])*s1+
		           (d[0][1][2*i]-d[2][1][2*i]);

		Y[2*i+1] = (d[1][2][2*i+1]-d[1][0][2*i+1])*c3 - (d[1][2][2*i]-d[1][0][2*i])*s3+
		           (d[2][1][2*i+1]-d[0][1][2*i+1])*c2 - (d[2][1][2*i]-d[0][1][2*i])*s2+
		           (d[1][0][2*i+1]-d[1][2][2*i+1])*c1 - (d[1][0][2*i]-d[1][2][2*i])*s1+
		           (d[0][1][2*i+1]-d[2][1][2*i+1]);

		Z[2*i]   = (d[2][0][2*i]-d[2][1][2*i])*c3 + (d[2][0][2*i+1]-d[2][1][2*i+1])*s3+
		           (d[0][2][2*i]-d[1][2][2*i])*c2 + (d[0][2][2*i+1]-d[1][2][2*i+1])*s2+
		           (d[2][1][2*i]-d[2][0][2*i])*c1 + (d[2][1][2*i+1]-d[2][0][2*i+1])*s1+
		           (d[1][2][2*i]-d[0][2][2*i]);

		Z[2*i+1] = (d[2][0][2*i+1]-d[2][1][2*i+1])*c3 - (d[2][0][2*i]-d[2][1][2*i])*s3+
		           (d[0][2][2*i+1]-d[1][2][2*i+1])*c2 - (d[0][2][2*i]-d[1][2][2*i])*s2+
		           (d[2][1][2*i+1]-d[2][0][2*i+1])*c1 - (d[2][1][2*i]-d[2][0][2*i])*s1+
		           (d[1][2][2*i+1]-d[0][2][2*i+1]);

		XX[2*i]   =  2.0*fonfs*(X[2*i]*cSL - X[2*i+1]*sSL);
		XX[2*i+1] =  2.0*fonfs*(X[2*i]*sSL + X[2*i+1]*cSL);
		YY[2*i]   =  2.0*fonfs*(Y[2*i]*cSL - Y[2*i+1]*sSL);
		YY[2*i+1] =  2.0*fonfs*(Y[2*i]*sSL + Y[2*i+1]*cSL);
		ZZ[2*i]   =  2.0*fonfs*(Z[2*i]*cSL - Z[2*i+1]*sSL);
		ZZ[2*i+1] =  2.0*fonfs*(Z[2*i]*sSL + Z[2*i+1]*cSL);
	
		FILE *ptrXSL=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/XSL.bin","ab");
		double tmpRSL=XX[2*i];
		double tmpISL=XX[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrXSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrXSL);
		fclose(ptrXSL);

		FILE *ptrYSL=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/YSL.bin","ab");
		tmpRSL=YY[2*i];
		tmpISL=YY[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrYSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrYSL);
		fclose(ptrYSL);

		FILE *ptrZSL=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/ZSL.bin","ab");
		tmpRSL=ZZ[2*i];
		tmpISL=ZZ[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrZSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrZSL);
		fclose(ptrZSL);

		FILE *ptrX=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/X.bin","ab");
		double tmpR=X[2*i];
		double tmpI=X[2*i+1];
		fwrite(&tmpR, sizeof(double), 1, ptrX);
		fwrite(&tmpI, sizeof(double), 1, ptrX);
		fclose(ptrX);

		FILE *ptrY=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/Y.bin","ab");
		tmpR=Y[2*i];
		tmpI=Y[2*i+1];
		fwrite(&tmpR, sizeof(double), 1, ptrY);
		fwrite(&tmpI, sizeof(double), 1, ptrY);
		fclose(ptrY);

		FILE *ptrZ=fopen("/home/franzi/home/franzi/Documents/01_Code/05_LDC/ldc_waveform_generation/tests/TDI/Z.bin","ab");
		tmpR=Z[2*i];
		tmpI=Z[2*i+1];
		fwrite(&tmpR, sizeof(double), 1, ptrZ);
		fwrite(&tmpI, sizeof(double), 1, ptrZ);
		fclose(ptrZ);
	}

	free(X); free(XX);
	free(Y); free(YY);
	free(Z); free(ZZ);

	return;
}

void init_lisa(struct LISA *lisa, double arm_length, double init_position, double init_rotation){
	
	lisa->number_of_spacecrafts=3;
	lisa->number_of_arms=6;
	lisa->inter_spacecraft_phase= 2 * M_PI / 3;
	lisa->init_position=init_position;
	lisa->init_rotation=init_rotation;
	lisa->spacecraft_separation = arm_length;
	lisa->eccentricity = lisa->spacecraft_separation/(2*sqrt(3)*AU_IN_M);
    lisa->init_time = lisa->init_position * YEAR / (2 * M_PI);
    for (int i = 0; i < lisa->number_of_spacecrafts; i += 1) {
      lisa->rotation[i] = i * lisa->inter_spacecraft_phase + lisa->init_rotation;
      lisa->c_rotation[i] = cos(lisa->rotation[i]);
      lisa->s_rotation[i] = sin(lisa->rotation[i]);
    }
  }

double alpha(double init_time, double time) {
  return BARYCENTER_ANGULAR_VELOCITY * (init_time + time);
} 
    
double AnalyticOrbits_get_position_x(struct LISA *lisa, int sci, double time) {
	double a = alpha(lisa->init_time, time);
	double x;
	x = AU_IN_M * (cos(a) + lisa->eccentricity * (sin(a) * cos(a) * lisa->s_rotation[sci-1] - (1 + sin(a)*sin(a)) * lisa->c_rotation[sci-1]));
	return x;
}

double AnalyticOrbits_get_position_y(struct LISA *lisa, int sci, double time) {
	double a = alpha(lisa->init_time, time);
	double y;
	y = AU_IN_M * (sin(a) + lisa->eccentricity * (sin(a) * cos(a) * lisa->c_rotation[sci-1] - (1 + cos(a)*cos(a)) * lisa->s_rotation[sci-1]));
	return y;
}

double AnalyticOrbits_get_position_z(struct LISA *lisa, int sci, double time) {
	double a = alpha(lisa->init_time, time);
	double z;
	z = -AU_IN_M * lisa->eccentricity * sqrt(3) * cos(a-lisa->rotation[sci-1]);
	return z;
}