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
	int n_waves; 
	n_waves=(int)atoi(argv[1]);

	int N =1024;
	
	for(int nn=0; nn<n_waves; nn++){
		printf("n_waves= %d\n ", nn);

		double *parameters, *XFFT, *YFFT, *ZFFT;
		parameters = (double *)malloc(NUM_PARAMETERS_GALBIN*sizeof(double));
		XFFT = (double *)malloc(2*N*sizeof(double)); YFFT = (double *)malloc(2*N*sizeof(double)); ZFFT = (double *)malloc(2*N*sizeof(double));
		struct Waveform *wfm = (struct Waveform*) malloc(sizeof(struct Waveform));
		
		generateParameters(parameters, 1e-3, TOBS);

		Fast_GB(parameters, wfm, N, TOBS, XFFT, YFFT, ZFFT, NUM_PARAMETERS_GALBIN);

		free(parameters); free(XFFT); free(YFFT); free(ZFFT);
	}

    return 0;
}

void generateParameters(double *parameters, double freq0, double T){
		double datatmp[NUM_PARAMETERS_GALBIN];
		datatmp[0]=2*PI*(double)(rand()/ (double)(RAND_MAX))-(PI);
		datatmp[1]=PI*(double)(rand()/ (double)(RAND_MAX))-PI/2;
		datatmp[2]=2*PI*(double)(rand()/ (double)(RAND_MAX));
		datatmp[3]=2*(double)(rand()/ (double)(RAND_MAX));
		datatmp[4]=2*PI*(double)(rand()/ (double)(RAND_MAX));
		datatmp[5]=freq0+2e-6*(double)(rand()/ (double)(RAND_MAX))-1e-6;
		datatmp[6]=2*PI*(double)(rand()/ (double)(RAND_MAX));
		datatmp[7]=2e-19*(double)(rand()/ (double)(RAND_MAX))-1e-19;

		parameters[0]=datatmp[5]*T; // frequency
		parameters[7]=datatmp[7]*T*T; // frequency dot
		parameters[1]=cos(PI/2-datatmp[1]); // cos(theta) i.e. cos(pi/2-beta)
		parameters[2]=datatmp[0]; // phi i.e. lambda
		parameters[3]=log(datatmp[3]); // amplitude
		parameters[4]=cos(datatmp[4]); // cos(iota)
		parameters[5]=datatmp[2]; // psi
		parameters[6]=-1*datatmp[6]; // phi_0

		FILE *ptrpara=fopen("matlab/parameters.bin","ab");
		for(int ii=0; ii<NUM_PARAMETERS_GALBIN; ii++){
			fwrite(&datatmp[ii], sizeof(double), 1, ptrpara);
		}
		fclose(ptrpara);
}

void Fast_GB(double *params, struct Waveform* wfm, long N, double T, double *XSL, double *YSL, double *ZSL, int NP)
{		
	long n; // iterator
	double t; // time
	double* orbit_params = (double*) malloc(3*sizeof(double));
  	orbit_params[0] = LARM; orbit_params[1] = LAMBDA; orbit_params[2] = KAPPA;

	struct LISA *lisa = (struct LISA*) malloc(sizeof(struct LISA));
	init_lisa(lisa, orbit_params[0], orbit_params[1], orbit_params[2]);
	double fstar = (CLIGHT/lisa->spacecraft_separation)/(2*PI);

	wfm->N  = N; // set number of samples
	wfm->T  = T; // set observation period
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

	XYZ(wfm->d, wfm->params[0]/wfm->T, wfm->q, N, wfm->T, lisa->spacecraft_separation, fstar, XSL, YSL, ZSL);

	free_waveform(wfm); // Deallocate memory
	free(wfm);

	return;
}

void calc_xi_f(struct Waveform *wfm, struct LISA* lisa, double fstar, double t)
{
	long i;

	double f0, dfdt_0;

	f0       = wfm->params[0]/wfm->T;
	dfdt_0   = wfm->params[7]/wfm->T/wfm->T;

	// Calculate position of each spacecraft at time t
	for(i=0; i<3; i++){
	  wfm->x[i] = AnalyticOrbits_get_position_x(lisa, i+1, t);
	  wfm->y[i] = AnalyticOrbits_get_position_y(lisa, i+1, t);
	  wfm->z[i] = AnalyticOrbits_get_position_z(lisa, i+1, t);
	}

	for(i=0; i<3; i++)
	{
		wfm->kdotx[i] = (wfm->x[i]*wfm->k[0] + wfm->y[i]*wfm->k[1] + wfm->z[i]*wfm->k[2])/CLIGHT;
		wfm->xi[i]    = t - wfm->kdotx[i]; //Wave arrival time at spacecraft i
		wfm->f[i]     = f0+dfdt_0*wfm->xi[i];; //First order approximation to frequency at spacecraft i
		wfm->fonfs[i] = wfm->f[i]/fstar; //Ratio of true frequency to transfer frequency
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
	double f0, dfdt_0, df, phi0;

	f0       = wfm->params[0]/wfm->T;
	dfdt_0   = wfm->params[7]/wfm->T/wfm->T;
	phi0     = wfm->params[6];

	q  = wfm->q;
	df = PI*2*(((double)q)/wfm->T);

	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if(i!=j)
			{
				
				arg1 = 0.5*wfm->fonfs[i]*(1. + wfm->kdotr[i][j]); //Argument of transfer function
				arg2 =  PI*2*f0*wfm->xi[i] + phi0 - df*t+M_PI*dfdt_0*wfm->xi[i]*wfm->xi[i]; //Argument of complex exponentials
				sinc = 0.25*sin(arg1)/arg1; //Transfer function
				aevol = 1.0; //Evolution of amplitude
				tran1r = aevol*(wfm->dplus[i][j]*wfm->DPr + wfm->dcross[i][j]*wfm->DCr); ///Real and imaginary pieces of time series (no complex exponential)
				tran1i = aevol*(wfm->dplus[i][j]*wfm->DPi + wfm->dcross[i][j]*wfm->DCi);
				tran2r = cos(arg1 + arg2); //Real and imaginry components of complex exponential
				tran2i = sin(arg1 + arg2);
				wfm->TR[i][j] = sinc*(tran1r*tran2r - tran1i*tran2i); //Real & Imaginary part of the slowly evolving signal
				wfm->TI[i][j] = sinc*(tran1r*tran2i + tran1i*tran2r);
			}
		}
	}

	return;
}

void XYZ(double ***d, double f0, long q, long M, double T, double Larm, double fstar,double *XSL, double *YSL, double *ZSL)
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
		f = ((double)(q + i - M/2))/T;
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
	
		FILE *ptrXSL=fopen("matlab/X.bin","ab");
		double tmpRSL=XX[2*i];
		double tmpISL=XX[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrXSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrXSL);
		fclose(ptrXSL);

		FILE *ptrYSL=fopen("matlab/Y.bin","ab");
		tmpRSL=YY[2*i];
		tmpISL=YY[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrYSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrYSL);
		fclose(ptrYSL);

		FILE *ptrZSL=fopen("matlab/Z.bin","ab");
		tmpRSL=ZZ[2*i];
		tmpISL=ZZ[2*i+1];
		fwrite(&tmpRSL, sizeof(double), 1, ptrZSL);
		fwrite(&tmpISL, sizeof(double), 1, ptrZSL);
		fclose(ptrZSL);
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