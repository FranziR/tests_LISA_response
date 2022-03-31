#define PI 3.1415926535897932384626433832795029L
#define LARM 2.5e9 /* Mean arm length of constellation (m) */
#define ec 0.0048241852 /* LISA orbital eccentricity */
#define KAPPA 0.000000 /* Initial azimuthal position of the guiding center */
#define LAMBDA 0.000000 /* Initial orientation of the LISA constellation */
#define fm 3.168753575e-8 /* LISA modulation frequency */
#define AU_IN_M 1.49597870660e11
#define SIDEREALYEAR_J2000DAY 365.256363004
#define YEAR SIDEREALYEAR_J2000DAY*24*60*60
#define CLIGHT 299792458.0
#define HALF_SCHWARZSCHILD_RADIUS 1.47664E3
#define POST_NEWTONIAN_CONSTANT 1.0
#define BARYCENTER_ANGULAR_VELOCITY (2.0*PI)/(SIDEREALYEAR_J2000DAY*24*60*60)
#define NUM_PARAMETERS_GALBIN 8
#define ECCENTRICITY 0.0048241852
#define INTER_SPACECRAFT_PHASE 2 * M_PI / 3
#define TOBS 33554432

void Fast_GB(double *params, struct Waveform* wfm, long N, double Tobs, double *XSL, double *YSL, double *ZSL, int NP);
void calc_xi_f(struct Waveform *wfm, struct LISA* lisa, double fstar, double t);
void copy_params(struct Waveform *wfm, double *params);
void fill_time_series(struct Waveform *wfm, int n);
void unpack_data(struct Waveform *wfm);
void fft_data(struct Waveform *wfm);
void alloc_waveform(struct Waveform *wfm);
void free_waveform(struct Waveform *wfm);
void calc_d_matrices(struct Waveform *wfm);
void calc_sep_vecs(struct Waveform *wfm, double Larm);
void calc_kdotr(struct Waveform *wfm);
void set_const_trans(struct Waveform *wfm);
void get_basis_tensors(struct Waveform *wfm);	
void get_basis_vecs(double *params, double *u, double *v, double *k);
void get_transfer(struct Waveform *wfm, double t);
void XYZ(double ***d, double f0, long q, long M,  double Tobs, double Larm, double fstar,double *XSL, double *YSL, double *ZSL);
void init_lisa(struct LISA *lisa, double arm_length, double init_position, double init_rotation);
void generateParameters(double *parameters, double freq0, double T);
double alpha(double init_time, double time);
double AnalyticOrbits_get_position_x(struct LISA *lisa, int sci, double time);
double AnalyticOrbits_get_position_y(struct LISA *lisa, int sci, double time);
double AnalyticOrbits_get_position_z(struct LISA *lisa, int sci, double time);
void get_sc_positions(double time, int sci, double *x);

struct LISA
{
	int number_of_spacecrafts;
  	int number_of_arms;
  	double inter_spacecraft_phase;

	double init_time;
	double init_position;
	double init_rotation;
	double eccentricity;
	double spacecraft_separation;

	double rotation[3];
	double c_rotation[3];
	double s_rotation[3];
};

struct Waveform
{
	long N;
	long q; // Fgw carrier bin

	int NP;		// number of parameters in signal

	double T; 		// observation period

	double *params;

	double **eplus, **ecross;
	double **dplus, **dcross;

	double DPr, DPi, DCr, DCi;

	// direction vector of GW
	double *k;

	// separation unit vectors between S/C
	double *r12, *r21;
	double *r13, *r31;
	double *r23, *r32;

	double **kdotr;
	double *kdotx;

	double *xi, *f, *fonfs;

	// Time series of slowly evolving terms at each vertex
	// dataij corresponds to fractional arm length difference yij
	double *data12, *data21;
	double *data13, *data31;
	double *data23, *data32;

	// Fourier coefficients of slowly evolving terms (numerical)
	double *a12, *a21;
	double *a13, *a31;
	double *a23, *a32;

	// S/C position
	double *x, *y, *z;

	// Time varrying quantities (Re & Im) broken up into convenient segments
	double **TR, **TI;

	//Package cij's into proper form for TDI subroutines
	double ***d;
};