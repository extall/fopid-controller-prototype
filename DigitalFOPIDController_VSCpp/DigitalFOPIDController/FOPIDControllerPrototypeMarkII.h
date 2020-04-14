/*
 * FOPIDControllerPrototypeMarkII.h
 *
 * Created: 21.08.2013 10:22:19
 *  Author: Alex
 */

 // Project standard library includes
#include <stdint.h>
#include <stdbool.h>
#include <math.h>



// ****************************
// Custom structure definitions
// ****************************

// IO -----------------------------------------------

typedef struct T_reading
{
	double T_tc;
	double T_internal;
	unsigned int state;
} T_reading;

// --------------------------------------------------

// Oustaloup's approximation parameters
typedef struct oustapp_params
{
	double wb;
	double wh;
	int8_t N;
	double Ts;
} oustapp_params;

// FO PID controller parameters
typedef struct fopid
{
	double Kp;
	double Ki;
	double Kd;
	double lambda;
	double mu;
} fopid;

// FO FOPDT model parameters
typedef struct fofopdt_model
{
	double K;
	double L;
	double T;
	double alpha;
} fofopdt_model;

// FO control system tuning parameters
typedef struct  design_specs
{
	double wc;
	double pm;
	double opt_norm;
} design_specs;

// Solution to a system of linear equations (3)
typedef struct sle_sol
{
	double x1;
	double x2;
	double x3;
	bool exists;
} sle_sol;



// Function prototypes

// **************************************
// Filtering and control module functions
// **************************************
void Compute_IIR_SOS_Oustaloup(volatile double zCoeffArray[][2], volatile double pCoeffArray[][2],
	volatile double *Kc, oustapp_params params, volatile double alpha);
int IIR_SOS_Stability_Test(volatile double pCoeffArray[][2]);
void Generate_FOPID_Controller(void);
void Clear_IIR_Memory(void);
double Do_IIR_Filtering(volatile double zCoeffArray[][2], volatile double pCoeffArray[][2],
	volatile double Kc, volatile double memArray[][2], double input);
double Do_FOPID_Control(double err);

void Generate_FOPID_Controller(void);
int IIR_SOS_Stability_Test(double arr[][2]);
void Clear_IIR_Memory();
double Set_Scaled_Output(double out);

// *********************************
// Additional mathematical functions
// *********************************
double sign(double x);
double double_scale_saturation(double x, double min, double max);

// Helper functions
double _sf(double x);
double _cf(double x);

// *****************************************
// Magnitude and phase response computations
// *****************************************
double magng(double w);
double magnfopid(double w);
double phg(double w);
double phfopid(double w);

// Derivatives
double dphg(double w);
double dphfopid(double w);

// The PSI functions
double psi_pm(double w);
double dpsi_pm(double w);
double psi_gm(double w);
double dpsi_gm(double w);

// Function for optimization
void Do_FOPID_Optimization(void);

// Functions for specifications
double kappa1(void);
double kappa2(void);
double kappa3(void);

// Function for Jacobian
void compute_specs_J(void);

// 3x3 matrix: determinant computation
double compute_det3(double A[3][3]);

// 3x3 matrix: Cramer rule
sle_sol compute_cramer3(double A[3][3], double b[3]);

// Vector norm
double norm2_v3(double b[3])
{
	return sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
}


// Newtons Method
double NRM_simple(double w0, double(*f)(), double(*df)());



