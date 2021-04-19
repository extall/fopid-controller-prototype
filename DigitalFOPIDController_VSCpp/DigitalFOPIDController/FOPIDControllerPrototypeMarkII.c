/*
 * FOPIDControllerPrototypeMarkII.c
 *
 * Created: 21.08.2013 10:21:28
 *  Author: Alex
 */

 // Project includes
#include "FOPIDControllerPrototypeMarkII.h"\

// Socket communication related
#include<stdio.h>

// NB! Socket now working on WINDOWS ONLY
#include<winsock2.h>
#include<ws2tcpip.h>

#pragma comment(lib,"ws2_32.lib") //Winsock Library

#pragma warning(disable:4996) 

// ---

#define MATLAB_SERVER "127.0.0.1" //IP address of MATLAB server
#define MATLAB_BUFLEN 128         //Max length of buffer
#define MATLAB_PORT 5110          //The port on which to listen for incoming data
#define MATLAB_FROM_PORT 5101     // Send to here


// Implementation specific defines

// ---------
// Constants
// ---------
#define NMAX  10			// Max order of Oustaloup approximation
#define NSEC  (NMAX+1) 
#define PZMAX (NMAX*2)+1    // Max number of poles (zeros)

#define EPS		 1.0e-10L  // Epsilon
#define NUMMAX	 1.0e+10L

#define M_PI 3.14159265358979323846

// **************
// I/O saturation
// **************
// Input scaling
#define MAX_SCALED_IN     1.0L
#define MIN_SCALED_IN	 -1.0L

#define INPUT_MARGIN      0.6L  // Input change rate limit before
								// sample memory is cleared

// Output scaling
#define MAX_SCALED_OUT    1.0L
#define MIN_SCALED_OUT   -1.0L

// Offset
#define MED_SCALED_IN     0.5L * (MAX_SCALED_IN + MIN_SCALED_IN)
#define MED_SCALED_OUT    0.5L * (MAX_SCALED_OUT + MIN_SCALED_OUT)

// Factor
#define SCALE_FACTOR_IN  (MAX_SCALED_IN  - MED_SCALED_IN)
#define SCALE_FACTOR_OUT (MAX_SCALED_OUT - MED_SCALED_OUT)


// ******* Required sample time *******
#define DT	0.01  		         // [s]
#define SAMPLE_RATE		1 / DT
// ************************************


// ----------------
// GLOBAL variables
// ----------------

// Coefficient generation temporary static storage
volatile double zz[PZMAX] = { 0 };		// Discrete-time zero array
volatile double zp[PZMAX] = { 0 };		// Discrete-time pole array

// IIR SO section coefficient storage for I and D components
volatile double KIc = 0;	// IIR filter gain
volatile double I_zsos[NSEC][2] = { 0 };	// Zero polynomial second-order sections
volatile double I_psos[NSEC][2] = { 0 };	// Pole polynomial second-order sections	
volatile double KDc = 0;
volatile double D_zsos[NSEC][2] = { 0 };
volatile double D_psos[NSEC][2] = { 0 };

// States
volatile double s_I[NSEC][2] = { 0 };
volatile double s_D[NSEC][2] = { 0 };
volatile double s_IntMem = 0;		 // Regular integrator memory
volatile double in_Mem = 0;		 // Previous sample

// FOPID controller parameters: defaults to regular PID controller with direct feed-through
volatile fopid the_fopid = { 1, 0, 0, 1, 1 };

// SOCKET RELATED
volatile int sock = -1;

// The FOPID for which parameters are sought
volatile fopid und_fopid = { 1, 0, 0, 1, 1 };
// Search for 

// Computation flags
volatile int8_t flag_FOPID_Ready = 0;					// FOPID coefficients are being computed
volatile int8_t flag_FOPID_Computing_Output = 0;		// FOPID is still computing the output sample
volatile int8_t flag_FOPID_Schedule_Generation = 0;		// FOPID generation has been scheduled and will
														// take place once the flag_FOPID_Computing_Output
														// will be cleared.

volatile oustapp_params params;							// Approximation parameters

volatile fofopdt_model the_fofopdt;						// The model

volatile design_specs dspecs;							// Design specifications

// Optimization norm
#define OPT_NORM     0.0001L
#define OPT_MAX_ITER 10

// Jacobian and specification function vector are shared
double Jac[3][3] = { 0 };
double kappa_vec[3] = { 0 };

// DEBUG: iterations
uint8_t numIters = 0;

#define ACTIVATE_TUNING 0


void on_exit(void) {
	if (sock != -1) {
		closesocket(sock);
		WSACleanup();
	}
}

// Basic example of a main function
int main(void)
{

	// Define the model
	the_fofopdt.K = 66.16;
	the_fofopdt.L = 1.93;
	the_fofopdt.T = 12.72;
	the_fofopdt.alpha = 0.5;

	the_fopid.Kp = -0.002934;
	the_fopid.Ki = 0.01030;
	the_fopid.Kd = 0.05335;
	the_fopid.lambda = 0.9;
	the_fopid.mu = 0.5;

	// Controller parameters to be tuned
	und_fopid.Kp = 1 / the_fofopdt.K;
	und_fopid.Ki = 1 / the_fofopdt.K;
	und_fopid.Kd = 1 / the_fofopdt.K;
	und_fopid.lambda = 0.9;
	und_fopid.mu = 0.5;

	// Specifications
	dspecs.wc = 0.1;
	dspecs.pm = (60 * M_PI) / 180; // Convert to radians
	dspecs.opt_norm = 0.001; // Optimization termination criterion

	// Set approximation parameters and generate a FOPID controller
	params.wb = 0.001;
	params.wh = 1000;
	params.N = 5;
	params.Ts = DT;

	// If ACTIVATE_TUNING is on, use do the tuning here
	if (ACTIVATE_TUNING) {
		Do_FOPID_Optimization();
	}

	// Generate the FOPID controller
	Generate_FOPID_Controller();

	// We assume UDP based communication with the controller
	// and that the sender of UDP data will take care of the timing
	// This can be considered suboptimal, but for the sake of simplicity
	// and because we are using MATLAB for SIL/HIL experiments, this
	// is good enough for initial tests.

	// Here we are using code from
    // https://www.binarytides.com/udp-socket-programming-in-winsock/
	// as an example implementation of a C UDP client. This is Windows only.
	// For unix based systems, consider this one:
	// https://www.cs.cmu.edu/afs/cs/academic/class/15213-f99/www/class26/udpclient.c

	/* **********************
	// BEGIN UDP CLIENT SETUP
	// ********************** */

	struct sockaddr_in si_other;
	int slen = sizeof(si_other);
	char buf[MATLAB_BUFLEN];
	char message[8];
	WSADATA wsa;

	//Initialise winsock
	printf("\nInitialising Winsock...");
	if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0)
	{
		printf("Failed. Error Code : %d", WSAGetLastError());
		exit(EXIT_FAILURE);
	}
	printf("Initialised.\n");

	//create socket
	if ((sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP)) == SOCKET_ERROR)
	{
		printf("socket() failed with error code : %d", WSAGetLastError());
		exit(EXIT_FAILURE);
	}

	//setup address structure
	memset((char*)&si_other, 0, sizeof(si_other));
	si_other.sin_family = AF_INET;
	si_other.sin_port = htons(MATLAB_PORT);
	si_other.sin_addr.S_un.S_addr = inet_addr(MATLAB_SERVER);

	// Other
	struct sockaddr_in si_out;
	memset((char*)&si_out, 0, sizeof(si_out));
	int olen = sizeof si_out;
	si_out.sin_family = AF_INET;
	si_out.sin_port = htons(MATLAB_FROM_PORT);
	si_out.sin_addr.S_un.S_addr = inet_addr(MATLAB_SERVER);

	// bind
	bind(sock, (SOCKADDR*)&si_other, sizeof(si_other));

	atexit(on_exit); // Need to clean up on exit

	/* ********************
	// END UDP CLIENT SETUP
	// ******************** */

	while (1)
	{
		// TODO: CODE BELOW UNTIL first "memset" is useless right now
		// If a FOPID generation has been scheduled,
		if (flag_FOPID_Schedule_Generation)
		{
			// check whether FOPID output sample is still being computed
			if (!flag_FOPID_Computing_Output)
			{
				// And if not, generate the new controller and clear flag
				Generate_FOPID_Controller();
				flag_FOPID_Schedule_Generation = 0;
			}
		}

		memset(buf, '\0', MATLAB_BUFLEN);
		//try to receive some data, this is a blocking call
		if (recvfrom(sock, buf, MATLAB_BUFLEN, 0, (struct sockaddr*)&si_other, &slen) == SOCKET_ERROR)
		{
			printf("recvfrom() failed with error code : %d", WSAGetLastError());
			exit(EXIT_FAILURE);
		}

		// Convert received value to double (ASSUME it is double and in Little Endian format)
		double err = 0;
		memcpy(&err, buf, sizeof err);

		// Do FOPID control on ERR
		double u = Do_FOPID_Control(err);

		printf("Cur value of u: %f\n", u);

		// Send value back
		memcpy(&message, &u, sizeof u);
		if (sendto(sock, message, sizeof message, 0, (struct sockaddr*)&si_out, olen) == SOCKET_ERROR)
		{
			printf("sendto() failed with error code : %d", WSAGetLastError());
			exit(EXIT_FAILURE);
		}

	}


}


// ********************************
// Coefficient computation function
// ********************************

void Compute_IIR_SOS_Oustaloup(volatile double zCoeffArray[][2], volatile double pCoeffArray[][2],
	volatile double *Kc, volatile oustapp_params params, double alpha)
{
	// Fetch approximation parameter values from params structure
	double T = params.Ts;
	double wb = params.wb;
	double wh = params.wh;
	int8_t N = params.N;

	// Check input: upper frequency bound
	if (wh > (2 / T)) wh = (2 / T);

	double omu = wh / wb;

	int8_t k = 0;		// Iteration variables

	// Calculate discrete-time transfer function zeros and poles
	for (k = -N; k <= N; k++)
	{
		// Continuous-time zeros and poles
		double w_kp = wb * pow(omu, ((k + N + 0.5 - 0.5*alpha) / (2 * N + 1)));
		double w_k = wb * pow(omu, ((k + N + 0.5 + 0.5*alpha) / (2 * N + 1)));

		// Discrete mapping
		zz[k + N] = exp(-T * w_kp);
		zp[k + N] = exp(-T * w_k);
	}

	// Compute center frequency and correct gain
	double wu = sqrt(wb*wh);
	volatile double Ks = pow(wu, alpha);

	// Theta
	volatile double theta = cos(wu*T);

	// Compute absolute value ||H(z)|| at wu rad/s
	volatile double nk, dk, Ku = 1;
	for (k = 0; k < (2 * N + 1); k++)
	{
		nk = (1 - 2 * zz[k] * theta + zz[k] * zz[k]); //nk = sqrt(1-2*zz[k]*theta+zz[k]*zz[k]);
		dk = (1 - 2 * zp[k] * theta + zp[k] * zp[k]); //dk = sqrt(1-2*zp[k]*theta+zp[k]*zp[k]);
		if ((nk > EPS) && (dk > EPS))
		{
			Ku = Ku * (nk / dk);
		}
	}

	Ku = sqrt(Ku);

	// Compute the correct gain
	*Kc = Ks / Ku;

	// Compute second-order section form
	zCoeffArray[0][0] = -zz[2 * N];
	zCoeffArray[0][1] = 0.0;
	pCoeffArray[0][0] = -zp[2 * N];
	pCoeffArray[0][1] = 0.0;
	for (k = N - 1; k >= 0; k--)
	{
		zCoeffArray[N - k][0] = -zz[2 * k + 1] - zz[2 * k];
		zCoeffArray[N - k][1] = zz[2 * k + 1] * zz[2 * k];
		pCoeffArray[N - k][0] = -zp[2 * k + 1] - zp[2 * k];
		pCoeffArray[N - k][1] = zp[2 * k + 1] * zp[2 * k];
	}

	// Check stability
	volatile int isStable = IIR_SOS_Stability_Test(pCoeffArray);
}

void Generate_FOPID_Controller()
{
	// Reset FOPID_Ready flag
	flag_FOPID_Ready = 0;

	// Clear the memory banks
	Clear_IIR_Memory();

	// Generate the integrator
	Compute_IIR_SOS_Oustaloup(I_zsos, I_psos, &KIc, params, 1.0L - the_fopid.lambda);

	// Generate the differentiator
	Compute_IIR_SOS_Oustaloup(D_zsos, D_psos, &KDc, params, the_fopid.mu);

	// All done: Set FOPID_Ready flag
	flag_FOPID_Ready = 1;
}

void Clear_IIR_Memory()
{
	for (uint8_t k = 0; k < NSEC; k++)
	{
		s_I[k][0] = 0; s_I[k][1] = 0;
		s_D[k][0] = 0; s_D[k][1] = 0;
	}
	s_IntMem = 0;
}

int IIR_SOS_Stability_Test(volatile double pCoeffArray[][2])
{
	// Use the "triangle rule" for SOS poles
	volatile double d1, d2, d1_abs, d2_abs, d1_abs_m_1, d2_p_1;
	for (uint8_t i = 0; i <= params.N; i++)
	{
		d1 = pCoeffArray[i][0];
		d2 = pCoeffArray[i][1];
		d1_abs = fabs(d1);
		d2_abs = fabs(d2);
		d1_abs_m_1 = d1_abs - 1.0L;
		d2_p_1 = d2 + 1.0L;

		if (!(fabs(d1) < d2 + 1.0L)) return 0; // Condition (1)
		if (!(d2_abs < 1)) return 0;           // Condition (2)
	}

	// All coefficients are stable
	return 1;
}

// ***********************
// IIR filtering functions
// ***********************
double Do_IIR_Filtering(volatile double zCoeffArray[NSEC][2],
	volatile double pCoeffArray[NSEC][2],
	volatile double Kc, volatile double memArray[][2], double input)
{
	// Assign local pointers
	volatile double(*b)[2] = zCoeffArray;
	volatile double(*a)[2] = pCoeffArray;
	volatile double(*s)[2] = memArray;

	// Local input/output
	volatile double u_n, y_n;

	// Filter signal
	u_n = Kc * input;
	for (uint8_t m = 0; m < params.N + 1; m++)
	{
		y_n = u_n + s[m][0];
		s[m][0] = b[m][0] * u_n - a[m][0] * y_n + s[m][1];
		s[m][1] = b[m][1] * u_n - a[m][1] * y_n;

		u_n = y_n;
	}

	// Check memory for under/overflow
	// NB! TODO: why do we do BOTH checks here? Should actually be part of DO_FOPID_Control function
	for (uint8_t k = 0; k < NSEC; k++)
	{
		s_I[k][0] = double_scale_saturation(s_I[k][0], EPS, NUMMAX);
		s_I[k][1] = double_scale_saturation(s_I[k][1], EPS, NUMMAX);
		s_D[k][0] = double_scale_saturation(s_D[k][0], EPS, NUMMAX);
		s_D[k][1] = double_scale_saturation(s_D[k][1], EPS, NUMMAX);
	}
	s_IntMem = double_scale_saturation(s_IntMem, EPS, NUMMAX);

	// Assign output
	return u_n;
}

double Do_FOPID_Control(double err)
{
	// Begin computing the output sample
	flag_FOPID_Computing_Output = 1;

	volatile double in, out;
	volatile double i_out, foi_out, fod_out;

	// Get the scaled ADC value
	in = err;

	// Check the input margin; clear memory if exceeded
	if (fabs(in - in_Mem) > INPUT_MARGIN)
	{
		Clear_IIR_Memory();
	}
	in_Mem = in;

	// FOPID computation for this sample
	foi_out = the_fopid.Ki*Do_IIR_Filtering(I_zsos, I_psos, KIc, s_I, in);
	fod_out = the_fopid.Kd*Do_IIR_Filtering(D_zsos, D_psos, KDc, s_D, in);
	s_IntMem += params.Ts * foi_out;
	i_out = s_IntMem;
	out = the_fopid.Kp*in + i_out + fod_out;

	// Done computing the output sample
	flag_FOPID_Computing_Output = 0;

	// Set the scaled value to DAC
	return Set_Scaled_Output(out);
}


double Set_Scaled_Output(double out)
{
	// Saturate at top values
	volatile double scaled_out = out;
	if (out > MAX_SCALED_OUT) scaled_out = MAX_SCALED_OUT;
	if (out < MIN_SCALED_OUT) scaled_out = MIN_SCALED_OUT;

	// Transform to normalized range -1...1
	scaled_out = (scaled_out - MED_SCALED_OUT) / (SCALE_FACTOR_OUT);

	// Return the scaled output
	return scaled_out;
}


// Additional mathematical functions
double sign(double x)
{
	return (x > 0) - (x < 0);
}

double double_scale_saturation(double x, double min, double max)
{
	double y = x;
	if (fabs(x) < min) y = 0;
	if (fabs(x) > max) y = sign(x) * max;
	return y;
}

// Helper trigonometric functions
double _sf(double x) {
	return sin((M_PI * x) / 2);
}

double _cf(double x) {
	return cos((M_PI * x) / 2);
}

// *****************************************
// Magnitude and phase response computations
// *****************************************

// Plant magnitude
double magng(double w) {
	return fabs(the_fofopdt.K) /
		(sqrt(1 + pow(the_fofopdt.T, 2)*pow(w, 2 * the_fofopdt.alpha)
			+ 2 * the_fofopdt.T*pow(w, the_fofopdt.alpha)*_cf(the_fofopdt.alpha)));
};

// Controller magnitude
double magnfopid(double w) {
	volatile double CR = und_fopid.Kp + pow(w, -und_fopid.lambda)*und_fopid.Ki*_cf(und_fopid.lambda)
		+ pow(w, und_fopid.mu) * und_fopid.Kd*_cf(und_fopid.mu);

	volatile double CI = -pow(w, -und_fopid.lambda)*und_fopid.Ki*_sf(und_fopid.lambda)
		+ pow(w, und_fopid.mu)*und_fopid.Kd*_sf(und_fopid.mu);

	return sqrt(pow(CR, 2) + pow(CI, 2));
};

// Plant phase
double phg(double w) {
	return -the_fofopdt.L*w - atan(((the_fofopdt.T)*_sf(the_fofopdt.alpha)) /
		(pow(w, -the_fofopdt.alpha) + the_fofopdt.T*_cf(the_fofopdt.alpha)));
};

// Controller phase
double phfopid(double w) {

	volatile double CN = pow(w, und_fopid.lambda + und_fopid.mu)*und_fopid.Kd*_sf(und_fopid.mu)
		- und_fopid.Ki * _sf(und_fopid.lambda);

	volatile double CD = und_fopid.Ki*_cf(und_fopid.lambda) + pow(w, und_fopid.lambda)
		* (pow(w, und_fopid.mu)*und_fopid.Kd*_cf(und_fopid.mu) + und_fopid.Kp);

	return atan(CN / CD);
};

// The PSI functions and their derivatives

// Phase margin
double psi_pm(double w) {
	return (magng(w) * magnfopid(w) - 1);
};

double dpsi_pm(double w) {

	// System magnitude responses at w
	double GM = magng(w);
	double CM = magnfopid(w);

	// Helper values
	double A11 = pow(w, -1 - 2 * und_fopid.lambda)*
		(
			und_fopid.mu*pow(w, 2 * (und_fopid.lambda + und_fopid.mu))*pow(und_fopid.Kd, 2) - und_fopid.lambda*und_fopid.Ki *
			(
				und_fopid.Ki + pow(w, und_fopid.lambda)*und_fopid.Kp*_cf(und_fopid.lambda)
				)
			+ pow(w, und_fopid.lambda + und_fopid.mu) * und_fopid.Kd *
			(
			(und_fopid.lambda - und_fopid.mu)*und_fopid.Ki*_sf(und_fopid.lambda + und_fopid.mu)
				+ und_fopid.mu * pow(w, und_fopid.lambda) * und_fopid.Kp * _cf(und_fopid.mu)
				)
			);

	double A1d = A11 / CM;

	double A2d = -((
		the_fofopdt.T*the_fofopdt.alpha*pow(w, the_fofopdt.alpha - 1)*
		(
			the_fofopdt.T * pow(w, the_fofopdt.alpha) + _cf(the_fofopdt.alpha)
			)
		) * GM) /
		(1 + pow(the_fofopdt.T, 2)*pow(w, 2 * the_fofopdt.alpha)
			+ 2 * the_fofopdt.T*pow(w, the_fofopdt.alpha)*_cf(the_fofopdt.alpha));

	return A1d * GM + CM * A2d;

};

// Gain margin
double psi_gm(double w) {
	return phg(w) + phfopid(w) + M_PI;
};

double dpsi_gm(double w) {

	double vd = dphg(w);
	double B1d = dphfopid(w);

	return B1d + vd;

};

// Simple Newton's method test
double NRM_simple(double w0, double(*f)(double), double(*df)(double))
{
	// Parameters
	uint8_t N = 25;
	uint8_t k = 0;

	volatile double x = w0;
	volatile double xo = x;
	volatile double NRM_eps = 0.001;
	volatile double gamma = 1.5;

	volatile double f_val = f(x);
	volatile double df_val = df(x);

	while (k++<N && f_val>NRM_eps)
	{
		f_val = f(x);
		df_val = df(x);
		xo = x;
		x = x - f_val / df_val;
		if (x < 0) {
			x = xo * gamma;
		}
	}

	return x;
}

// Derivative of plant phase response 
double dphg(double w)
{
	return -the_fofopdt.L -
		(the_fofopdt.alpha*the_fofopdt.T*_sf(the_fofopdt.alpha)) /
		(w*(
			2 * the_fofopdt.T*_cf(the_fofopdt.alpha) + pow(w, -the_fofopdt.alpha)
			+ pow(the_fofopdt.T, 2)*pow(w, the_fofopdt.alpha)
			));
}

// Derivative of controller phase response
double dphfopid(double w)
{
	double B11 = (und_fopid.lambda + und_fopid.mu)*und_fopid.Ki*_cf(und_fopid.lambda + und_fopid.mu - 1)
		+ und_fopid.mu * pow(w, und_fopid.lambda) * und_fopid.Kp * _sf(und_fopid.mu);

	double B20 = w * (
		pow(w, und_fopid.lambda + 2 * und_fopid.mu)*pow(und_fopid.Kd, 2) + pow(w, -und_fopid.lambda)*pow(und_fopid.Ki, 2)
		+ 2 * und_fopid.Kp*und_fopid.Ki*_cf(und_fopid.lambda)
		+ pow(w, und_fopid.lambda) * pow(und_fopid.Kp, 2) - 2 * pow(w, und_fopid.mu)*und_fopid.Kd*(
			und_fopid.Ki*_sf(und_fopid.lambda + und_fopid.mu - 1) - pow(w, und_fopid.lambda) * und_fopid.Kp * _cf(und_fopid.mu)
			)
		);

	double B10 = und_fopid.lambda * und_fopid.Kp * und_fopid.Ki * _sf(und_fopid.lambda) +
		pow(w, und_fopid.mu) * und_fopid.Kd * B11;

	return B10 / B20;
}

// FOPID optimization based on frequency-domain specifications
void Do_FOPID_Optimization(void)
{
	uint8_t k;

	for (k = 0; k < OPT_MAX_ITER; k++)
	{
		numIters = k;									// Number of used iterations
		double b[3] = { -kappa1(), -kappa2(), -kappa3() };	// Cost functions
		compute_specs_J();								// Update Jacobian
		sle_sol dx = compute_cramer3(Jac, b);			// Solve the system of equations
		double xn[3] = { und_fopid.Kp + dx.x1,			// Update the solution
					   und_fopid.Ki + dx.x2,
					   und_fopid.Kd + dx.x3 };
		und_fopid.Kp = xn[0];							// Set the new solution
		und_fopid.Ki = xn[1];
		und_fopid.Kd = xn[2];

		// Compute norm and check it
		if (norm2_v3(b) < OPT_NORM) return;
	}
}

// Critical frequency specification
double kappa1(void)
{
	return magnfopid(dspecs.wc) * magng(dspecs.wc) - 1;
}

// Phase margin specification
double kappa2(void)
{
	return phfopid(dspecs.wc) + phg(dspecs.wc) + M_PI - dspecs.pm;
}

// Phase response flatness at wc
double kappa3(void)
{
	return dphg(dspecs.wc) + dphfopid(dspecs.wc);
}

// Jacobian computation
void compute_specs_J(void)
{
	// Enhance readability
	// |
	// Controller parameters
	double Kp = und_fopid.Kp;
	double Ki = und_fopid.Ki;
	double Kd = und_fopid.Kd;
	double lam = und_fopid.lambda;
	double mu = und_fopid.mu;
	// |
	// Critical frequency
	double w = dspecs.wc;

	double A12 = pow(w, -2 * lam)*(-pow(w, lam + mu)*_sf(lam + mu - 1)*Kd + Ki + pow(w, lam)*_cf(lam)*Kp);
	double A13 = pow(w, -lam + mu)*(pow(w, lam + mu)*Kd - _sf(lam + mu - 1)*Ki + pow(w, lam)*_cf(mu)*Kp);
	double A2 = pow(w, 2 * (lam + mu))*pow(Kd, 2) + pow(Ki, 2) + 2 * pow(w, lam)*_cf(lam)*Ki*Kp + pow(w, 2 * lam)*pow(Kp, 2)
		+ 2 * pow(w, lam + mu)*Kd*(-_sf(lam + mu - 1)*Ki + pow(w, lam)*_cf(mu)*Kp);
	double A3 = pow(A2, 2);

	double A31 = pow(w, lam - 1)*(
		mu*pow(w, 3 * (lam + mu))*_sf(mu)*pow(Kd, 3) - pow(w, 2 * (lam + mu))*(2 * mu*_sf(lam)
			+ lam * _sf(lam + 2 * mu))*pow(Kd, 2)*Ki + lam * _sf(lam)*Ki*(pow(Ki, 2) - pow(w, 2 * lam)*pow(Kp, 2)) - pow(w, lam + mu)*Kd*
			((2 * lam*_sf(mu) + mu * _sf(2 * lam + mu))*pow(Ki, 2) + 2 * (lam + mu)*pow(w, lam)*_sf(lam + mu - 1)*Ki*Kp + mu * pow(w, 2 * lam)*_sf(mu)*pow(Kp, 2)));

	double A32 = pow(w, lam - 1)*(
		(lam + mu)*pow(w, 2 * lam + 3 * mu)*_cf(lam + mu - 1)*pow(Kd, 3) + pow(w, 2 * (lam + mu))*
		(2 * (lam + mu)*_sf(lam) + lam * _sf(lam + 2 * mu))*pow(Kd, 2)*Kp +
		lam * _sf(lam)*Kp*(-pow(Ki, 2) + pow(w, 2 * lam)*pow(Kp, 2)) + pow(w, mu)*Kd*
		(-(lam + mu)*_cf(mu + lam - 1)*pow(Ki, 2) - 2 * mu*pow(w, lam)*_sf(mu)*Ki*Kp
			+ pow(w, 2 * lam)*(2 * lam*_cf(lam + mu - 1) + (lam + mu)*_sf(lam - mu))*pow(Kp, 2))
		);

	double A33 = pow(w, lam + mu - 1)*(
		(lam + mu)*_cf(lam + mu - 1)*pow(Ki, 3) - 2 * lam*pow(w, 2 * lam + mu)*
		_sf(lam)*Kd*Ki*Kp + pow(w, lam) * (2 * (lam + mu)*_sf(mu) + mu * _sf(2 * lam + mu))*
		pow(Ki, 2)*Kp + pow(w, 2 * lam)*(2 * mu*_cf(lam + mu - 1) - (lam + mu)*_sf(lam - mu))*Ki*pow(Kp, 2)
		+ mu * pow(w, 3 * lam)*_sf(mu)*pow(Kp, 3) - pow(w, 2 * (lam + mu))*pow(Kd, 2)*((lam + mu)*_cf(lam + mu - 1)*Ki
			+ mu * pow(w, lam)*_sf(mu)*Kp));

	double ACR = Kp + pow(w, -lam)*_cf(lam)*Ki + pow(w, mu)*_cf(mu)*Kd;

	double magndiv = magng(w) / magnfopid(w);

	Jac[0][0] = magndiv * ACR;
	Jac[0][1] = magndiv * A12;
	Jac[0][2] = magndiv * A13;
	Jac[1][0] = (pow(w, lam)*(-pow(w, lam + mu)*_sf(mu)*Kd + _sf(lam)*Ki)) / A2;
	Jac[1][1] = -(pow(w, lam)*(pow(w, mu)*_cf(mu + lam - 1)*Kd + _sf(lam)*Kp)) / A2;
	Jac[1][2] = (pow(w, lam + mu)*(_cf(lam + mu - 1)*Ki + pow(w, lam)*_sf(mu)*Kp)) / A2;
	Jac[2][0] = A31 / A3;
	Jac[2][1] = A32 / A3;
	Jac[2][2] = A33 / A3;
}

// 3x3 matrix: determinant computation
double compute_det3(double A[3][3])
{
	return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
		A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
		A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

// 3x3 matrix: Cramer rule
sle_sol compute_cramer3(double A[3][3], double b[3])
{
	// EPS
	double my_eps = 1e-7;

	// The solution
	sle_sol system_solution;
	system_solution.exists = false;

	// Compute the determinant of A and check it
	double detA = compute_det3(A);
	if (fabs(detA) < my_eps) {
		return system_solution;
	}

	// Solve the system
	double mat_x1[3][3] = { {b[0],A[0][1],A[0][2]},{b[1],A[1][1],A[1][2]},{b[2],A[2][1],A[2][2]} };
	double mat_x2[3][3] = { {A[0][0],b[0],A[0][2]},{A[1][0],b[1],A[1][2]},{A[2][0],b[2],A[2][2]} };
	double mat_x3[3][3] = { {A[0][0],A[0][1],b[0]},{A[1][0],A[1][1],b[1]},{A[2][0],A[2][1],b[2]} };

	system_solution.x1 = compute_det3(mat_x1) / detA;
	system_solution.x2 = compute_det3(mat_x2) / detA;
	system_solution.x3 = compute_det3(mat_x3) / detA;
	system_solution.exists = true;

	return system_solution;
}


