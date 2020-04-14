/*
 * rt_test.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "rt_test".
 *
 * Model version              : 1.8
 * Simulink Coder version : 9.0 (R2018b) 24-May-2018
 * C source code generated on : Mon Apr 13 12:36:53 2020
 *
 * Target selection: sldrt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "rt_test.h"
#include "rt_test_private.h"
#include "rt_test_dt.h"

/* options for Simulink Desktop Real-Time board 0 */
static double SLDRTBoardOptions0[] = {
  0.0,
  5110.0,
  49.0,
  50.0,
  55.0,
  46.0,
  48.0,
  46.0,
  48.0,
  46.0,
  49.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};

/* options for Simulink Desktop Real-Time board 1 */
static double SLDRTBoardOptions1[] = {
  0.0,
  5111.0,
  49.0,
  50.0,
  55.0,
  46.0,
  48.0,
  46.0,
  48.0,
  46.0,
  49.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
};

/* list of Simulink Desktop Real-Time timers */
const int SLDRTTimerCount = 1;
const double SLDRTTimers[2] = {
  0.01, 0.0,
};

/* list of Simulink Desktop Real-Time boards */
const int SLDRTBoardCount = 2;
SLDRTBOARD SLDRTBoards[2] = {
  { "Standard_Devices/UDP_Protocol", 5100U, 256, SLDRTBoardOptions0 },

  { "Standard_Devices/UDP_Protocol", 5101U, 256, SLDRTBoardOptions1 },
};

/* Block signals (default storage) */
B_rt_test_T rt_test_B;

/* Continuous states */
X_rt_test_T rt_test_X;

/* Block states (default storage) */
DW_rt_test_T rt_test_DW;

/* Real-time model */
RT_MODEL_rt_test_T rt_test_M_;
RT_MODEL_rt_test_T *const rt_test_M = &rt_test_M_;

/*
 * This function updates continuous states using the ODE5 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE5_A[6] = {
    1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0
  };

  static const real_T rt_ODE5_B[6][6] = {
    { 1.0/5.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    { 3.0/40.0, 9.0/40.0, 0.0, 0.0, 0.0, 0.0 },

    { 44.0/45.0, -56.0/15.0, 32.0/9.0, 0.0, 0.0, 0.0 },

    { 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0, 0.0 },

    { 9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0,
      0.0 },

    { 35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE5_IntgData *id = (ODE5_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T *f3 = id->f[3];
  real_T *f4 = id->f[4];
  real_T *f5 = id->f[5];
  real_T hB[6];
  int_T i;
  int_T nXc = 2;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  rt_test_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE5_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE5_A[0]);
  rtsiSetdX(si, f1);
  rt_test_output();
  rt_test_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE5_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE5_A[1]);
  rtsiSetdX(si, f2);
  rt_test_output();
  rt_test_derivatives();

  /* f(:,4) = feval(odefile, t + hA(3), y + f*hB(:,3), args(:)(*)); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE5_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, t + h*rt_ODE5_A[2]);
  rtsiSetdX(si, f3);
  rt_test_output();
  rt_test_derivatives();

  /* f(:,5) = feval(odefile, t + hA(4), y + f*hB(:,4), args(:)(*)); */
  for (i = 0; i <= 3; i++) {
    hB[i] = h * rt_ODE5_B[3][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2] +
                   f3[i]*hB[3]);
  }

  rtsiSetT(si, t + h*rt_ODE5_A[3]);
  rtsiSetdX(si, f4);
  rt_test_output();
  rt_test_derivatives();

  /* f(:,6) = feval(odefile, t + hA(5), y + f*hB(:,5), args(:)(*)); */
  for (i = 0; i <= 4; i++) {
    hB[i] = h * rt_ODE5_B[4][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2] +
                   f3[i]*hB[3] + f4[i]*hB[4]);
  }

  rtsiSetT(si, tnew);
  rtsiSetdX(si, f5);
  rt_test_output();
  rt_test_derivatives();

  /* tnew = t + hA(6);
     ynew = y + f*hB(:,6); */
  for (i = 0; i <= 5; i++) {
    hB[i] = h * rt_ODE5_B[5][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2] +
                   f3[i]*hB[3] + f4[i]*hB[4] + f5[i]*hB[5]);
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void rt_test_output(void)
{
  if (rtmIsMajorTimeStep(rt_test_M)) {
    /* set solver stop time */
    if (!(rt_test_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&rt_test_M->solverInfo,
                            ((rt_test_M->Timing.clockTickH0 + 1) *
        rt_test_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&rt_test_M->solverInfo,
                            ((rt_test_M->Timing.clockTick0 + 1) *
        rt_test_M->Timing.stepSize0 + rt_test_M->Timing.clockTickH0 *
        rt_test_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(rt_test_M)) {
    rt_test_M->Timing.t[0] = rtsiGetT(&rt_test_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(rt_test_M)) {
    /* S-Function (sldrtpo): '<Root>/Packet Output' */
    /* S-Function Block: <Root>/Packet Output */

    /* no code required */
  }

  /* TransferFcn: '<Root>/Transfer Fcn1' */
  rt_test_B.TransferFcn1 = 0.0;
  rt_test_B.TransferFcn1 += rt_test_P.TransferFcn1_C *
    rt_test_X.TransferFcn1_CSTATE;

  /* TransferFcn: '<Root>/Transfer Fcn' */
  rt_test_B.TransferFcn = 0.0;
  rt_test_B.TransferFcn += rt_test_P.TransferFcn_C *
    rt_test_X.TransferFcn_CSTATE;
  if (rtmIsMajorTimeStep(rt_test_M)) {
  }

  /* Sum: '<Root>/Sum' incorporates:
   *  Constant: '<Root>/Constant'
   */
  rt_test_B.Sum = rt_test_P.Constant_Value - rt_test_B.TransferFcn1;
  if (rtmIsMajorTimeStep(rt_test_M)) {
    /* DiscreteZeroPole: '<S1>/Discrete Zero-Pole' */
    {
      {
        static const int_T colCidxRow0[23] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 };

        const int_T *pCidx = &colCidxRow0[0];
        const real_T *pC0 = rt_test_P.DiscreteZeroPole_C;
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *y0 = &rt_test_B.DiscreteZeroPole;
        int_T numNonZero = 22;
        *y0 = (*pC0++) * xd[*pCidx++];
        while (numNonZero--) {
          *y0 += (*pC0++) * xd[*pCidx++];
        }
      }

      rt_test_B.DiscreteZeroPole += rt_test_P.DiscreteZeroPole_D*rt_test_B.Sum;
    }

    /* S-Function (sldrtpi): '<Root>/Packet Input' */
    /* S-Function Block: <Root>/Packet Input */
    {
      uint8_T indata[8U];
      int status = RTBIO_DriverIO(1, STREAMINPUT, IOREAD, 8U,
        &rt_test_P.PacketInput_PacketID, (double*) indata, NULL);
      if (status & 0x1) {
        RTWin_ANYTYPEPTR indp;
        indp.p_uint8_T = indata;
        rt_test_B.PacketInput = *indp.p_real_T++;
      }
    }
  }

  /* Sum: '<Root>/Sum1' incorporates:
   *  Constant: '<Root>/Constant'
   */
  rt_test_B.Sum1 = rt_test_P.Constant_Value - rt_test_B.TransferFcn;
}

/* Model update function */
void rt_test_update(void)
{
  if (rtmIsMajorTimeStep(rt_test_M)) {
    /* Update for S-Function (sldrtpo): '<Root>/Packet Output' */

    /* S-Function Block: <Root>/Packet Output */
    {
      uint8_T outdata[8U];
      RTWin_ANYTYPEPTR outdp;
      outdp.p_uint8_T = outdata;

      {
        real_T pktout = rt_test_B.Sum1;
        *outdp.p_real_T++ = pktout;
      }

      RTBIO_DriverIO(0, STREAMOUTPUT, IOWRITE, 8U,
                     &rt_test_P.PacketOutput_PacketID, (double*) outdata, NULL);
    }

    /* Update for DiscreteZeroPole: '<S1>/Discrete Zero-Pole' */
    {
      real_T xnew[23];
      xnew[0] = (rt_test_P.DiscreteZeroPole_A[0])*
        rt_test_DW.DiscreteZeroPole_DSTATE[0];
      xnew[0] += (rt_test_P.DiscreteZeroPole_B[0])*rt_test_B.Sum;
      xnew[1] = (rt_test_P.DiscreteZeroPole_A[1])*
        rt_test_DW.DiscreteZeroPole_DSTATE[0]
        + (rt_test_P.DiscreteZeroPole_A[2])*rt_test_DW.DiscreteZeroPole_DSTATE[1]
        + (rt_test_P.DiscreteZeroPole_A[3])*rt_test_DW.DiscreteZeroPole_DSTATE[2];
      xnew[1] += (rt_test_P.DiscreteZeroPole_B[1])*rt_test_B.Sum;
      xnew[2] = (rt_test_P.DiscreteZeroPole_A[4])*
        rt_test_DW.DiscreteZeroPole_DSTATE[1];

      {
        static const int_T colAidxRow3[5] = { 0, 1, 2, 3, 4 };

        const int_T *pAidx = &colAidxRow3[0];
        const real_T *pA5 = &rt_test_P.DiscreteZeroPole_A[5];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew3 = &xnew[3];
        int_T numNonZero = 4;
        *pxnew3 = (*pA5++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew3 += (*pA5++) * xd[*pAidx++];
        }
      }

      xnew[3] += (rt_test_P.DiscreteZeroPole_B[2])*rt_test_B.Sum;
      xnew[4] = (rt_test_P.DiscreteZeroPole_A[10])*
        rt_test_DW.DiscreteZeroPole_DSTATE[3];

      {
        static const int_T colAidxRow5[7] = { 0, 1, 2, 3, 4, 5, 6 };

        const int_T *pAidx = &colAidxRow5[0];
        const real_T *pA11 = &rt_test_P.DiscreteZeroPole_A[11];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew5 = &xnew[5];
        int_T numNonZero = 6;
        *pxnew5 = (*pA11++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew5 += (*pA11++) * xd[*pAidx++];
        }
      }

      xnew[5] += (rt_test_P.DiscreteZeroPole_B[3])*rt_test_B.Sum;
      xnew[6] = (rt_test_P.DiscreteZeroPole_A[18])*
        rt_test_DW.DiscreteZeroPole_DSTATE[5];

      {
        static const int_T colAidxRow7[9] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

        const int_T *pAidx = &colAidxRow7[0];
        const real_T *pA19 = &rt_test_P.DiscreteZeroPole_A[19];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew7 = &xnew[7];
        int_T numNonZero = 8;
        *pxnew7 = (*pA19++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew7 += (*pA19++) * xd[*pAidx++];
        }
      }

      xnew[7] += (rt_test_P.DiscreteZeroPole_B[4])*rt_test_B.Sum;
      xnew[8] = (rt_test_P.DiscreteZeroPole_A[28])*
        rt_test_DW.DiscreteZeroPole_DSTATE[7];

      {
        static const int_T colAidxRow9[11] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        };

        const int_T *pAidx = &colAidxRow9[0];
        const real_T *pA29 = &rt_test_P.DiscreteZeroPole_A[29];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew9 = &xnew[9];
        int_T numNonZero = 10;
        *pxnew9 = (*pA29++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew9 += (*pA29++) * xd[*pAidx++];
        }
      }

      xnew[9] += (rt_test_P.DiscreteZeroPole_B[5])*rt_test_B.Sum;
      xnew[10] = (rt_test_P.DiscreteZeroPole_A[40])*
        rt_test_DW.DiscreteZeroPole_DSTATE[9];

      {
        static const int_T colAidxRow11[13] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12 };

        const int_T *pAidx = &colAidxRow11[0];
        const real_T *pA41 = &rt_test_P.DiscreteZeroPole_A[41];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew11 = &xnew[11];
        int_T numNonZero = 12;
        *pxnew11 = (*pA41++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew11 += (*pA41++) * xd[*pAidx++];
        }
      }

      xnew[11] += (rt_test_P.DiscreteZeroPole_B[6])*rt_test_B.Sum;
      xnew[12] = (rt_test_P.DiscreteZeroPole_A[54])*
        rt_test_DW.DiscreteZeroPole_DSTATE[11];

      {
        static const int_T colAidxRow13[15] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14 };

        const int_T *pAidx = &colAidxRow13[0];
        const real_T *pA55 = &rt_test_P.DiscreteZeroPole_A[55];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew13 = &xnew[13];
        int_T numNonZero = 14;
        *pxnew13 = (*pA55++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew13 += (*pA55++) * xd[*pAidx++];
        }
      }

      xnew[13] += (rt_test_P.DiscreteZeroPole_B[7])*rt_test_B.Sum;
      xnew[14] = (rt_test_P.DiscreteZeroPole_A[70])*
        rt_test_DW.DiscreteZeroPole_DSTATE[13];

      {
        static const int_T colAidxRow15[17] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16 };

        const int_T *pAidx = &colAidxRow15[0];
        const real_T *pA71 = &rt_test_P.DiscreteZeroPole_A[71];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew15 = &xnew[15];
        int_T numNonZero = 16;
        *pxnew15 = (*pA71++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew15 += (*pA71++) * xd[*pAidx++];
        }
      }

      xnew[15] += (rt_test_P.DiscreteZeroPole_B[8])*rt_test_B.Sum;
      xnew[16] = (rt_test_P.DiscreteZeroPole_A[88])*
        rt_test_DW.DiscreteZeroPole_DSTATE[15];

      {
        static const int_T colAidxRow17[19] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16, 17, 18 };

        const int_T *pAidx = &colAidxRow17[0];
        const real_T *pA89 = &rt_test_P.DiscreteZeroPole_A[89];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew17 = &xnew[17];
        int_T numNonZero = 18;
        *pxnew17 = (*pA89++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew17 += (*pA89++) * xd[*pAidx++];
        }
      }

      xnew[17] += (rt_test_P.DiscreteZeroPole_B[9])*rt_test_B.Sum;
      xnew[18] = (rt_test_P.DiscreteZeroPole_A[108])*
        rt_test_DW.DiscreteZeroPole_DSTATE[17];

      {
        static const int_T colAidxRow19[21] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

        const int_T *pAidx = &colAidxRow19[0];
        const real_T *pA109 = &rt_test_P.DiscreteZeroPole_A[109];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew19 = &xnew[19];
        int_T numNonZero = 20;
        *pxnew19 = (*pA109++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew19 += (*pA109++) * xd[*pAidx++];
        }
      }

      xnew[19] += (rt_test_P.DiscreteZeroPole_B[10])*rt_test_B.Sum;
      xnew[20] = (rt_test_P.DiscreteZeroPole_A[130])*
        rt_test_DW.DiscreteZeroPole_DSTATE[19];

      {
        static const int_T colAidxRow21[23] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
          11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22 };

        const int_T *pAidx = &colAidxRow21[0];
        const real_T *pA131 = &rt_test_P.DiscreteZeroPole_A[131];
        const real_T *xd = &rt_test_DW.DiscreteZeroPole_DSTATE[0];
        real_T *pxnew21 = &xnew[21];
        int_T numNonZero = 22;
        *pxnew21 = (*pA131++) * xd[*pAidx++];
        while (numNonZero--) {
          *pxnew21 += (*pA131++) * xd[*pAidx++];
        }
      }

      xnew[21] += (rt_test_P.DiscreteZeroPole_B[11])*rt_test_B.Sum;
      xnew[22] = (rt_test_P.DiscreteZeroPole_A[154])*
        rt_test_DW.DiscreteZeroPole_DSTATE[21];
      (void) memcpy(&rt_test_DW.DiscreteZeroPole_DSTATE[0], xnew,
                    sizeof(real_T)*23);
    }
  }

  if (rtmIsMajorTimeStep(rt_test_M)) {
    rt_ertODEUpdateContinuousStates(&rt_test_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++rt_test_M->Timing.clockTick0)) {
    ++rt_test_M->Timing.clockTickH0;
  }

  rt_test_M->Timing.t[0] = rtsiGetSolverStopTime(&rt_test_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.01s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++rt_test_M->Timing.clockTick1)) {
      ++rt_test_M->Timing.clockTickH1;
    }

    rt_test_M->Timing.t[1] = rt_test_M->Timing.clockTick1 *
      rt_test_M->Timing.stepSize1 + rt_test_M->Timing.clockTickH1 *
      rt_test_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void rt_test_derivatives(void)
{
  XDot_rt_test_T *_rtXdot;
  _rtXdot = ((XDot_rt_test_T *) rt_test_M->derivs);

  /* Derivatives for TransferFcn: '<Root>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += rt_test_P.TransferFcn1_A *
    rt_test_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += rt_test_B.DiscreteZeroPole;

  /* Derivatives for TransferFcn: '<Root>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += rt_test_P.TransferFcn_A *
    rt_test_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += rt_test_B.PacketInput;
}

/* Model initialize function */
void rt_test_initialize(void)
{
  /* Start for S-Function (sldrtpo): '<Root>/Packet Output' */

  /* S-Function Block: <Root>/Packet Output */
  /* no initial value should be set */

  /* InitializeConditions for TransferFcn: '<Root>/Transfer Fcn1' */
  rt_test_X.TransferFcn1_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<Root>/Transfer Fcn' */
  rt_test_X.TransferFcn_CSTATE = 0.0;
}

/* Model terminate function */
void rt_test_terminate(void)
{
  /* Terminate for S-Function (sldrtpo): '<Root>/Packet Output' */

  /* S-Function Block: <Root>/Packet Output */
  /* no initial value should be set */
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  rt_test_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  rt_test_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  rt_test_initialize();
}

void MdlTerminate(void)
{
  rt_test_terminate();
}

/* Registration function */
RT_MODEL_rt_test_T *rt_test(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)rt_test_M, 0,
                sizeof(RT_MODEL_rt_test_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&rt_test_M->solverInfo, &rt_test_M->Timing.simTimeStep);
    rtsiSetTPtr(&rt_test_M->solverInfo, &rtmGetTPtr(rt_test_M));
    rtsiSetStepSizePtr(&rt_test_M->solverInfo, &rt_test_M->Timing.stepSize0);
    rtsiSetdXPtr(&rt_test_M->solverInfo, &rt_test_M->derivs);
    rtsiSetContStatesPtr(&rt_test_M->solverInfo, (real_T **)
                         &rt_test_M->contStates);
    rtsiSetNumContStatesPtr(&rt_test_M->solverInfo,
      &rt_test_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&rt_test_M->solverInfo,
      &rt_test_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&rt_test_M->solverInfo,
      &rt_test_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&rt_test_M->solverInfo,
      &rt_test_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&rt_test_M->solverInfo, (&rtmGetErrorStatus(rt_test_M)));
    rtsiSetRTModelPtr(&rt_test_M->solverInfo, rt_test_M);
  }

  rtsiSetSimTimeStep(&rt_test_M->solverInfo, MAJOR_TIME_STEP);
  rt_test_M->intgData.y = rt_test_M->odeY;
  rt_test_M->intgData.f[0] = rt_test_M->odeF[0];
  rt_test_M->intgData.f[1] = rt_test_M->odeF[1];
  rt_test_M->intgData.f[2] = rt_test_M->odeF[2];
  rt_test_M->intgData.f[3] = rt_test_M->odeF[3];
  rt_test_M->intgData.f[4] = rt_test_M->odeF[4];
  rt_test_M->intgData.f[5] = rt_test_M->odeF[5];
  rt_test_M->contStates = ((real_T *) &rt_test_X);
  rtsiSetSolverData(&rt_test_M->solverInfo, (void *)&rt_test_M->intgData);
  rtsiSetSolverName(&rt_test_M->solverInfo,"ode5");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = rt_test_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    rt_test_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    rt_test_M->Timing.sampleTimes = (&rt_test_M->Timing.sampleTimesArray[0]);
    rt_test_M->Timing.offsetTimes = (&rt_test_M->Timing.offsetTimesArray[0]);

    /* task periods */
    rt_test_M->Timing.sampleTimes[0] = (0.0);
    rt_test_M->Timing.sampleTimes[1] = (0.01);

    /* task offsets */
    rt_test_M->Timing.offsetTimes[0] = (0.0);
    rt_test_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(rt_test_M, &rt_test_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = rt_test_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    rt_test_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(rt_test_M, 20.0);
  rt_test_M->Timing.stepSize0 = 0.01;
  rt_test_M->Timing.stepSize1 = 0.01;

  /* External mode info */
  rt_test_M->Sizes.checksums[0] = (1421901802U);
  rt_test_M->Sizes.checksums[1] = (132779391U);
  rt_test_M->Sizes.checksums[2] = (2264187680U);
  rt_test_M->Sizes.checksums[3] = (721584221U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    rt_test_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(rt_test_M->extModeInfo,
      &rt_test_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(rt_test_M->extModeInfo, rt_test_M->Sizes.checksums);
    rteiSetTPtr(rt_test_M->extModeInfo, rtmGetTPtr(rt_test_M));
  }

  rt_test_M->solverInfoPtr = (&rt_test_M->solverInfo);
  rt_test_M->Timing.stepSize = (0.01);
  rtsiSetFixedStepSize(&rt_test_M->solverInfo, 0.01);
  rtsiSetSolverMode(&rt_test_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  rt_test_M->blockIO = ((void *) &rt_test_B);
  (void) memset(((void *) &rt_test_B), 0,
                sizeof(B_rt_test_T));

  /* parameters */
  rt_test_M->defaultParam = ((real_T *)&rt_test_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &rt_test_X;
    rt_test_M->contStates = (x);
    (void) memset((void *)&rt_test_X, 0,
                  sizeof(X_rt_test_T));
  }

  /* states (dwork) */
  rt_test_M->dwork = ((void *) &rt_test_DW);
  (void) memset((void *)&rt_test_DW, 0,
                sizeof(DW_rt_test_T));

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    rt_test_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 14;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  rt_test_M->Sizes.numContStates = (2);/* Number of continuous states */
  rt_test_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  rt_test_M->Sizes.numY = (0);         /* Number of model outputs */
  rt_test_M->Sizes.numU = (0);         /* Number of model inputs */
  rt_test_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  rt_test_M->Sizes.numSampTimes = (2); /* Number of sample times */
  rt_test_M->Sizes.numBlocks = (10);   /* Number of blocks */
  rt_test_M->Sizes.numBlockIO = (6);   /* Number of block outputs */
  rt_test_M->Sizes.numBlockPrms = (202);/* Sum of parameter "widths" */
  return rt_test_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
