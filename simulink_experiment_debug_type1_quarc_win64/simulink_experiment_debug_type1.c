/*
 * simulink_experiment_debug_type1.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.0
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Mon Apr  7 14:57:34 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "simulink_experiment_debug_type1.h"
#include "simulink_experiment_debug_type1_types.h"
#include <math.h>
#include "rt_nonfinite.h"
#include <string.h>
#include <emmintrin.h>
#include "rtwtypes.h"
#include "simulink_experiment_debug_type1_private.h"
#include "simulink_experiment_debug_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_debug_t_T simulink_experiment_debug_typ_B;

/* Block states (default storage) */
DW_simulink_experiment_debug__T simulink_experiment_debug_ty_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_debug_ty_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_debug_ty_M =
  &simulink_experiment_debug_ty_M_;
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(simulink_experiment_debug_ty_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(simulink_experiment_debug_ty_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *         This function updates active task flag for each subrate
 *         and rate transition flags for tasks that exchange data.
 *         The function assumes rate-monotonic multitasking scheduler.
 *         The function must be called at model base rate so that
 *         the generated code self-manages all its subrates and rate
 *         transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[1] == 0) {
    simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2 =
      (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits[5] =
      simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2])++;
  if ((simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  print_discrepancies_simulink__T *obj_0;
  studentControllerInterface_si_T *obj;
  real_T a[8];
  real_T dx_hat[4];
  real_T dx_hat_0[4];
  real_T x_hat_prev[4];
  real_T y_0[4];
  real_T y[2];
  real_T K;
  real_T amplitude;
  real_T c;
  real_T c_0;
  real_T c_1;
  real_T c_2;
  real_T omega;
  real_T t_prev;
  real_T u0;
  static const real_T tmp[4] = { 40.0, 41.774, 9.123, 1.731 };

  __m128d tmp_0;
  __m128d tmp_1;
  real_T y_idx_0;
  int32_T i;
  static const int8_T tmp_2[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

  static const real_T tmp_3[8] = { 25.132, 151.9416, -1.2948, 32.6879, -0.5331,
    -6.4568, -10.132, 619.3871 };

  {                                    /* Sample time: [0.0s, 0.0s] */
    rate_monotonic_scheduler();
  }

  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    } else {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase =
        simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_debug_ty_DW.HILInitialize_Card,
       &simulink_experiment_debug_typ_P.HILReadAnalog_channels, 1,
       &simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }

    simulink_experiment_debug_typ_B.HILReadAnalog =
      simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_debug_typ_B.BB01SensorGainmV =
    simulink_experiment_debug_typ_P.BB01SensorGainmV_Gain *
    simulink_experiment_debug_typ_B.HILReadAnalog;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_debug_typ_B.EncoderCalibrationradcount =
    simulink_experiment_debug_typ_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_debug_typ_B.Bias =
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount +
    simulink_experiment_debug_typ_P.Bias_Bias;

  /* Clock: '<Root>/Clock' */
  simulink_experiment_debug_typ_B.Clock =
    simulink_experiment_debug_ty_M->Timing.t[0];

  /* MATLABSystem: '<Root>/MATLAB System' */
  u0 = simulink_experiment_debug_typ_B.Clock;
  amplitude = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  omega = simulink_experiment_debug_typ_B.Bias;
  obj = &simulink_experiment_debug_ty_DW.obj;

  /*  function setupImpl(obj) */
  /*     disp("You can use this function for initializaition."); */
  /*  end */
  /*  This is the main function called every iteration. You have to implement */
  /*  the controller in this function, bu you are not allowed to */
  /*  change the signature of this function.  */
  /*  Input arguments: */
  /*    t: current time */
  /*    p_ball: position of the ball provided by the ball position sensor (m) */
  /*  */
  /*    theta: servo motor angle provided by the encoder of the motor (rad) */
  /*  Output: */
  /*    V_servo: voltage to the servo input.         */
  /*             %% Sample Controller: Simple Proportional Controller */
  t_prev = obj->t_prev;
  x_hat_prev[0] = obj->x_hat_prev[0];
  x_hat_prev[1] = obj->x_hat_prev[1];
  x_hat_prev[2] = obj->x_hat_prev[2];
  x_hat_prev[3] = obj->x_hat_prev[3];
  y_idx_0 = amplitude;

  /*  System parameters */
  /*  Extract reference trajectory at the current timestep. */
  /*  [p_ref, v_ref, a_ref] = get_ref_traj(t) */
  /*  Reference trajectory at time t: */
  /*    Inputs */
  /*        t: query time */
  /*    Outputs */
  /*        p_ref: reference position of the ball */
  /*        v_ref: reference velocity of the ball */
  /*        a_ref: reference acceleration of the ball */
  /*  m */
  /*  sec */
  /*  %% Sine wave. */
  /*      p_ref = amplitude * sin(omega * t); */
  /*      v_ref = amplitude * omega * cos(omega * t); */
  /*      a_ref = - amplitude * omega^2 * sin(omega * t); */
  /*  Square wave. */
  K = 0.25132741228718347 * u0;
  K = sin(K);
  if (rtIsNaN(K)) {
    K = (rtNaN);
  } else if (K < 0.0) {
    K = -1.0;
  } else {
    K = (K > 0.0);
  }

  amplitude = 0.04 * K;

  /*  state_estimate -- luenberger observer */
  t_prev = u0 - t_prev;

  /*  x_hat_prev -- 4x1 vector */
  /*  y -- 2x1 vector */
  /*  x_op = [0; 0; 0; 0]; % [p_ball; v_ball; theta; dtheta] */
  /*  Input at equilibrium */
  /*  Compute A and B matrices at equilibrium */
  /*  Compute the Jacobians */
  /*      [~, A] = jaccsd(@(x) ball_and_beam_dynamics(0, x, u_op), x_op); */
  /*      [~, B] = jaccsd(@(u) ball_and_beam_dynamics(0, x_op, u), u_op); */
  /*  Output equations (assuming we measure p_ball and theta) */
  /* h = @(x) [x(1); x(3)]; */
  /* [~, C] = jaccsd(h, x_op); */
  /* D = zeros(2, 1); % No direct feedthrough */
  /* G = eye(4); */
  /*  Covariance matrices */
  /* Q = [10, 0, 0, 0; 0, 100, 0, 0; 0, 0, 10, 0; 0, 0, 0, 100];%eye(4) * 10; % Process noise covariance */
  /* R = eye(2) * 0.1; % Measurement noise covariance */
  /* P = eye(4); */
  /* sys = ss(A,[B G],C,[D zeros(2, 4)]); */
  /* [kalmf,Ll,P] = kalman(sys,Q,R); */
  /* Luenberger gain */
  /* LO = Ll; */
  /* selects position and theta */
  /* observer_poles = [-10,-12,-15,-18]; */
  /* LO = place(A',C',observer_poles); */
  /*  x_hat dynamics */
  dx_hat[0] = x_hat_prev[1];
  K = x_hat_prev[3];
  c = K * K;
  K = x_hat_prev[2];
  K = cos(K);
  c_0 = K * K;
  K = x_hat_prev[3];
  c_1 = K * K;
  K = x_hat_prev[2];
  K = cos(K);
  c_2 = K * K;
  K = x_hat_prev[2];
  K = sin(K);
  dx_hat[1] = (0.41828772872251135 * K - 0.00054151418499244583 * c * c_0) +
    0.0025453075675320605 * x_hat_prev[0] * c_1 * c_2;
  dx_hat[2] = x_hat_prev[3];
  dx_hat[3] = -x_hat_prev[3] / 0.025;

  /*  Luenberger part */
  for (i = 0; i < 8; i++) {
    a[i] = tmp_2[i];
  }

  for (i = 0; i <= 0; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&a[i]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_hat_prev[0]));
    tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
    tmp_1 = _mm_loadu_pd(&a[i + 2]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[1]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_1 = _mm_loadu_pd(&a[i + 4]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[2]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_1 = _mm_loadu_pd(&a[i + 6]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[3]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&y[i], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  memcpy(&a[0], &tmp_3[0], sizeof(real_T) << 3U);
  K = y_idx_0;
  K -= y[0];
  y_idx_0 = K;
  K = omega;
  K -= y[1];
  omega = K;

  /*  x_hat update step */
  /* k_p = 3; */
  /* theta_d = - k_p * (p_ball - p_ball_ref); */
  /* theta_saturation = 56 * pi / 180;     */
  /*  Compute A and B matrices at equilibrium */
  /*  Solve LQR */
  /*  sine -- 0.95, square -- 4 */
  /*  Q = diag([1200, 10, 10, 10]); */
  /*  Klqr = lqr(A, B, Q, R); */
  /* Klqr = [10, 25.1525, 13.0233, 2.6315]; */
  /* Klqr = [10, 45, 11, 2.3]; % sine cost: 0.97, square -- 4.4 */
  /* Klqr = [24.5, 33.5, 9.6, 2.7]; % sine cost: 0.85, square --3.5 */
  /* Klqr = [67.08, 66.6, 13.534, 2.7]; % sine cost: , square --  */
  /*  sine cost: 0.8081, square --  */
  for (i = 0; i <= 2; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_loadu_pd(&a[i]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y_idx_0));
    tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
    tmp_1 = _mm_loadu_pd(&a[i + 4]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(omega));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_1 = _mm_loadu_pd(&dx_hat[i]);
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(t_prev));
    tmp_1 = _mm_loadu_pd(&x_hat_prev[i]);
    tmp_1 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System' */
    _mm_storeu_pd(&x_hat_prev[i], tmp_1);
    _mm_storeu_pd(&dx_hat[i], tmp_0);
    tmp_0 = _mm_loadu_pd(&tmp[i]);
    _mm_storeu_pd(&y_0[i], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System' */
  dx_hat[0] = amplitude;
  omega = dx_hat[0];
  omega = x_hat_prev[0] - omega;
  amplitude = y_0[0] * omega;
  omega = x_hat_prev[1];
  amplitude += y_0[1] * omega;
  omega = x_hat_prev[2];
  amplitude += y_0[2] * omega;
  omega = x_hat_prev[3];
  amplitude += y_0[3] * omega;
  omega = 0.0 - amplitude;

  /*             %% nonlinear input control ( works but at a higher energy cost) */
  /*  V_servo = 1*sign(u_eq - Klqr  * (x_hat - x_ref)); */
  /*             %% saturate V_servo */
  /*  lb = -1 perform better for square */
  /*  ub = 1 perform better for square */
  if (!(omega >= -1.0)) {
    omega = -1.0;
  }

  if (!(omega <= 3.0)) {
    omega = 3.0;
  }

  /*  % Decide desired servo angle based on simple proportional feedback. */
  /*  k_p = 3; */
  /*  theta_d = - k_p * (p_ball - p_ball_ref); */
  /*  % Make sure that the desired servo angle does not exceed the physical */
  /*  % limit. This part of code is not necessary but highly recommended */
  /*  % because it addresses the actual physical limit of the servo motor. */
  /*  theta_saturation = 56 * pi / 180;     */
  /*  theta_d = min(theta_d, theta_saturation); */
  /*  theta_d = max(theta_d, -theta_saturation); */
  /*  % Simple position control to control servo angle to the desired */
  /*  % position. */
  /*  k_servo = 10; */
  /*  V_servo = k_servo * (theta_d - theta); */
  /*  Update class properties if necessary. */
  obj->t_prev = u0;
  obj->theta_d = 0.0;

  /* obj.theta_d = x_hat(3); */
  obj->x_hat_prev[0] = x_hat_prev[0];
  obj->x_hat_prev[1] = x_hat_prev[1];
  obj->x_hat_prev[2] = x_hat_prev[2];
  obj->x_hat_prev[3] = x_hat_prev[3];

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem = omega;

  /* Saturate: '<Root>/+//-10V' */
  u0 = simulink_experiment_debug_typ_B.MATLABSystem;
  amplitude = simulink_experiment_debug_typ_P.u0V_LowerSat;
  omega = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (u0 > omega) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = omega;
  } else if (u0 < amplitude) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = amplitude;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u0;
  }

  /* End of Saturate: '<Root>/+//-10V' */

  /* Gain: '<S1>/Motor  Gain (V//V)' */
  simulink_experiment_debug_typ_B.MotorGainVV =
    simulink_experiment_debug_typ_P.MotorGainVV_Gain *
    simulink_experiment_debug_typ_B.u0V;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_debug_ty_DW.HILInitialize_Card,
      &simulink_experiment_debug_typ_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_debug_typ_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }
  }

  /* MATLAB Function: '<Root>/MATLAB Function' */
  /* MATLAB Function 'MATLAB Function': '<S2>:1' */
  u0 = sin(0.25132741228718347 * simulink_experiment_debug_typ_B.Clock);
  if (rtIsNaN(u0)) {
    K = (rtNaN);
  } else if (u0 < 0.0) {
    K = -1.0;
  } else {
    K = (u0 > 0.0);
  }

  simulink_experiment_debug_typ_B.p_ref = 0.04 * K;
  simulink_experiment_debug_typ_B.v_ref = 0.0;
  simulink_experiment_debug_typ_B.a_ref = 0.0;

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* Gain: '<Root>/m to cm' */
  /* '<S2>:1:3' */
  simulink_experiment_debug_typ_B.mtocm[0] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.p_ref;
  simulink_experiment_debug_typ_B.mtocm[1] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.BB01SensorGainmV;

  /* MATLABSystem: '<Root>/MATLAB System1' */
  u0 = simulink_experiment_debug_typ_B.Clock;
  amplitude = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  omega = simulink_experiment_debug_typ_B.Bias;
  obj_0 = &simulink_experiment_debug_ty_DW.obj_i;

  /*  function setupImpl(obj) */
  /*     disp("You can use this function for initializaition."); */
  /*  end */
  /*  This is the main function called every iteration. You have to implement */
  /*  the controller in this function, bu you are not allowed to */
  /*  change the signature of this function.  */
  /*  Input arguments: */
  /*    t: current time */
  /*    p_ball: position of the ball provided by the ball position sensor (m) */
  /*  */
  /*    theta: servo motor angle provided by the encoder of the motor (rad) */
  /*  Output: */
  /*    V_servo: voltage to the servo input.         */
  /*             %% Sample Controller: Simple Proportional Controller */
  t_prev = obj_0->t_prev;
  x_hat_prev[0] = obj_0->x_hat_prev[0];
  x_hat_prev[1] = obj_0->x_hat_prev[1];
  x_hat_prev[2] = obj_0->x_hat_prev[2];
  x_hat_prev[3] = obj_0->x_hat_prev[3];
  y_idx_0 = amplitude;

  /*  System parameters */
  /*  Extract reference trajectory at the current timestep. */
  /*  [p_ref, v_ref, a_ref] = get_ref_traj(t) */
  /*  Reference trajectory at time t: */
  /*    Inputs */
  /*        t: query time */
  /*    Outputs */
  /*        p_ref: reference position of the ball */
  /*        v_ref: reference velocity of the ball */
  /*        a_ref: reference acceleration of the ball */
  /*  m */
  /*  sec */
  /*  %% Sine wave. */
  /*      p_ref = amplitude * sin(omega * t); */
  /*      v_ref = amplitude * omega * cos(omega * t); */
  /*      a_ref = - amplitude * omega^2 * sin(omega * t); */
  /*  Square wave. */
  /*  state_estimate -- luenberger observer */
  t_prev = u0 - t_prev;

  /*  x_hat_prev -- 4x1 vector */
  /*  y -- 2x1 vector */
  /*  x_op = [0; 0; 0; 0]; % [p_ball; v_ball; theta; dtheta] */
  /*  Input at equilibrium */
  /*  Compute A and B matrices at equilibrium */
  /*  Compute the Jacobians */
  /*      [~, A] = jaccsd(@(x) ball_and_beam_dynamics(0, x, u_op), x_op); */
  /*      [~, B] = jaccsd(@(u) ball_and_beam_dynamics(0, x_op, u), u_op); */
  /*  Output equations (assuming we measure p_ball and theta) */
  /* h = @(x) [x(1); x(3)]; */
  /* [~, C] = jaccsd(h, x_op); */
  /* D = zeros(2, 1); % No direct feedthrough */
  /* G = eye(4); */
  /*  Covariance matrices */
  /* Q = [10, 0, 0, 0; 0, 100, 0, 0; 0, 0, 10, 0; 0, 0, 0, 100];%eye(4) * 10; % Process noise covariance */
  /* R = eye(2) * 0.1; % Measurement noise covariance */
  /* P = eye(4); */
  /* sys = ss(A,[B G],C,[D zeros(2, 4)]); */
  /* [kalmf,Ll,P] = kalman(sys,Q,R); */
  /* Luenberger gain */
  /* LO = Ll; */
  /* selects position and theta */
  /* observer_poles = [-10,-12,-15,-18]; */
  /* LO = place(A',C',observer_poles); */
  /*  x_hat dynamics */
  dx_hat_0[0] = x_hat_prev[1];
  K = x_hat_prev[3];
  c = K * K;
  K = x_hat_prev[2];
  K = cos(K);
  c_0 = K * K;
  K = x_hat_prev[3];
  c_1 = K * K;
  K = x_hat_prev[2];
  K = cos(K);
  c_2 = K * K;
  K = x_hat_prev[2];
  K = sin(K);
  dx_hat_0[1] = (0.41828772872251135 * K - 0.00054151418499244583 * c * c_0) +
    0.0025453075675320605 * x_hat_prev[0] * c_1 * c_2;
  dx_hat_0[2] = x_hat_prev[3];
  dx_hat_0[3] = -x_hat_prev[3] / 0.025;

  /*  Luenberger part */
  for (i = 0; i < 8; i++) {
    a[i] = tmp_2[i];
  }

  for (i = 0; i <= 0; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_0 = _mm_loadu_pd(&a[i]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_hat_prev[0]));
    tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
    tmp_1 = _mm_loadu_pd(&a[i + 2]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[1]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_1 = _mm_loadu_pd(&a[i + 4]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[2]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_1 = _mm_loadu_pd(&a[i + 6]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(x_hat_prev[3]));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    _mm_storeu_pd(&y[i], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System1' */
  memcpy(&a[0], &tmp_3[0], sizeof(real_T) << 3U);
  K = y_idx_0;
  K -= y[0];
  y_idx_0 = K;
  K = omega;
  K -= y[1];
  omega = K;

  /*  x_hat update step */
  for (i = 0; i <= 2; i += 2) {
    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_0 = _mm_loadu_pd(&a[i]);
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y_idx_0));
    tmp_0 = _mm_add_pd(tmp_0, _mm_set1_pd(0.0));
    tmp_1 = _mm_loadu_pd(&a[i + 4]);
    tmp_1 = _mm_mul_pd(tmp_1, _mm_set1_pd(omega));
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_1 = _mm_loadu_pd(&dx_hat_0[i]);
    tmp_0 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(t_prev));
    tmp_1 = _mm_loadu_pd(&x_hat_prev[i]);
    tmp_1 = _mm_add_pd(tmp_1, tmp_0);

    /* MATLABSystem: '<Root>/MATLAB System1' */
    _mm_storeu_pd(&x_hat_prev[i], tmp_1);
    _mm_storeu_pd(&dx_hat_0[i], tmp_0);
  }

  /* MATLABSystem: '<Root>/MATLAB System1' */
  /* k_p = 3; */
  /* theta_d = - k_p * (p_ball - p_ball_ref); */
  /* theta_saturation = 56 * pi / 180;     */
  /*  Compute A and B matrices at equilibrium */
  /*  Solve LQR */
  /*  sine -- 0.95, square -- 4 */
  /*  Q = diag([1200, 10, 10, 10]); */
  /*  Klqr = lqr(A, B, Q, R); */
  /* Klqr = [10, 25.1525, 13.0233, 2.6315]; */
  /* Klqr = [10, 45, 11, 2.3]; % sine cost: 0.97, square -- 4.4 */
  /* Klqr = [24.5, 33.5, 9.6, 2.7]; % sine cost: 0.85, square --3.5 */
  /* Klqr = [67.08, 66.6, 13.534, 2.7]; % sine cost: , square --  */
  /*  sine cost: 0.8081, square --  */
  /*             %% nonlinear input control ( works but at a higher energy cost) */
  /*  V_servo = 1*sign(u_eq - Klqr  * (x_hat - x_ref)); */
  /*             %% saturate V_servo */
  /*  lb = -1 perform better for square */
  /*  ub = 1 perform better for square */
  amplitude = x_hat_prev[0] - amplitude;

  /*  % Decide desired servo angle based on simple proportional feedback. */
  /*  k_p = 3; */
  /*  theta_d = - k_p * (p_ball - p_ball_ref); */
  /*  % Make sure that the desired servo angle does not exceed the physical */
  /*  % limit. This part of code is not necessary but highly recommended */
  /*  % because it addresses the actual physical limit of the servo motor. */
  /*  theta_saturation = 56 * pi / 180;     */
  /*  theta_d = min(theta_d, theta_saturation); */
  /*  theta_d = max(theta_d, -theta_saturation); */
  /*  % Simple position control to control servo angle to the desired */
  /*  % position. */
  /*  k_servo = 10; */
  /*  V_servo = k_servo * (theta_d - theta); */
  /*  Update class properties if necessary. */
  obj_0->t_prev = u0;
  obj_0->theta_d = 0.0;

  /* obj.theta_d = x_hat(3); */
  obj_0->x_hat_prev[0] = x_hat_prev[0];
  obj_0->x_hat_prev[1] = x_hat_prev[1];
  obj_0->x_hat_prev[2] = x_hat_prev[2];
  obj_0->x_hat_prev[3] = x_hat_prev[3];

  /* MATLABSystem: '<Root>/MATLAB System1' */
  simulink_experiment_debug_typ_B.MATLABSystem1 = amplitude;

  /* Gain: '<S3>/Gain' */
  simulink_experiment_debug_typ_B.Gain =
    simulink_experiment_debug_typ_P.Gain_Gain *
    simulink_experiment_debug_typ_B.Bias;

  /* RateTransition: '<Root>/Rate Transition' */
  if (simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2) {
    simulink_experiment_debug_ty_DW.RateTransition_Buffer =
      simulink_experiment_debug_typ_B.Clock;

    /* RateTransition: '<Root>/Rate Transition1' */
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer =
      simulink_experiment_debug_typ_B.p_ref;

    /* RateTransition: '<Root>/Rate Transition2' */
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer =
      simulink_experiment_debug_typ_B.MATLABSystem;

    /* RateTransition: '<Root>/Rate Transition3' */
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer =
      simulink_experiment_debug_typ_B.BB01SensorGainmV;

    /* RateTransition: '<Root>/Rate Transition4' */
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer =
      simulink_experiment_debug_typ_B.Bias;
  }

  /* End of RateTransition: '<Root>/Rate Transition' */
}

/* Model update function for TID0 */
void simulink_experiment_debug_type1_update0(void) /* Sample time: [0.0s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick0)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH0;
  }

  simulink_experiment_debug_ty_M->Timing.t[0] =
    simulink_experiment_debug_ty_M->Timing.clockTick0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 +
    simulink_experiment_debug_ty_M->Timing.clockTickH0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 * 4294967296.0;

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick1)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH1;
  }

  simulink_experiment_debug_ty_M->Timing.t[1] =
    simulink_experiment_debug_ty_M->Timing.clockTick1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 +
    simulink_experiment_debug_ty_M->Timing.clockTickH1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 * 4294967296.0;
}

/* Model output function for TID2 */
void simulink_experiment_debug_type1_output2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* RateTransition: '<Root>/Rate Transition2' */
  simulink_experiment_debug_typ_B.RateTransition2 =
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer;

  /* RateTransition: '<Root>/Rate Transition1' */
  simulink_experiment_debug_typ_B.RateTransition1 =
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer;

  /* RateTransition: '<Root>/Rate Transition3' */
  simulink_experiment_debug_typ_B.RateTransition3 =
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer;

  /* RateTransition: '<Root>/Rate Transition4' */
  simulink_experiment_debug_typ_B.RateTransition4 =
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer;

  /* RateTransition: '<Root>/Rate Transition' */
  simulink_experiment_debug_typ_B.RateTransition =
    simulink_experiment_debug_ty_DW.RateTransition_Buffer;
}

/* Model update function for TID2 */
void simulink_experiment_debug_type1_update2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick2)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH2;
  }

  simulink_experiment_debug_ty_M->Timing.t[2] =
    simulink_experiment_debug_ty_M->Timing.clockTick2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 +
    simulink_experiment_debug_ty_M->Timing.clockTickH2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 * 4294967296.0;
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_output(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_output0();
    break;

   case 2 :
    simulink_experiment_debug_type1_output2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_update(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_update0();
    break;

   case 2 :
    simulink_experiment_debug_type1_update2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Model initialize function */
void simulink_experiment_debug_type1_initialize(void)
{
  {
    print_discrepancies_simulink__T *b_obj_0;
    studentControllerInterface_si_T *b_obj;

    /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
    {
      t_int result;
      t_boolean is_switching;
      result = hil_open("q2_usb", "0",
                        &simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      is_switching = false;
      result = hil_set_card_specific_options
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      result = hil_watchdog_clear
        (simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        result = hil_set_analog_input_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        result = hil_set_analog_output_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        result = hil_write_analog
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_AOReset) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        result = hil_watchdog_set_analog_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      result = hil_set_digital_directions
        (simulink_experiment_debug_ty_DW.HILInitialize_Card, NULL, 0U,
         simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_DOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_DOEnter && is_switching))
      {
        {
          int_T i1;
          boolean_T *dw_DOBits =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOBits[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOInitial;
          }
        }

        result = hil_write_digital
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U,
           (t_boolean *) &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_DOReset) {
        {
          int_T i1;
          int32_T *dw_DOStates =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOStates[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOWatchdog;
          }
        }

        result = hil_watchdog_set_digital_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U, (const
            t_digital_state *)
           &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        result = hil_set_encoder_quadrature_mode
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           (t_encoder_quadrature_mode *)
           &simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        result = hil_set_encoder_counts
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }
    }

    /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_create_encoder_reader
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         simulink_experiment_debug_typ_P.HILReadEncoderTimebase_SamplesI,
         &simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Channels, 1,
         &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task);
      if (result >= 0) {
        result = hil_task_set_buffer_overflow_mode
          (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task,
           (t_buffer_overflow_mode)
           (simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Overflow - 1));
      }

      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
      }
    }

    /* Start for MATLABSystem: '<Root>/MATLAB System' */
    b_obj = &simulink_experiment_debug_ty_DW.obj;
    b_obj->t_prev = -1.0;
    b_obj->x_hat_prev[0] = 0.0;
    b_obj->x_hat_prev[1] = 0.0;
    b_obj->x_hat_prev[2] = -1.0;
    b_obj->x_hat_prev[3] = 0.0;
    b_obj->theta_d = 0.0;
    simulink_experiment_debug_ty_DW.objisempty_g = true;

    /* Start for MATLABSystem: '<Root>/MATLAB System1' */
    b_obj_0 = &simulink_experiment_debug_ty_DW.obj_i;
    b_obj_0->t_prev = -1.0;
    b_obj_0->x_hat_prev[0] = 0.0;
    b_obj_0->x_hat_prev[1] = 0.0;
    b_obj_0->x_hat_prev[2] = -1.0;
    b_obj_0->x_hat_prev[3] = 0.0;
    b_obj_0->theta_d = 0.0;
    simulink_experiment_debug_ty_DW.objisempty = true;
  }
}

/* Model terminate function */
void simulink_experiment_debug_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_debug_typ_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_debug_typ_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_debug_typ_P.HILInitialize_DOFinal;
        }
      }

      num_final_digital_outputs = 8U;
    } else {
      num_final_digital_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_digital_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(simulink_experiment_debug_ty_DW.HILInitialize_Card
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_close(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    simulink_experiment_debug_ty_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_update(tid);
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
  simulink_experiment_debug_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_debug_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_debug_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_debug_ty_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          &simulink_experiment_debug_ty_M->Timing.simTimeStep);
    rtsiSetTPtr(&simulink_experiment_debug_ty_M->solverInfo, &rtmGetTPtr
                (simulink_experiment_debug_ty_M));
    rtsiSetStepSizePtr(&simulink_experiment_debug_ty_M->solverInfo,
                       &simulink_experiment_debug_ty_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          (&rtmGetErrorStatus(simulink_experiment_debug_ty_M)));
    rtsiSetRTModelPtr(&simulink_experiment_debug_ty_M->solverInfo,
                      simulink_experiment_debug_ty_M);
  }

  rtsiSetSimTimeStep(&simulink_experiment_debug_ty_M->solverInfo,
                     MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange
    (&simulink_experiment_debug_ty_M->solverInfo, false);
  rtsiSetSolverName(&simulink_experiment_debug_ty_M->solverInfo,
                    "FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_debug_ty_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_debug_ty_M->Timing.sampleTimes =
      (&simulink_experiment_debug_ty_M->Timing.sampleTimesArray[0]);
    simulink_experiment_debug_ty_M->Timing.offsetTimes =
      (&simulink_experiment_debug_ty_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_debug_ty_M->Timing.sampleTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[1] = (0.002);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    simulink_experiment_debug_ty_M->Timing.offsetTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[1] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_debug_ty_M,
             &simulink_experiment_debug_ty_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_debug_ty_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      simulink_experiment_debug_ty_M->Timing.perTaskSampleHitsArray;
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits =
      (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    simulink_experiment_debug_ty_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_debug_ty_M, 100.0);
  simulink_experiment_debug_ty_M->Timing.stepSize0 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize1 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (3133384680U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (1489600088U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (2937111737U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (1314312827U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[4];
    simulink_experiment_debug_ty_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    systemRan[3] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_debug_ty_M->extModeInfo,
      &simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_debug_ty_M->extModeInfo,
                        simulink_experiment_debug_ty_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_debug_ty_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_debug_ty_M));
  }

  simulink_experiment_debug_ty_M->solverInfoPtr =
    (&simulink_experiment_debug_ty_M->solverInfo);
  simulink_experiment_debug_ty_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_debug_ty_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_debug_ty_M->solverInfo,
                    SOLVER_MODE_MULTITASKING);

  /* block I/O */
  simulink_experiment_debug_ty_M->blockIO = ((void *)
    &simulink_experiment_debug_typ_B);

  {
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_debug_typ_B.HILReadAnalog = 0.0;
    simulink_experiment_debug_typ_B.BB01SensorGainmV = 0.0;
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_debug_typ_B.Bias = 0.0;
    simulink_experiment_debug_typ_B.Clock = 0.0;
    simulink_experiment_debug_typ_B.u0V = 0.0;
    simulink_experiment_debug_typ_B.MotorGainVV = 0.0;
    simulink_experiment_debug_typ_B.mtocm[0] = 0.0;
    simulink_experiment_debug_typ_B.mtocm[1] = 0.0;
    simulink_experiment_debug_typ_B.Gain = 0.0;
    simulink_experiment_debug_typ_B.RateTransition2 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition1 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition3 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition4 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem1 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem = 0.0;
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* parameters */
  simulink_experiment_debug_ty_M->defaultParam = ((real_T *)
    &simulink_experiment_debug_typ_P);

  /* states (dwork) */
  simulink_experiment_debug_ty_M->dwork = ((void *)
    &simulink_experiment_debug_ty_DW);
  (void) memset((void *)&simulink_experiment_debug_ty_DW, 0,
                sizeof(DW_simulink_experiment_debug__T));
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition1_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition2_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition3_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition4_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 23;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_debug_ty_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_debug_ty_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_debug_ty_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_debug_ty_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_debug_ty_M->Sizes.numSampTimes = (3);/* Number of sample times */
  simulink_experiment_debug_ty_M->Sizes.numBlocks = (32);/* Number of blocks */
  simulink_experiment_debug_ty_M->Sizes.numBlockIO = (20);/* Number of block outputs */
  simulink_experiment_debug_ty_M->Sizes.numBlockPrms = (88);/* Sum of parameter "widths" */
  return simulink_experiment_debug_ty_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
