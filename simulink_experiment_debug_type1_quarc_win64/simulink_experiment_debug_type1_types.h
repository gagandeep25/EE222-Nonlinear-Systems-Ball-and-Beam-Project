/*
 * simulink_experiment_debug_type1_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.2
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Mon Apr 21 15:07:09 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_simulink_experiment_debug_type1_types_h_
#define RTW_HEADER_simulink_experiment_debug_type1_types_h_
#include "rtwtypes.h"
#ifndef struct_tag_nMbvhz1rQO3iurRnLnfOEH
#define struct_tag_nMbvhz1rQO3iurRnLnfOEH

struct tag_nMbvhz1rQO3iurRnLnfOEH
{
  real_T t_prev;
  real_T x_hat_prev[4];
  real_T e_prev;
  real_T e_integral;
  real_T u_prev;
};

#endif                                 /* struct_tag_nMbvhz1rQO3iurRnLnfOEH */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_nMbvhz1rQO3iurRnLnfOEH studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

#ifndef struct_tag_7s6bAsLxjGRG4VeLPjkpZH
#define struct_tag_7s6bAsLxjGRG4VeLPjkpZH

struct tag_7s6bAsLxjGRG4VeLPjkpZH
{
  real_T t_prev;
  real_T x_hat_prev[4];
  real_T theta_d;
};

#endif                                 /* struct_tag_7s6bAsLxjGRG4VeLPjkpZH */

#ifndef typedef_print_discrepancies_simulink__T
#define typedef_print_discrepancies_simulink__T

typedef struct tag_7s6bAsLxjGRG4VeLPjkpZH print_discrepancies_simulink__T;

#endif                             /* typedef_print_discrepancies_simulink__T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_debug_t_T_ P_simulink_experiment_debug_t_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_d_T RT_MODEL_simulink_experiment__T;

#endif                 /* RTW_HEADER_simulink_experiment_debug_type1_types_h_ */
