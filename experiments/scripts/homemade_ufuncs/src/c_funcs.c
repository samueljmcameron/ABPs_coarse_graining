#include <math.h>
#include "s_and_q_funcs/header.h"
#include "s_and_q_funcs/chbevl.h"



double i1_V_int(double r_in, double Pi, double r0)
{
  /*
    computing function

    1/{\epsilon}*int_{r_0}^{r_in} i1(\sqrt{Pi} u) V'(u) du,

    using the fact that

    s_1(x) = int_{0.5}^x i1(u)/u^12,

    s_2(x) = int_{0.5}^x i1(u)/u^6.

  */

  double sqrtPi = sqrt(Pi);
    
  double s1_tmp = 2*Pi*Pi*Pi*( s1(r_in*sqrtPi) - s1(r0*sqrtPi) );
    
  double s2_tmp = ( s2(r_in*sqrtPi) - s2(r0*sqrtPi) );
    
  return -24*Pi*Pi*sqrtPi*(s1_tmp - s2_tmp);

}

double k1_V_int(double r_in, double Pi, double r0)
{
  /*
    computing function

    1/{\epsilon}*int_{r_0}^{r_in} k1(\sqrt{Pi} u) V'(u) du,

    using the fact that

    q_1(x) = int_{0.5}^x k1(u)/u^12,

    q_2(x) = int_{0.5}^x k1(u)/u^6.

  */
  
  double sqrtPi = sqrt(Pi);
    
  double q1_tmp = 2*Pi*Pi*Pi*( q1(r_in*sqrtPi) - q1(r0*sqrtPi) );
    
  double q2_tmp = ( q2(r_in*sqrtPi) - q2(r0*sqrtPi) );
    
  return -24*Pi*Pi*sqrtPi*(q1_tmp - q2_tmp);
}
    

