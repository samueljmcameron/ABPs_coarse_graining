#include <math.h>
#include "chbevl.h"

/* Starting with
 *	 s_1(x) = int_{0.5}^x I_1(u)/u^12 du
 * function. */

/* coeffs for s_1(x) with 0.5 < x < 1.0 */
static double s_1_A[] = 
  {
    -1.702735175079795e-14,
    -1.433107105208720e-14,
    5.463511587588954e-15,
    3.795713743315332e-14,
    3.794568825821188e-14,
    1.290981210821940e-14,
    -3.794430047943109e-14,
    2.514377595019823e-14,
    2.032402024454427e-15,
    -5.018832571757059e-14,
    2.333966353518235e-14,
    -1.280427153194097e-13,
    7.421153969122685e-13,
    -3.126746084269882e-12,
    1.357821910463741e-11,
    -5.857593923797033e-11,
    2.494124684560006e-10,
    -1.050430034621463e-09,
    4.369853330576934e-09,
    -1.794073322819445e-08,
    7.261357310023309e-08,
    -2.893813460733141e-07,
    1.133966887122140e-06,
    -4.362379903161129e-06,
    1.644580808164708e-05,
    -6.063019851911683e-05,
    2.180567356646554e-04,
    -7.628835018099328e-04,
    2.587546472553325e-03,
    -8.474300622463541e-03,
    2.666738107779117e-02,
    -8.015149801972403e-02,
    2.283779522356360e-01,
    -6.110863892047467e-01,
    1.516914187104991e+00,
    -3.437567942059640e+00,
    6.959376368387966e+00,
    -1.221773986490679e+01,
    1.785070871181493e+01,
    8.588143208035940e+01,
  };
/* coeffs for s_1(x) with 1.0 < x < 5.0 */
static double s_1_B[] = 
  {
    9.405607605892576e-14,
    -2.582580614009864e-13,
    8.140558934015147e-13,
    -1.779021374659351e-12,
    3.992281253059363e-12,
    -8.685981227203724e-12,
    1.851504808056151e-11,
    -3.936568243096355e-11,
    8.341933204808976e-11,
    -1.757966701502896e-10,
    3.691209762089122e-10,
    -7.703991787945018e-10,
    1.598620064607365e-09,
    -3.297103519499119e-09,
    6.756686361923544e-09,
    -1.375332155233040e-08,
    2.779594316331730e-08,
    -5.575443395855473e-08,
    1.109456468513641e-07,
    -2.189079621411400e-07,
    4.280593285893081e-07,
    -8.290586328159448e-07,
    1.589385052498797e-06,
    -3.013918355130520e-06,
    5.648806076795733e-06,
    -1.045522108452960e-05,
    1.909188676404741e-05,
    -3.435931158564513e-05,
    6.087054479115233e-05,
    -1.060137750534616e-04,
    1.812439487681223e-04,
    -3.036559465085469e-04,
    4.976134942799733e-04,
    -7.958999393364055e-04,
    1.239413365353734e-03,
    -1.873928850457351e-03,
    2.742139641159466e-03,
    -3.869545770247643e-03,
    5.244460346610938e-03,
    -6.796014304558572e-03,
    8.379251479095341e-03,
    -9.780908388536473e-03,
    1.075827726566235e-02,
    1.064448816293857e+02,
  };
/* coeffs for s_1(x) with 5.0 < x < 30.0 */
static double s_1_C[] = 
  {
    -5.314940406978250e-14,
    1.164926740929414e-13,
    -2.684721132275378e-14,
    7.141761929316007e-14,
    -2.842170943040401e-14,
    -7.896713586061113e-14,
    5.240252676230739e-14,
    -4.553932988280642e-14,
    4.198661620400592e-14,
    1.023827487436144e-13,
    -1.324193280280187e-14,
    -4.521635591200637e-15,
    1.808654236480255e-14,
    -5.006096547400706e-15,
    -1.695613346700239e-14,
    -6.669412497020941e-14,
    8.946378991161262e-14,
    -3.988728539380562e-14,
    1.782816318816251e-13,
    -2.067033413120291e-13,
    5.219259368128736e-13,
    -9.153082332473290e-13,
    1.594522493839825e-12,
    -2.954565884878817e-12,
    5.426608657382365e-12,
    -9.617195928512956e-12,
    1.733498193475084e-11,
    -2.939612190030775e-11,
    5.593037144535629e-11,
    -7.851529527546187e-11,
    1.985569688465396e-10,
    -1.128703595234335e-10,
    9.046116583084023e-10,
    7.142976311446215e-10,
    5.160859860676614e-09,
    8.834179125601209e-09,
    2.891191383272681e-08,
    5.501721809444213e-08,
    1.293205271353641e-07,
    2.229013595922424e-07,
    4.046901619470976e-07,
    5.688908230958780e-07,
    7.876809827419865e-07,
    1.064559931482405e+02,
  };

double s1(double x)
{
  double y,z;
  if ( x< 0.5) {

    z = sqrt(-1);
    
  } else if (x < 1.0) {

    y = (2*x - 1.5)/0.25;
    z = chbevl(y, s_1_A, 40);

  } else if (x < 5.0) {
    
    y = (2*x - 6)/2.0;
    z = chbevl(y, s_1_B, 44);

  } else if (x <= 30.0) {
    
    y = (2*x - 35)/(12.5);
    z = chbevl(y,s_1_C,44);

  } else {

    z = sqrt(-1);

  }

  return z;
}
