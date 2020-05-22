Computing pressure with the more sophisticated formula in Winkler et al,
and then comparing to theoretical predictions.

All data is computed from the lammps executable built with commit:

commit 88506a4415247f78d32c982169008ca18a83371a
Author: sam <samuel.j.m.cameron@gmail.com>
Date:   Wed Apr 8 10:51:18 2020 -0300

    Created new pair style to compute active pressure (less the swim pressure) using an v_tally style command.