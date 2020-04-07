Checking that the implementation I have put into LAMMPS to calculate the
pressure of the ABPs is doing what I think it's doing (it is). The
pressure formula in the current (as of this folder date) version of
my lammps fork is

/       P = (fp*\sum_{i} r_i dot mu_i
/           +\sum_{i,j>i} r_{ij} dot F_{ij}^{lj})/(2*Area).

This is not the correct pressure calculation for ABPs I don't think?
But I will work on implementing the correct version next.

All data is computed from the lammps executable built with commit:

commit 9af076e952cd121be18a95630b2a6e7e675f24dd
(HEAD -> activityFIX, origin/activityFIX)
Author: sam <samuel.j.m.cameron@gmail.com>
Date:   Mon Apr 6 10:59:51 2020 -0300

