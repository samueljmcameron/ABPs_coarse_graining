this whole experiment was just testing how lammps calculates force
between two particles. I was confused because the force
calculations I was doing weren't matching up with the dump file.
The reason for this is because my lammps executable was built
with fix_bd_activity.cpp built as having an integration step
after the forces were computed (i.e. I used final_integrate()
in the fix). Changed fix to use initial_integrate() and now
things seem fine.

All data is computed from the lammps executable built with commit:

commit 9af076e952cd121be18a95630b2a6e7e675f24dd
(HEAD -> activityFIX, origin/activityFIX)
Author: sam <samuel.j.m.cameron@gmail.com>
Date:   Mon Apr 6 10:59:51 2020 -0300

