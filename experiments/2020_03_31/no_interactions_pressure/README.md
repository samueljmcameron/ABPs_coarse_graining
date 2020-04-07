NOTE THAT THIS DATA IS NOT FULLY EQUILIBRATED as the diffusion
coefficient is not yet approaching its asymptotic value. This
is seen easily by running the script single_pressure.py for a
specified fp and activity, and you'll see the time behaviour of
the pressure and the diffusion coefficient. It is also notice-
able in the forcefree.py output, as the diffusion coefficients
are slightly below the expectated value.


This data was produced by samslammps github repo, but on the
branch activityFIX, where I am starting to transition from
having the activity as a pair potential to having it be a
separate fix. The commit number is

013dd0fcb44312cf730cc6ae1f01756168f93880

