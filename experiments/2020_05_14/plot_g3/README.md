This data was generated using samslammps repo, commit (this commit was buggy):

commit d5c83d3fbe8a619c6b7f4a651fcd26792b75a965
(HEAD -> steady_press_ABP, origin/steady_press_ABP)
Author: sam <samuel.j.m.cameron@gmail.com>
Date:   Fri May 1 09:07:02 2020 -0300

    added options in three body calculation to minimize data size.

The commit was buggy because it still did not handle the binning of
theta values well, in particular since theta*delinvtheta created
roundoff errors. The data generated in this folder is fine, since
by chance there were no round off errors occuring in 4 of the
14 attempted calculations, but that's why
only 4 of a possible 14 three body correlation functions have been
calculated.