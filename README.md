# thomson_py
The classical Thomson problem is to find the lowest energy (ground) state of N electrons on the surface of a sphere.
https://en.wikipedia.org/wiki/Thomson_problem

The number of competing minima grows exponentially with N.  This means for even modest N (~100), advanced optimization techniques are required
to determine the global minimum.  Solving this problem is of theoretical interest, but also has practical application to viral capsid and colloidal self assembly problems.

This code extends the Thomson problem to consider points with different charges.  The program random configurations, then uses the basin hopping algorithm to find the global minimum.
