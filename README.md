# thomson_py
The classical Thomson problem is to find the lowest energy (ground) state of N electrons on the surface of a sphere.
https://en.wikipedia.org/wiki/Thomson_problem

The number of metastable/competing ground states grows exponentially with N.  This means for even modest N=200, advanced optimization techniques are required
to determine the ground state.

This small python program generates a random configuration, then uses the basin hopping algorithm to find the global minimum.
