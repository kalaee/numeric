Exam in Numerical Methods 2015
===========================================
Student: Alex Arash Sand Kalaee
Student no.: 20112369
Examination date: 2015/06/30 - 2015/07/03
Examination project: no. 24 -- Two-dimensional integrator


Project Description
===========================================
Task A: Two-dimensional iterated integral
-------------------------------------------
1. Implement a two-dimensional integrator for integrals in the form
\int_a^b dx  \int_{d(x)}^{u(x)} dy f(x,y)
which consecutively applies your favorite adaptive one-dimensional
integrator along each of the two dimensions
2. Solve some interesting two-dimensional integrals of the above shape

Task B: Generalisation to higher dimensions
-------------------------------------------
1. Reimplement your integrator from A, this time it must handle arbitrary
number of dimensions, where for each dimension, the limits are functions
of the above dimensions, i.e.:
\int_a^b dx_0 \int_{d(x_0)}^{u(x_0)} dx_1 ...
\int_{d(x_0,...,x_{n-2})}^{u(x_0,...,x_{n-2})} dx_{n-1} f(x)
in n dimensions, where x = (x_0, ..., x_{n-1})
2. Implement a Monte Carlo routine which also solves the above system.
3. Solve some interesting integrals of the above shape in 2, 3 and 4 dimensions
and note have many function calls each routine uses

Task C: Volume of an n-ball
-------------------------------------------
1. Use your routines from B to estimate the volume of a ball in n-dimensions,
where n = 1, ..., 6, and compare both the error of the estimates and number
of function calls each routine needs
2. Further: integrate the function exp(\sum_i x_i) in the n-dimensional ball
and compare the number of function calls


My solution
===========================================
Task A
-------------------------------------------
In my answer to task A, see mainA.c and adapt_2d_speclim.c, i have
applied the adaptive routine with trapezoid quadrature of order 24 to the
very concrete case of two-dimensional integral.
The application of the routine is performed in mainA and the results,
including a description of the considered integrals, are outputtet
to the file A.out

Task B
-------------------------------------------
The n-dim adpative quadrature is generalised by having the
integrator recursively calling itself to the integrator of the next
coordinate, all the way to the last coordinate axis. After integration
in the last coordinate axis the part-integral values propagate back
one level as a function value, until the first coordinate is
finally integrated.

The n-dim Monte Carlo routine has to take account of the fact that
the volume of the integration region is not readily available.
Instead a random point in the first coordinate axis is chosen,
from a uniform distribution between the limits. After this
the limits of the subsequent axis is estimated and a random
point from a uniform distribution is again chosen. Yet, this means
that the overall probability distribution for a point to fall in 
a given point in the integration region is non-uniform, instead it
is governed by the distribution
\rho(x) = V/(h_0*h(x_0)*...*h(x_0,...,x_{n-2}))
where V is the actual volume of the integration region and h(...,x_i)
describes the width of the i+1'th coordinate at that point in the
above coordinate
h(...,x_i) = u(...,x_i) - d(...,x_i)
This Monte Carlo integration is thus an instance of importance sampling
and the estimated integral value, after N samples, is thus
I = h_0 <h(x_0)*...*h(...,x_{i-2})*f(x)> := h_0 <hf>
with variance
\sigma = <hf^2> - <hf>^2
and error estimate like in the book: eq. (10.4)

The application of the above two routines, adapt_nd_speclim and
monte_carlo_speclim is done in mainB and the output is given in B.out

Task C
-------------------------------------------
Both integrators from B are applied the task of estimating the volume
of an n-dimensional ball, radius = 1, by integrating the function
f(x) = 1
over the entire volume, with the limits
d(...,x_i) = -sqrt(1- x_i^2 - ... - x_0^2)
and
u(...,x_i) = -d(...,x_i)
The true volume in n-dimenions is given by
V = \pi^{n/2}/\Gamma(n/2+1)
The number of function calls are shown in C_calls.pdf
and the relative error of the estimates are shown in C_error.pdf
Note that the neccesary number of function calls increases
much faster for adaptive quadratures compared to Monte Carlo
integration. This pattern is the same whether we consider the
integration over the constant function or over the exponential
function.
Lastly we note that the number of samples neccesary in Monte
Carlo integration of the exponential function begins considerably
higher than that of the constant function, yet remains of the same
order of magnitude in dimensions 1 through 6, differently from
the behaviour of the other integrations which increase at
a higher rate.
