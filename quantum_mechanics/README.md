# Julia code to visualise Hydrogen orbitals

- ```plotP(n, l, m, rmax, path, steps, c_lim, legend)``` plots a heatmap for the hydrogen orbital with quantum numbers n,l,m up to a distance from the centre rmax (in terms of aₒ). Plot output to given path. The number of steps at which the probability is calculated between the cntre and the radius is defined by the steps argument (default 100). Optionally, the colormap upper limit can be specified with c_lim and colorbar can be toggled on or off through the legend option (default off).
- ```showAll(max, path, steps)``` runs the plotP function for all combinations of n,l,m until n reaches max. Args path and steps same as in plotP.
- ```nlm(max)``` creates an array of all possible (n,l,m) tuples for n up to max.
- ```Hydrogen_ψ(n, l, m, r, θ, ϕ)``` returns the wavefunction for hydrogen for the specified conditions.
- ```Hydrogen_P(n, l, m, r, θ, ϕ)``` returns the probability density for hydrogen for the specified conditions.
- ```Harmonic(l, m, θ, ϕ)``` returns the spherical harmonic function of degree l and order m.
- ```Laguerre(n, a, x)``` general Laguerre polynomial of degree n, alpha a.
- ```Binomial(n, k)``` binomial coefficient of (n,k), accepting n=1/2.
- ```Legendre(l, m, x)``` associated Legendre polynomial of order m.
- ```polar2cartesian(r,θ)``` and ```cartesian2polar(x,y)``` 2D transformation from polar to cartesian coordinates and vice versa.