# ECM-factorisation
Python 3 program to factorise a number using the Elliptic Curve Method (ECM)

 This is an experiment to see how quickly a number could be factorised into primes using straightforward readable Python 3. Python 3.10 or higher is needed as the code uses pow(a,-1,n) to perform modular inversion.
 
 In tests with moderately sized factors up to 20 digits this algorithm outperforms pyecm, although I haven't experimented with tuning the parameters of pyecm.

The function factorECM tests primality of a number using Miller-Rabin and factorises a number using trial division and Elliptic Curve Method (ECM). Factorisation may take some time if the number has two large prime factors.


The ECM algorithm uses Montgomery elliptic curves with Suyama parameterisation to construct the curves, and has 2 phases based on the descriptions in [1] and [2]

Phase 1 multiplies the initial point on the elliptic curve by all prime powers pn < B1 to produce a point Q on the curve

Phase 2 calculates multiples pQ for each prime B1 < p ≤ B2
Multiples of Q : [Q,2Q,3Q,..,(D-1)Q] and [DQ,2DQ,3DQ,...], where D is set to √B2 , are pre-calculated and then used to calculate pQ = kDQ + rQ where 0 ≤ r < D
The pre-calculated multiples mQ=(x,z) are converted to (x/z,1) to simplify the calculation of pQ

Bounds B1 and B2 are calculated from the size of the number to be factored, based on the table in [3], and assumes that the number to be factored is a semiprime with 2 similarly sized factors
References

[1] P. Zimmermann. 20 years of ECM. 7th Algorithmic Number Theory Symposium (ANTS VII),2006, Berlin, pp.525–542. inria-00070192v1 https://hal.inria.fr/inria-00070192v1/document

[2] R.Brent, Factorization Of The Tenth Fermat Number, 1999, Mathematics of Computation, Vol. 68, p. 429-451 https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00992-8/S0025-5718-99-00992-8.pdf

[3] README file for GMP-ECM https://github.com/sethtroisi/gmp-ecm/blob/master/README

Typical use:

    from ecmAV import factorECM

    f = factorECM(n) # returns a list of prime factors of n

A website hosting an interactive Javascript implemetation of the algorithm can be found at https://factorisation.netlify.app/
