# Factorisation program with trial division and Elliptic Curve Method (ECM)
# extensively adapted from
# https://www.codingame.com/playgrounds/54090/factorisation-ecm-version-montgomery
# (which is why some comments are in French)

# arthur.vause@gmail.com 

from time import time, perf_counter
from random import randint, randrange
from math import exp, prod
import math
import bisect
from functools import reduce

try:
  #raise Exception("test without GMP")
  from gmpy2 import mpz, invert, gcd, xmpz, gcdext
except:
  from math import gcd
  def invert(a,n): return pow(a,-1,n)
  def mpz(x): return x
  def xmpz(x): return x

# implementations of sieve(n), generator of primes up to and including n
#  with and without gmpy2.
# The gmpy2 version is an order of magnitude faster.

try:
  #raise Exception("test without GMP")
  import gmpy2
  SMALLPRIMELIM = 1_000_000

  # from gmpy2 docs https://gmpy2.readthedocs.io/en/latest/advmpz.html
  def sieve(n):
    '''Returns a generator that yields the prime numbers up to n.'''

    # Increment by 1 to account for the fact that slices  do not include
    # the last index value but we do want to include the last value for
    # calculating a list of primes.
    sieve_limit = gmpy2.isqrt(n) + 1
    n += 1

    # Mark bit positions 0 and 1 as not prime.
    bitmap = gmpy2.xmpz(3)

    # Process 2 separately. This allows us to use p+p for the step size
    # when sieving the remaining primes.
    bitmap[4 : n : 2] =  ~0

    # Sieve the remaining primes.
    for p in bitmap.iter_clear(3, sieve_limit):
        bitmap[p*p : n : p+p] = -1

    return bitmap.iter_clear(2, n)  # clear bits denote a prime
  
except:  
  SMALLPRIMELIM = 200_000
  # Eratosthenes -  generator
  # adapted from https://codereview.stackexchange.com/questions/219998/implementing-sieve-of-eratosthenes-faster
  def sieve(n):
    a = bytearray(n+1)
    a[1:n+1:2] = [1]*len(a[1:n+1:2])
    a[1:3] = [0, 1]
    for i in range(3, int(n**0.5)+1, 2):
      if a[i]:  a[i*i:n+1:2*i] = bytes(len(a[i*i:n+1:2*i]))
    yield 2  
    for i,x in enumerate(a[3:n+1:2]):
      if x: yield 2*i+3


def millerTest(a,d,n,r):
    # test de Miller pour un témoin a
    # Retourne faux si n est composé et vrai si n est probablement premier
    # d et r doivent vérifier n = 2^r * d + 1 avec d impair   
           
    x = pow(a, d, n) # Calcule a^d % n   
    if (x == 1  or x == n-1): 
       return True

    for _ in range(r):    
        x = (x * x) % n 
        if (x == 1):
            
            return False 
        if (x == n-1):
            return True    
    
    return False 

def isPrime(n, k=25): 
    # Test de primalité de Miller Rabin  
    # Si faux alors n est composé et si vrai alors n est probablement premier 
    # k determine le niveau de certitude : P(erreur) < 1/4**k
    
    if (n <= 1 or n == 4):
        return False 
    if (n <= 5):
        return True   
    
    # Trouver d et r tels que n = 2^r * d + 1 avec d impair 
    d = n - 1 
    r = 0
    while (d&1 == 0): 
        d  >>= 1 
        r += 1 
    
    # Effectuer k tests de Miller
    for i in range(k):
        a = randint(2,n-2) 
        if (not millerTest(a, d, n, r)):
              return False  
    return True

def nextPrime(n):
    # premier suivant n
    while not isPrime(n):
        n += 1
    return n


# Addition de deux points d'une courbe de Montgomery
# Calcule P+Q sur base de P, Q et P-Q
def addPoints(xp, zp, xq, zq, xpq, zpq, n):
    u = (xp + zp) * (xq - zq) 
    v = (xp - zp) * (xq + zq) 
    w = (u + v) 
    t = (u - v) 
    w = (w * w) 
    t = (t * t) 
    X = (w * zpq) %n
    Z = (t * xpq) %n
    return X, Z

# add points that are represented as (x::1), i.e. z coordinate = 1
# returns (x,g) an x coordinate and gcd of inversion step
# if gcd > 1, then the inversion has failed, but we have found a factor of n
def addPointsX(xp, xq, xpq, n):
  try:
    Zinv = invert((xq-xp)**2 * xpq, n)
    return ( ((xp*xq -1)**2 * Zinv ) %n, 1 )
  except:
    return (0, gcd((xq-xp)**2 * xpq , n) )


def addPoints2( P, Q, PmQ, n):
    return addPoints(P[0],P[1],Q[0],Q[1],PmQ[0],PmQ[1],n)

# Duplication d'un point d'une courbe de Montgomery
def duplicatePoint(x, z, a, n):
    x2 = x*x 
    z2 = z*z 
    xz = x*z 
    u = (x2-z2)
    Z = ((xz*(x2+a*xz+z2))<<2) %n
    X = u*u %n
    return X, Z

    # (X,Z) <-   ( (x^2-z^2)^2,   4xz(x^2 + axz + z^2))
    #          = ( (x^2-z^2)^2,   4xz( (x-z)^2 + 2xz + axz))
    #          = ( (x^2-z^2)^2,   4xz( (x-z)^2 + 4xz((a+2)/4)))
    

# this isn't used but is included for future use
def duplicatePointA24(x, z, a24, n):
    u, v = x+z, x-z
    u2, v2 = u*u, v*v
    t = u2 - v2
    X = (u2 * v2) %n
    Z = (t * (v2 + a24*t)) %n
    return X, Z

# for testing and development
def samePoint(P,Q,n):
    return P[0]*Q[1] %n == Q[0]*P[1] %n


def montgomeryLadder(k, px, pz, a, n):
    qx, qz = px, pz                                        # Q <= P
    rx, rz = duplicatePoint(px, pz, a, n)                  # R <= 2P
    for c in bin(k)[3:]:
        if c == '1':
            qx, qz = addPoints(rx, rz, qx, qz, px, pz, n)  # Q <= R+Q
            rx, rz = duplicatePoint(rx, rz, a, n)          # R <= 2R
        else:
            rx, rz = addPoints(qx, qz, rx, rz, px, pz, n)  # R <= Q+R
            qx, qz = duplicatePoint(qx, qz, a, n)          # Q <= 2Q
    return qx, qz

#---------------------------------------------------------------------------------
verbose = False # switch for information about the search
def setVerbose( x ):
    global verbose
    verbose=x

timing = False # boolean to switch timing messages
def setTiming( x ):
    global timing
    timing=x

#---------------------------------------------------------------------------------

def stage2preparation(B2,Q,a,n):
  
  start=perf_counter()
  D = int((2*B2)**0.5)+1
  while (D%4 != 0): D+=1 # make sure D is even and a multiple of 4
  
  # array T holds points [None,Q,2Q,3Q,None,5Q,None, ...., (D-1)Q] (odd numbers)
  # array S holds [None, DQ, 2DQ, 3DQ, ..........[(B2//D)*D]Q ]
  # this is messy due to the need to check whether modular inversion
  # is successful at each point, and to calculate the gcd if not

  T = [None]*D
  S = [None]*(B2//D+1)

  try: T[1] = Q[0]*invert(Q[1],n)%n
  except:  return gcd(Q[1],n),D,S,T

  _2Q = duplicatePoint(Q[0],Q[1],a,n)  # = 2Q
  try: T[2] = _2Q[0] * invert(_2Q[1],n)%n
  except:  return gcd(_2Q[1],n),D,S,T

  T[3],g = addPointsX(T[2],T[1],T[1],n)       # = 3Q
  if g>1:  return g,D,S,T

  for d in range(5,D,2):
      T[d],g =  addPointsX(T[d-2], T[2], T[d-4], n)  # = dQ
      if g>1: return g,D,S,T

  assert samePoint( (T[D-1],1), montgomeryLadder(D-1, Q[0],Q[1],a,n),n )

  
  DQ,g = addPointsX(T[D//2 + 1], T[D//2 -1 ], T[2], n) # need D=0 mod 4 for this to work
  if g>1:  return g,D,S,T

  S[1] = DQ

  _2DQ =  duplicatePoint(DQ,1,a,n)  # = 2DQ
  try : S[2] = _2DQ[0] * invert(_2DQ[1],n)%n
  except: return gcd(_2DQ[1],n),D,S,T

  for d in range(3,len(S)):
      S[d],g =  addPointsX(S[d-1], DQ, S[d-2], n)
      if g>1:  return g,D,S,T

  assert samePoint( (S[B2//D],1), montgomeryLadder( (B2//D)*D, Q[0],Q[1],a,n),n )

  if(timing):print("phase 2 S,T construction took",round(perf_counter()-start,3),f"sec for {B2//D+1} S, {D//2} T values");

  return g,D,S,T
#---------------------------------------------------------------------------------
def montgomery(n, B1, B2, primes, count=0, primeqr=[]):

    if verbose: print(f"montgomery n={n}")
    assert(n&1)
    
    n = mpz(n)
    
    # construct a curve and initial point
    # suyama curve from https://members.loria.fr/PZimmermann/papers/ecm.pdf
    start=perf_counter()
    sigma = randint(6, 2**31-1)
    u = (sigma*sigma-5)%n
    v= (4*sigma)%n
    x0 = u*u*u %n
    z0 = v*v*v %n
    try:
      a = (  ((v-u)**3 * (3*u+v)) * pow(4*u*u*u*v,-1,n) - 2 ) %n
    except:
      # pow must have failed, so extract gcd and return
      return int(gcd(4*u*u*u*v,n))
    
    # sanity check - don't understand this. Can't we just set y=1, b=z^-1 mod n ?
    b = (u * pow(z0,-1,n) )%n
    g = gcd(( a * a - 4)* b, n)
    if (g>1):
      return int(g)
    y0 = (sigma*sigma-1)*(sigma*sigma-25)*(sigma**4-25) % n
    assert (b*y0*y0*z0) %n == (x0**3 + a*x0*x0*z0 + x0*z0*z0) %n
    Q = (mpz(x0),mpz(z0))
    if(timing): print(f"Curve construction for n={n} sigma={sigma} took",round(perf_counter()-start,3),"sec");


    # Phase 1
    # it would be marginally faster to pre-calculate the product of primes first
    # than calculating montgomeryLadder for each prime,
    # but that might miss some factors
    start=perf_counter()
    primes_to_B1 = bisect.bisect_right(primes, B1)
    check_interval = int(primes_to_B1**0.5)
    i=0
    for p in primes:
        if p>B1:
            break
        pp = p
        while pp < B1:
            Q = montgomeryLadder(p, Q[0], Q[1], a, n)
            i += 1    
            if i%check_interval==0 :                
                g = gcd(Q[1],n)
                if g>1:
                    if (verbose): print('montgomery phase 1 intermediate found factor g='
                          ,g, "for n=",n,f"after {i} or fewer prime factors out of {primes_to_B1} ",
                          f"with check interval {check_interval}")   
                    return int(g)                   
            pp *= p
    
    g = gcd(Q[1], n)
    
    if 1 < g < n:
        if (verbose): print('montgomery phase 1 found factor g=',g, "for n=",n)    
        return int(g)

    if(timing): print(f"phase 1 for B1={B1} took",round(perf_counter()-start,3),"sec");

    if Q[1]==0 or g!=1:
        print("we have missed a chance to spot a factor somewhere in phase 1\n",
        f"so can't proceed with phase 2, n={n} sigma={sigma} Q={Q} g={g}")
        return int(n)

    """ Phase 2 
    Modifying ideas from

    [1] Paul Zimmermann. 20 years of ECM.
    7th Algorithmic Number Theory Symposium (ANTS VII),2006, Berlin, pp.525–542. inria-00070192v1
    https://hal.inria.fr/inria-00070192v1/document 

    [2] Factorization Of The Tenth Fermat Number, Brent 1999,
    https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00992-8/S0025-5718-99-00992-8.pdf

    Define D = sqrt(2*B2)+1 , and adjust upwards to ensure D is even
    Construct arrays of points S = [0, DQ, 2DQ, 3DQ,.. [B2-B2%D]Q] and T = [0,Q,2Q,3Q 5Q,7Q,...., (D-1)Q]
     ~ sqrt(2*B2) points.

    Represent a prime p as p=kD+r (0<=r<D)
    pQ = S[k]+T[r], so we extract the relevant part of the z coordinate from  S[k]+T[r]

    In a similar manner, [2] pre-calculates [1,2D+1,4D+1,... ] and [2,4,...D], 
    """

    g,D,S,T = stage2preparation(B2,Q,a,n)
    if g>1:
      return g

    # construct list containing (q,r)  for primes where p = q*D+r, 0<=r<D
    if len(primeqr)==0:
      start=perf_counter()
      primeqr += [divmod(p,D) for p in primes[bisect.bisect_left(primes, B1):] ]
      if(timing): print(f"construction of primeqr for n={n} took",perf_counter()-start,f"sec for {len(primeqr)} primes")
    

    # main loop of stage 2
    start=perf_counter()
    
    z = xmpz(1)
    for q,r in primeqr:
      z *= S[q] - T[r] # with xmpz, this is marginally faster than combining * and %
      z %= n

    # Experiments that turned out to be slower than a simple "for" loop
    #z = reduce( (lambda x, y: x * y %n), map(lambda x:Sx[x[0]] - Tx[x[1]] , primeqr))
    #z = reduce( lambda x, y: x*(Sx[y[0]] - Tx[y[1]]) %n , [xmpz(1)]+primeqr )
    if(timing):print("phase 2 main loop took",round(perf_counter()-start,3),f"sec for {len(primeqr)} primes");

    g = gcd(z , n)
    if g>1:
      return int(g)

    return int(n)        
    
    """
    unexpected items to investigate further:
    
    Phase 2 found factor 14240327145353635433586810673273
    for n= 14240327145353635433586810673273
    sigma= 1401904696

    Factor 1 found in construction of Tx
    t=(14663344385599636251700342797170066, 15827680063744479454684689712136160)
    n=89152619386056937386516162741183619
    sigma=543520334 i=706
    T[i]=(14663344385599636251700342797170066, 15827680063744479454684689712136160)
    """
#------------------------------------------------------------------------------

def iterativeECM(factors) :
  # adds prime factors found to the array factors and replaces composites
  
  i=0;
  while ( i < len(factors) ):
    if (isPrime( factors[i] ) ):
      i+=1
    else:
      start = perf_counter()
      p = factors[i]
      
      # B1 based on https://members.loria.fr/PZimmermann/records/ecm/params.html
      bits_in_factor = p.bit_length()//2
      B1 = int(math.exp(.075*bits_in_factor + 5.332))
      B2=200*B1  # B2 set so Stage 2 runs for similar length of time as Stage 1
      
      if (verbose): print(f"factorECM, sieving up to B2={B2}")
      primes = list(sieve(B2))
      if (verbose): print(f"factorECM, finished sieving")
      primesqr = []
      count = 1;
      while (p==factors[i]):
        p = montgomery(factors[i], B1, B2, primes, count, primesqr);
        if (verbose): print("ECM tried", count, "curves, for", factors[i] );
        count+=1
 
      if (timing): print("ECM found factor", p, "for", factors[i],"in", (perf_counter()-start),"sec");
      factors[i]//=p;
      factors.append(p);
 
#--------------------------------------------------------------------------------------   

def factorECM(n):
    global timing, verbose

    if isPrime(n): 
        return [n]
             
    if (verbose): print(f"factorECM, n={n}")

    start=perf_counter()
    factors=[]
    for p in sieve(SMALLPRIMELIM):
        while n%p==0:
            factors.append(p)
            n//=p
    if(timing):print("factorECM small prime trial division took",round(perf_counter()-start,3),"sec");
    if (verbose): print(f"factorECM, factors",factors)

    if n==1:
      return factors
    
    factors.append(n)
    if isPrime(n): 
        return factors
      
    iterativeECM(factors)
    return sorted(factors)

#--------------------------------------------------------------------------------------------------

if __name__ =="__main__":

    """
    n=1010101010101010101010101010101010101
    f = factorECM( n )
    print(n,f)
    """
    
    """
    setTiming(0)
    setVerbose(0)
    trials = [  randint(10**19,10**20)*randint(10**19,10**20)*randint(10**19,10**20)   for i in range(3) ]
    print("\nProducts of random 20 digit numbers:")
    for n in trials:    
        print(f"\n{n} =")
        start= perf_counter()
        f = factorECM( n )
        end= perf_counter()
        assert all(isPrime(x) for x in f) and prod(f) == n
        print(f, round(end-start,1),"sec")
    """
    # multiple rounds of the same size factors
    setVerbose(0)
    setTiming(0)
    start0 = perf_counter()
    print("\nTime trial semiprimes:")
    for _ in range(5):
        i = 20  # number of digits in each prime factor
        p = nextPrime( randrange(10**(i-1), 10**i ) )
        q = nextPrime( randrange(10**(i-1), 10**i ) )
        print(f"{i} digits: ", end='')
        start= perf_counter()
        f = factorECM( p*q )
        end= perf_counter()
        print(p*q,"=",f, round(end-start,1),"sec")
    print("\ntotal time:", perf_counter()-start0,"\n")

    for i in range(17,31):
        p = nextPrime( randrange(10**(i-1),10**i ) )
        q = nextPrime( randrange(10**(i-1),10**i ) )
        print(f"{i} digits: ", end='')
        start= perf_counter()
        f= factorECM(p*q)
        end= perf_counter()
        print(p*q,"=",f, round(end-start,1),"sec")

    trials = [ # example from https://stackoverflow.com/questions/54176869/how-can-i-improve-this-code-of-elliptic-curve-factorization/69762662#69762662
              45733304697523851846830687775886905451041736136239,
              95185297426160521695860779427995281935275,
              97310386288595303632429456851128058791548,
              532123308939657437480746203701315959374550822002194588179897, 
              194684546363820970462807862979880002747247709894380352721688, #small *5601227* p20*p21
              939397082339643433664464465819696119447505323804397628709627  # 3*p25*p35 - leave this for now
              ]



