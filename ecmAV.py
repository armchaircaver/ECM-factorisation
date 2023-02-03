# extensively adapted from
# https://www.codingame.com/playgrounds/54090/factorisation-ecm-version-montgomery
# (which is why some comments are in French)
# arthur.vause@gmail.com

from time import time, perf_counter
from random import randint, randrange
from math import exp
from array import array
import bisect
try:
  #raise Exception("test without GMP")
  from gmpy2 import mpz, invert, gcd
except:
  from math import gcd
  def invert(a,n): return pow(a,-1,n)
  def mpz(x): return x

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

# Eratosthenes - returns an array b where b[i]=1 iff i is prime
def sieve(n):
    start=perf_counter()
    b = array('b',[1])*(n+5)
    b[0]=0
    b[1]=0
    for p in range(2, int(n**0.5)+1):
        if b[p]:
            b[p*p:n+1:p] = array('b',[0])*len(b[p*p:n+1:p])
    return b 


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
    qx, qz = px, pz                                        # Q = P
    rx, rz = duplicatePoint(px, pz, a, n)                  # R = 2P
    for c in bin(k)[3:]:
        if c == '1':
            qx, qz = addPoints(rx, rz, qx, qz, px, pz, n)  # Q = R+Q
            rx, rz = duplicatePoint(rx, rz, a, n)          # R = 2R
        else:
            rx, rz = addPoints(qx, qz, rx, rz, px, pz, n)  # R = Q+R
            qx, qz = duplicatePoint(qx, qz, a, n)          # Q = 2Q
    return qx, qz

verbose = False # switch for information about the search
def setVerbose( x ):
    global verbose
    verbose=x

timing = False # boolean to switch timing messages
def setTiming( x ):
    global timing
    timing=x

def montgomery(n, B1, B2, primes, count=0, primeqr=[]):

    n = mpz(n)
    
    # construct a curve and initial point
    # suyama curve from https://members.loria.fr/PZimmermann/papers/ecm.pdf
    start=perf_counter()
    if count==1:
        sigma=11
    else:
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
    if(timing): print("Curve construction took",round(perf_counter()-start,3),"sec");


    # Phase 1
    # it would be marginally faster to pre-calculate the product of primes first
    # than calculating montgomeryLadder for each prime,
    # but that might miss some factors
    start=perf_counter()
    primes_to_B1 = bisect.bisect_left(primes, B1)
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

    if(timing): print("phase 1 took",round(perf_counter()-start,3),"sec");

    if Q[1]==0 or g!=1:
        print("we have missed a chance to spot a factor somewhere in phase 1")
        print("so can't proceed with phase 2")
        print(f"n={n} sigma={sigma} Q={Q} g={g}")
        return int(n)

    def normalise(P):
        # convert point (x,z) to (x/z,1)
        return ( P[0]*invert(P[1],n)%n, 1 )

    """ Phase 2 
    Modifying ideas from [1] Paul Zimmermann. 20 years of ECM.
    7th Algorithmic Number Theory Symposium (ANTS VII),2006, Berlin, pp.525–542. inria-00070192v1
    https://hal.inria.fr/inria-00070192v1/document

    Define D = sqrt(B2)+1
    Construct arrays of points S = [0, DQ, 2DQ, 3DQ,.. ] and T = [0, Q, 2Q, 3Q, ...., (D-1)Q]

    represent a prime p as p=kD+r (0<=r<D)
    pQ = S[k]+T[r], so we extract the x coordinate from  S[k]+T[r]

    Using this strategy, we need to calculate more points than in the improved standard
    continuation in [2] Factorization Of The Tenth Fermat Number, Brent 1999,
    https://www.ams.org/journals/mcom/1999-68-225/S0025-5718-99-00992-8/S0025-5718-99-00992-8.pdf
    which pre-calculates [1,2D+1,4D+1,... ] and [2,4,...D],
    We calculate 2*sqrt(B2) vs sqrt(2*B2) of Brent.

    However the main loop is better here, as we don't calculate the main loop for composites as in [1],
    or have a more complex calculation of to determine the elements of arrays to add, i.e. calculate
    (p-1)/2 to determine the array items as in [2]
    """
    start=perf_counter()
    D = int(B2**0.5)+1
    while (D%4 != 0): D+=1 # make sure D is even and a multiple of 4
    
    # array T holds points [-, Q, 2Q, 3Q,-,5Q, ...., (D-1)Q] (odd numbers)
    T = ["not a number"]*D
    T[1]= [mpz(Q[0]),mpz(Q[1])]
    T[2]= duplicatePoint(Q[0],Q[1],a,n)  # = 2Q
    T[3]= addPoints2(T[2],Q,Q,n)       # = 3Q

    for d in range(5,D,2):
        T[d] =  addPoints2(T[d-2], T[2], T[d-4], n)  # = dQ


    Tx = ["not a number"]*D
    for i in range(1,D,2):
        try:
            Tx[i] = ( (T[i][0]*invert(T[i][1],n))% n )
        except:
            g = gcd(T[i][1],n)
            if (verbose): print(f"Factor {g} found in construction of Tx")
            if (verbose): print(f"t={t} n={n} sigma={sigma} i={i} T[i]={T[i]}")
            return int(g)


    DQ = addPoints2(T[D//2 + 1],T[D//2 -1 ],T[2],n) # need D=0 mod 4 for this to work
    #DQ2 = montgomeryLadder(D, Q[0], Q[1], a, n)
    #assert samePoint(DQ,DQ2,n)
  
    #array S holds [0, DQ, 2DQ, 3DQ, .......... ]
    S = [False]*(D+2)
    S[0] = (0,0)
    S[1] = ( mpz(DQ[0]),mpz(DQ[1]) )
    S[2] =  duplicatePoint(DQ[0],DQ[1],a,n)  # = 2DQ
    for d in range(3,len(S)):
        S[d] =  addPoints2(S[d-1], DQ, S[d-2], n)  

    Sx = ["not a number"]
    for s in S[1:]:
        try:
            Sx.append( (s[0]*invert(s[1],n))%n )
        except:
            g = gcd(s[1],n)
            if (verbose): print(f"Factor {g} found in construction of Sx")
            if (verbose): print(f"s={s} n={n} sigma={sigma}")
            return int(g)


    if(timing):print("phase 2 S,T construction took",round(perf_counter()-start,3),"sec");

    # construct list containing (q,r)  for primes where p = q*D+r, 0<=r<D
    start=perf_counter()
    if len(primeqr)==0:
      startindex = bisect.bisect_left(primes, B1)
      for p in primes[startindex:] :
        primeqr.append ( divmod(p,D) )
      
      if(timing): print("construction of primeqr took",perf_counter()-start,"sec")
    

    # main loop
    start=perf_counter()
    
    z = mpz(1)
    for q,r in primeqr:
      z = (z * (Sx[q] - Tx[r])) %n
      
    g = gcd(z , n)
    if g>1:
      if (verbose): print("Phase 2 found factor",g,"for n=",n,"sigma=",sigma)
      return int(g)

    if(timing):print("phase 2 main loop took",round(perf_counter()-start,3),"sec");
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

    

def factorECM(n, k=25, primes=None):
    global timing, verbose
    
    if isPrime(n,k): 
        return [n]

    # curve fitting B1,B2,curves from data in
    # https://github.com/sethtroisi/gmp-ecm/blob/master/README
    # assuming n is a semiprime with 2 similarly sized factors
    B1 = int( 28.6676* exp(0.147729* len(str(n))))
    B2 =int( 326.260 * exp(0.211736*len(str(n))))
    curves = int( 3.76045 * exp(0.080348*len(str(n))))
             
    if (verbose): print(f"montgomery, n={n} B1={B1} B2={B2}")

    start=perf_counter()         
    if primes is None:
        sieveB2 = sieve(B2)
        primes = [p for p in range(B2) if sieveB2[p]]
    if(timing): print("sieve took",round(perf_counter()-start,3),"sec");

    start=perf_counter()
    for p in primes:
        if n%p==0:
            return [p]+factorECM(n//p, k, B1, B2, primes)
    if(timing):print("small prime trial division took",round(perf_counter()-start,3),"sec");

    g = 1
    i = 0
    start=perf_counter()
    primeqr = []
    while g==1 or g==n: # tant qu'on n'a pas de facteur
        # on essaye avec d'autres valeurs
        i +=1
        g=montgomery(n, B1, B2, primes, i, primeqr)  
    if (verbose): print(f"{i} curves needed")
                       
    return factorECM(g, k, primes)+factorECM(n//g, k, primes)

if __name__ =="__main__":
    """
    for i in range(2):
        a=nextPrime(randint(1e8,1e12))    
        b=nextPrime(randint(1e8,1e10))
        c=randint(1e8,1e10) 
        n=a*b*c
        print(f'n={a}*{b}*{c}=',n)
        t=time()
        f=factorECM(n)  
        t=time()-t  
        print('les facteurs de n sont',f,'t=',round(t,1),'s')
        print()
    """

    # multiple rounds of the same size factors
    setVerbose(False)
    setTiming(True)
    start0 = perf_counter()
    for _ in range(10):
        i = 18  # number of digits in each prime factor
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




