"""Translations of functions from Release 2.3 of the Cephes Math Library,
(c) Stephen L. Moshier 1984, 1995.
"""

from numpy import exp, floor, log, sin, sqrt


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Sandra Smit", "Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Production"

log_epsilon = 1e-6  # for threshold in log/exp close to 1
# For IEEE arithmetic (IBMPC):
MACHEP = 1.11022302462515654042e-16  # 2**-53
MAXLOG = 7.09782712893383996843e2  # log(2**1024)
MINLOG = -7.08396418532264106224e2  # log(2**-1022)
MAXNUM = 1.7976931348623158e308  # 2**1024

PI = 3.14159265358979323846  # pi
PIO2 = 1.57079632679489661923  # pi/2
PIO4 = 7.85398163397448309616e-1  # pi/4
SQRT2 = 1.41421356237309504880  # sqrt(2)
SQRTH = 7.07106781186547524401e-1  # sqrt(2)/2
LOG2E = 1.4426950408889634073599  # 1/log(2)
SQ2OPI = 7.9788456080286535587989e-1  # sqrt( 2/pi )
LOGE2 = 6.93147180559945309417e-1  # log(2)
LOGSQ2 = 3.46573590279972654709e-1  # log(2)/2
THPIO4 = 2.35619449019234492885  # 3*pi/4
TWOOPI = 6.36619772367581343075535e-1  # 2/pi

ROUND_ERROR = 1e-14  # fp rounding error: causes some tests to fail
# will round to 0 if smaller in magnitude than this


def fix_rounding_error(x):
    """If x is almost in the range 0-1, fixes it.

    Specifically, if x is between -ROUND_ERROR and 0, returns 0.
    If x is between 1 and 1+ROUND_ERROR, returns 1.
    """
    if -ROUND_ERROR < x < 0:
        return 0
    elif 1 < x < 1 + ROUND_ERROR:
        return 1
    else:
        return x


def log_one_minus(x):
    """Returns natural log of (1-x). Useful for probability calculations."""
    if abs(x) < log_epsilon:
        return -x
    else:
        return log(1 - x)


def one_minus_exp(x):
    """Returns 1-exp(x). Useful for probability calculations."""
    if abs(x) < log_epsilon:
        return -x
    else:
        return 1 - exp(x)


def permutations(n, k):
    """Returns the number of ways of choosing k items from n, in order.

    Defined as n!/(n-k)!.
    """
    # Validation: k must be be between 0 and n (inclusive), and n must be >=0.
    if k > n:
        raise IndexError(f"Can't choose {k} items from {n}")
    elif k < 0:
        raise IndexError("Can't choose negative number of items")
    elif n < 0:
        raise IndexError("Can't choose from negative number of items")
    if min(n, k) < 20 and isinstance(n, int) and isinstance(k, int):
        return permutations_exact(n, k)
    else:
        return exp(ln_permutations(n, k))


def permutations_exact(n, k):
    """Calculates permutations by integer division.

    Preferred method for small permutations, but slow on larger ones.

    Note: no error checking (expects to be called through permutations())
    """
    product = 1
    for i in range(n - k + 1, n + 1):
        product *= i
    return product


def ln_permutations(n, k):
    """Calculates permutations by difference in log of gamma function.

    Preferred method for large permutations, but slow on smaller ones.

    Note: no error checking (expects to be called through permutations())
    """
    return lgam(n + 1) - lgam(n - k + 1)


def combinations(n, k):
    """Returns the number of ways of choosing k items from n, in order.

    Defined as n!/(k!(n-k)!).
    """
    # Validation: k must be be between 0 and n (inclusive), and n must be >=0.
    if k > n:
        raise IndexError(f"Can't choose {k} items from {n}")
    elif k < 0:
        raise IndexError("Can't choose negative number of items")
    elif n < 0:
        raise IndexError("Can't choose from negative number of items")
    # if min(n, k) < 20:
    if min(n, k) < 20 and isinstance(n, int) and isinstance(k, int):
        return combinations_exact(n, k)
    else:
        return exp(ln_combinations(n, k))


def combinations_exact(n, k):
    """Calculates combinations by integer division.

    Preferred method for small combinations, but slow on larger ones.

    Note: no error checking (expects to be called through combinations())
    """
    # permutations(n, k) = permutations(n, n-k), so reduce computation by
    # figuring out which requires calculation of fewer terms.
    if k > (n - k):
        larger = k
        smaller = n - k
    else:
        larger = n - k
        smaller = k

    product = 1
    # compute n!/(n-larger)! by multiplying terms from n to (n-larger+1)
    for i in range(larger + 1, n + 1):
        product *= i

    # divide by (smaller)! by multiplying terms from 2 to smaller
    for i in range(2, smaller + 1):  # no need to divide by 1...
        product /= i  # ok to use integer division: should always be factor

    return product


def ln_combinations(n, k):
    """Calculates combinations by difference in log of gamma function.

    Preferred method for large combinations, but slow on smaller ones.

    Note: no error checking (expects to be called through combinations())
    """
    return lgam(n + 1) - lgam(k + 1) - lgam(n - k + 1)


def ln_binomial(successes, trials, prob):
    """Returns the natural log of the binomial distribution.

    successes: number of successes
    trials: number of trials
    prob: probability of success

    Works for int and float values. Approximated by the gamma function.

    Note: no error checking (expects to be called through binomial_exact())
    """
    prob = fix_rounding_error(prob)
    return (
        ln_combinations(trials, successes)
        + successes * log(prob)
        + (trials - successes) * log(1.0 - prob)
    )


# Translations of functions from Cephes Math Library, by Stephen L. Moshier


def polevl(x, coef):
    """evaluates a polynomial y = C_0 + C_1x + C_2x^2 + ... + C_Nx^N

    Coefficients are stored in reverse order, i.e. coef[0] = C_N
    """
    result = 0
    for c in coef:
        result = result * x + c
    return result


# Coefficients for zdist follow:
ZP = [
    2.46196981473530512524e-10,
    5.64189564831068821977e-1,
    7.46321056442269912687e0,
    4.86371970985681366614e1,
    1.96520832956077098242e2,
    5.26445194995477358631e2,
    9.34528527171957607540e2,
    1.02755188689515710272e3,
    5.57535335369399327526e2,
]

ZQ = [
    1.0,
    1.32281951154744992508e1,
    8.67072140885989742329e1,
    3.54937778887819891062e2,
    9.75708501743205489753e2,
    1.82390916687909736289e3,
    2.24633760818710981792e3,
    1.65666309194161350182e3,
    5.57535340817727675546e2,
]

ZR = [
    5.64189583547755073984e-1,
    1.27536670759978104416e0,
    5.01905042251180477414e0,
    6.16021097993053585195e0,
    7.40974269950448939160e0,
    2.97886665372100240670e0,
]
ZS = [
    1.00000000000000000000e0,
    2.26052863220117276590e0,
    9.39603524938001434673e0,
    1.20489539808096656605e1,
    1.70814450747565897222e1,
    9.60896809063285878198e0,
    3.36907645100081516050e0,
]
ZT = [
    9.60497373987051638749e0,
    9.00260197203842689217e1,
    2.23200534594684319226e3,
    7.00332514112805075473e3,
    5.55923013010394962768e4,
]
ZU = [
    1.00000000000000000000e0,
    3.35617141647503099647e1,
    5.21357949780152679795e2,
    4.59432382970980127987e3,
    2.26290000613890934246e4,
    4.92673942608635921086e4,
]


def erf(a):
    """Returns the error function of a: see Cephes docs."""
    if abs(a) > 1:
        return 1 - erfc(a)
    z = a * a
    return a * polevl(z, ZT) / polevl(z, ZU)


def erfc(a):
    """Returns the complement of the error function of a: see Cephes docs."""
    if a < 0:
        x = -a
    else:
        x = a

    if x < 1:
        return 1 - erf(a)

    z = -a * a
    if z < -MAXLOG:  # underflow
        if a < 0:
            return 2
        else:
            return 0
    z = exp(z)

    if x < 8:
        p = polevl(x, ZP)
        q = polevl(x, ZQ)
    else:
        p = polevl(x, ZR)
        q = polevl(x, ZS)

    y = z * p / q

    if a < 0:
        y = 2 - y

    if y == 0:  # underflow
        if a < 0:
            return 2
        else:
            return 0
    else:
        return y


# Coefficients for Gamma follow:
GA = [
    8.11614167470508450300e-4,
    -5.95061904284301438324e-4,
    7.93650340457716943945e-4,
    -2.77777777730099687205e-3,
    8.33333333333331927722e-2,
]

GB = [
    -1.37825152569120859100e3,
    -3.88016315134637840924e4,
    -3.31612992738871184744e5,
    -1.16237097492762307383e6,
    -1.72173700820839662146e6,
    -8.53555664245765465627e5,
]

GC = [
    1.00000000000000000000e0,
    -3.51815701436523470549e2,
    -1.70642106651881159223e4,
    -2.20528590553854454839e5,
    -1.13933444367982507207e6,
    -2.53252307177582951285e6,
    -2.01889141433532773231e6,
]

GP = [
    1.60119522476751861407e-4,
    1.19135147006586384913e-3,
    1.04213797561761569935e-2,
    4.76367800457137231464e-2,
    2.07448227648435975150e-1,
    4.94214826801497100753e-1,
    9.99999999999999996796e-1,
]

GQ = [
    -2.31581873324120129819e-5,
    5.39605580493303397842e-4,
    -4.45641913851797240494e-3,
    1.18139785222060435552e-2,
    3.58236398605498653373e-2,
    -2.34591795718243348568e-1,
    7.14304917030273074085e-2,
    1.00000000000000000320e0,
]

STIR = [
    7.87311395793093628397e-4,
    -2.29549961613378126380e-4,
    -2.68132617805781232825e-3,
    3.47222221605458667310e-3,
    8.33333333333482257126e-2,
]

MAXSTIR = 143.01608
MAXLGM = 2.556348e305
MAXGAM = 171.624376956302725
LOGPI = 1.14472988584940017414
SQTPI = 2.50662827463100050242e0
LS2PI = 0.91893853320467274178

# Generally useful constants
SQRTH = 7.07106781186547524401e-1
SQRT2 = 1.41421356237309504880
MAXLOG = 7.09782712893383996843e2
MINLOG = -7.08396418532264106224e2
MACHEP = 1.11022302462515654042e-16
PI = 3.14159265358979323846

big = 4.503599627370496e15
biginv = 2.22044604925031308085e-16


def igamc(a, x):
    """Complemented incomplete Gamma integral: see Cephes docs."""
    if x <= 0 or a <= 0:
        return 1
    if x < 1 or x < a:
        return 1 - igam(a, x)
    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:  # underflow
        return 0
    ax = exp(ax)
    # continued fraction
    y = 1 - a
    z = x + y + 1
    c = 0
    pkm2 = 1
    qkm2 = x
    pkm1 = x + 1
    qkm1 = z * x
    ans = pkm1 / qkm1

    while 1:
        c += 1
        y += 1
        z += 2
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk / qk
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if abs(pk) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv
        if t <= MACHEP:
            break
    return ans * ax


def igam(a, x):
    """Left tail of incomplete gamma function: see Cephes docs for details"""
    if x <= 0 or a <= 0:
        return 0
    if x > 1 and x > a:
        return 1 - igamc(a, x)

    # Compute x**a * exp(x) / Gamma(a)

    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:  # underflow
        return 0.0
    ax = exp(ax)

    # power series
    r = a
    c = 1
    ans = 1
    while 1:
        r += 1
        c *= x / r
        ans += c
        if c / ans <= MACHEP:
            break

    return ans * ax / a


def lgam(x):
    """Natural log of the gamma fuction: see Cephes docs for details"""
    if x < -34:
        q = -x
        w = lgam(q)
        p = floor(q)
        if p == q:
            raise OverflowError("lgam returned infinity.")

        z = q - p
        if z > 0.5:
            p += 1
            z = p - q
        z = q * sin(PI * z)
        if z == 0:
            raise OverflowError("lgam returned infinity.")
        z = LOGPI - log(z) - w
        return z

    if x < 13:
        z = 1
        p = 0
        u = x
        while u >= 3:
            p -= 1
            u = x + p
            z *= u
        while u < 2:
            if u == 0:
                raise OverflowError("lgam returned infinity.")
            z /= u
            p += 1
            u = x + p
        if z < 0:
            z = -z
        if u == 2:
            return log(z)
        p -= 2
        x = x + p
        p = x * polevl(x, GB) / polevl(x, GC)
        return log(z) + p
    if x > MAXLGM:
        raise OverflowError("Too large a value of x in lgam.")
    q = (x - 0.5) * log(x) - x + LS2PI
    if x > 1.0e8:
        return q
    p = 1 / (x * x)
    if x >= 1000:
        q += (
            (7.9365079365079365079365e-4 * p - 2.7777777777777777777778e-3) * p
            + 0.0833333333333333333333
        ) / x
    else:
        q += polevl(p, GA) / x
    return q


def betai(aa, bb, xx):
    """Returns integral of the incomplete beta density function, from 0 to x.

    See Cephes docs for details.
    """
    if aa <= 0 or bb <= 0:
        raise ValueError("betai: a and b must both be > 0.")
    if xx == 0:
        return 0
    if xx == 1:
        return 1
    if xx < 0 or xx > 1:
        raise ValueError("betai: x must be between 0 and 1.")
    flag = 0
    if (bb * xx <= 1) and (xx <= 0.95):
        t = pseries(aa, bb, xx)
        return betai_result(t, flag)
    w = 1 - xx
    # reverse a and b if x is greater than the mean
    if xx > (aa / (aa + bb)):
        flag = 1
        a = bb
        b = aa
        xc = xx
        x = w
    else:
        a = aa
        b = bb
        xc = w
        x = xx
    if (flag == 1) and ((b * x) <= 1) and (x <= 0.95):
        t = pseries(a, b, x)
        return betai_result(t, flag)
    # choose expansion for better convergence
    y = x * (a + b - 2) - (a - 1)
    if y < 0:
        w = incbcf(a, b, x)
    else:
        w = incbd(a, b, x) / xc
    y = a * log(x)
    t = b * log(xc)
    if ((a + b) < MAXGAM) and (abs(y) < MAXLOG) and (abs(t) < MAXLOG):
        t = pow(xc, b)
        t *= pow(x, a)
        t /= a
        t *= w
        t *= Gamma(a + b) / (Gamma(a) * Gamma(b))
        return betai_result(t, flag)
    # resort to logarithms
    y += t + lgam(a + b) - lgam(a) - lgam(b)
    y += log(w / a)
    if y < MINLOG:
        t = 0
    else:
        t = exp(y)
    return betai_result(t, flag)


def betai_result(t, flag):
    if flag == 1:
        if t <= MACHEP:
            t = 1 - MACHEP
        else:
            t = 1 - t
    return t


incbet = betai  # shouldn't have renamed in first place...


def incbcf(a, b, x):
    """Incomplete beta integral, first continued fraction representation.

    See Cephes docs for details."""

    k1 = a
    k2 = a + b
    k3 = a
    k4 = a + 1
    k5 = 1
    k6 = b - 1
    k7 = k4
    k8 = a + 2

    pkm2 = 0
    qkm2 = 1
    pkm1 = 1
    qkm1 = 1
    ans = 1
    r = 1
    n = 0
    thresh = 3 * MACHEP

    while 1:
        xk = -(x * k1 * k2) / (k3 * k4)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        xk = (x * k5 * k6) / (k7 * k8)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        if qk != 0:
            r = pk / qk
        if r != 0:
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1
        if t < thresh:
            return ans
        k1 += 1
        k2 += 1
        k3 += 2
        k4 += 2
        k5 += 1
        k6 -= 1
        k7 += 2
        k8 += 2

        if (abs(qk) + abs(pk)) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv
        if (abs(qk) < biginv) or (abs(pk) < biginv):
            pkm2 *= big
            pkm1 *= big
            qkm2 *= big
            qkm1 *= big
        n += 1
        if n >= 300:
            return ans


def incbd(a, b, x):
    """Incomplete beta integral, second continued fraction representation.

    See Cephes docs for details."""

    k1 = a
    k2 = b - 1.0
    k3 = a
    k4 = a + 1.0
    k5 = 1.0
    k6 = a + b
    k7 = a + 1.0
    k8 = a + 2.0

    pkm2 = 0.0
    qkm2 = 1.0
    pkm1 = 1.0
    qkm1 = 1.0
    z = x / (1.0 - x)
    ans = 1.0
    r = 1.0
    n = 0
    thresh = 3 * MACHEP

    while 1:
        xk = -(z * k1 * k2) / (k3 * k4)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        xk = (z * k5 * k6) / (k7 * k8)
        pk = pkm1 + pkm2 * xk
        qk = qkm1 + qkm2 * xk
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk

        if qk != 0:
            r = pk / qk
        if r != 0:
            t = abs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        if t < thresh:
            return ans
        k1 += 1
        k2 -= 1
        k3 += 2
        k4 += 2
        k5 += 1
        k6 += 1
        k7 += 2
        k8 += 2

        if (abs(qk) + abs(pk)) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv

        if (abs(qk) < biginv) or (abs(pk) < biginv):
            pkm2 *= big
            pkm1 *= big
            qkm2 *= big
            qkm1 *= big
        n += 1
        if n >= 300:
            return ans


def Gamma(x):
    """Returns the gamma function, a generalization of the factorial.

    See Cephes docs for details."""
    if hasattr(x, "item"):
        # avoid issue of x being a limited precision numpy type
        # use item() method casts to the nearest Python type
        x = x.item()

    sgngam = 1
    q = abs(x)
    if q > 33:
        if x < 0:
            p = floor(q)
            if p == q:
                raise OverflowError("Bad value of x in Gamma function.")
            i = p
            if (i & 1) == 0:
                sgngam = -1
            z = q - p
            if z > 0.5:
                p += 1
                z = q - p
            z = q * sin(PI * z)
            if z == 0:
                raise OverflowError("Bad value of x in Gamma function.")
            z = abs(z)
            z = PI / (z * stirf(q))
        else:
            z = stirf(x)
        return sgngam * z
    z = 1
    while x >= 3:
        x -= 1
        z *= x
    while x < 0:
        if x > -1e9:
            return Gamma_small(x, z)
    while x < 2:
        if x < 1e-9:
            return Gamma_small(x, z)
        z /= x
        x += 1
    if x == 2:
        return float(z)
    x -= 2
    p = polevl(x, GP)
    q = polevl(x, GQ)
    return z * p / q


def Gamma_small(x, z):
    if x == 0:
        raise OverflowError("Bad value of x in Gamma function.")
    else:
        return z / ((1 + 0.5772156649015329 * x) * x)


def stirf(x):
    """Stirling's approximation for the Gamma function.

    Valid for 33 <= x <= 162.

    See Cephes docs for details.
    """
    w = 1.0 / x
    w = 1 + w * polevl(w, STIR)
    y = exp(x)
    if x > MAXSTIR:
        # avoid overflow in pow()
        v = pow(x, 0.5 * x - 0.25)
        y = v * (v / y)
    else:
        y = pow(x, x - 0.5) / y
    return SQTPI * y * w


def pseries(a, b, x):
    """Power series for incomplete beta integral.

    Use when b * x is small and x not too close to 1.

    See Cephes docs for details.
    """
    ai = 1 / a
    u = (1 - b) * x
    v = u / (a + 1)
    t1 = v
    t = u
    n = 2
    s = 0
    z = MACHEP * ai
    while abs(v) > z:
        u = (n - b) * x / n
        t *= u
        v = t / (a + n)
        s += v
        n += 1
    s += t1
    s += ai

    u = a * log(x)
    if ((a + b) < MAXGAM) and (abs(u) < MAXLOG):
        t = Gamma(a + b) / (Gamma(a) * Gamma(b))
        s = s * t * pow(x, a)
    else:
        t = lgam(a + b) - lgam(a) - lgam(b) + u + log(s)
        if t < MINLOG:
            s = 0
        else:
            s = exp(t)
    return s


def log1p(x):
    """Log for values close to 1: from Cephes math library"""
    z = 1 + x
    if (z < SQRTH) or (z > SQRT2):
        return log(z)
    z = x * x
    z = -0.5 * z + x * (z * polevl(x, LP) / polevl(x, LQ))
    return x + z


LP = [
    4.5270000862445199635215e-5,
    4.9854102823193375972212e-1,
    6.5787325942061044846969e0,
    2.9911919328553073277375e1,
    6.0949667980987787057556e1,
    5.7112963590585538103336e1,
    2.0039553499201281259648e1,
]
LQ = [
    1,
    1.5062909083469192043167e1,
    8.3047565967967209469434e1,
    2.2176239823732856465394e2,
    3.0909872225312059774938e2,
    2.1642788614495947685003e2,
    6.0118660497603843919306e1,
]


def expm1(x):
    """Something to do with exp? From Cephes."""
    if (x < -0.5) or (x > 0.5):
        return exp(x) - 1.0
    xx = x * x
    r = x * polevl(xx, EP)
    r /= polevl(xx, EQ) - r
    return r + r


EP = [
    1.2617719307481059087798e-4,
    3.0299440770744196129956e-2,
    9.9999999999999999991025e-1,
]

EQ = [
    3.0019850513866445504159e-6,
    2.5244834034968410419224e-3,
    2.2726554820815502876593e-1,
    2.0000000000000000000897e0,
]


def igami(a, y0):
    # bound the solution
    x0 = MAXNUM
    yl = 0
    x1 = 0
    yh = 1.0
    dithresh = 5.0 * MACHEP

    # handle easy cases
    if (y0 < 0.0) or (y0 > 1.0) or (a <= 0):
        raise ZeroDivisionError("y0 must be between 0 and 1; a >= 0")
    elif y0 == 0.0:
        return MAXNUM
    elif y0 == 1.0:
        return 0.0
    # approximation to inverse function
    d = 1.0 / (9.0 * a)
    y = 1.0 - d - ndtri(y0) * sqrt(d)
    x = a * y * y * y

    lgm = lgam(a)

    for i in range(10):
        # this loop is just to eliminate gotos
        while 1:
            if x > x0 or x < x1:
                break
            y = igamc(a, x)
            if y < yl or y > yh:
                break
            if y < y0:
                x0 = x
                yl = y
            else:
                x1 = x
                yh = y
            # compute the derivative of the function at this point
            d = (a - 1.0) * log(x) - x - lgm
            if d < -MAXLOG:
                break
            d = -exp(d)
            # compute the step to the next approximation of x
            d = (y - y0) / d
            if abs(d / x) < MACHEP:
                return x
            x -= d
            break

    # Resort to interval halving if Newton iteration did not converge.
    d = 0.0625
    if x0 == MAXNUM:
        if x <= 0.0:
            x = 1.0
        while x0 == MAXNUM:
            x = (1.0 + d) * x
            y = igamc(a, x)
            if y < y0:
                x0 = x
                yl = y
                break
            d += d
    d = 0.5
    dir = 0

    for i in range(400):
        x = x1 + d * (x0 - x1)
        y = igamc(a, x)
        lgm = (x0 - x1) / (x1 + x0)
        if abs(lgm) < dithresh:
            break
        lgm = (y - y0) / y0
        if abs(lgm) < dithresh:
            break
        if x <= 0.0:
            break
        if y >= y0:
            x1 = x
            yh = y
            if dir < 0:
                dir = 0
                d = 0.5
            elif dir > 1:
                d = 0.5 * d + 0.5
            else:
                d = (y0 - yl) / (yh - yl)
            dir += 1
        else:
            x0 = x
            yl = y
            if dir > 0:
                dir = 0
                d = 0.5
            elif dir < -1:
                d *= 0.5
            else:
                d = (y0 - yl) / (yh - yl)
            dir -= 1
    if x == 0.0:
        return 0
    return x


P0 = [
    -5.99633501014107895267e1,
    9.80010754185999661536e1,
    -5.66762857469070293439e1,
    1.39312609387279679503e1,
    -1.23916583867381258016e0,
]

Q0 = [
    1.00000000000000000000e0,
    1.95448858338141759834e0,
    4.67627912898881538453e0,
    8.63602421390890590575e1,
    -2.25462687854119370527e2,
    2.00260212380060660359e2,
    -8.20372256168333339912e1,
    1.59056225126211695515e1,
    -1.18331621121330003142e0,
]

s2pi = 2.50662827463100050242e0

P1 = [
    4.05544892305962419923e0,
    3.15251094599893866154e1,
    5.71628192246421288162e1,
    4.40805073893200834700e1,
    1.46849561928858024014e1,
    2.18663306850790267539e0,
    -1.40256079171354495875e-1,
    -3.50424626827848203418e-2,
    -8.57456785154685413611e-4,
]

Q1 = [
    1.00000000000000000000e0,
    1.57799883256466749731e1,
    4.53907635128879210584e1,
    4.13172038254672030440e1,
    1.50425385692907503408e1,
    2.50464946208309415979e0,
    -1.42182922854787788574e-1,
    -3.80806407691578277194e-2,
    -9.33259480895457427372e-4,
]

P2 = [
    3.23774891776946035970e0,
    6.91522889068984211695e0,
    3.93881025292474443415e0,
    1.33303460815807542389e0,
    2.01485389549179081538e-1,
    1.23716634817820021358e-2,
    3.01581553508235416007e-4,
    2.65806974686737550832e-6,
    6.23974539184983293730e-9,
]

Q2 = [
    1.00000000000000000000e0,
    6.02427039364742014255e0,
    3.67983563856160859403e0,
    1.37702099489081330271e0,
    2.16236993594496635890e-1,
    1.34204006088543189037e-2,
    3.28014464682127739104e-4,
    2.89247864745380683936e-6,
    6.79019408009981274425e-9,
]

exp_minus_2 = 0.13533528323661269189


def ndtri(y0):
    """Inverse normal distribution function.

    This is here and not in distributions because igami depends on it..."""
    y0 = fix_rounding_error(y0)
    # handle easy cases
    if y0 <= 0.0:
        return -MAXNUM
    elif y0 >= 1.0:
        return MAXNUM
    code = 1
    y = y0
    if y > (1.0 - exp_minus_2):
        y = 1.0 - y
        code = 0

    if y > exp_minus_2:
        y -= 0.5
        y2 = y * y
        x = y + y * (y2 * polevl(y2, P0) / polevl(y2, Q0))
        x = x * s2pi
        return x

    x = sqrt(-2.0 * log(y))
    x0 = x - log(x) / x

    z = 1.0 / x
    if x < 8.0:  # y > exp(-32) = 1.2664165549e-14
        x1 = z * polevl(z, P1) / polevl(z, Q1)
    else:
        x1 = z * polevl(z, P2) / polevl(z, Q2)
    x = x0 - x1
    if code != 0:
        x = -x
    return x


def incbi(aa, bb, yy0):
    """Incomplete beta inverse function. See Cephes for docs."""
    # handle easy cases first
    if yy0 <= 0:
        return 0.0
    elif yy0 >= 1.0:
        return 1.0

    # define inscrutable parameters
    x0 = 0.0
    yl = 0.0
    x1 = 1.0
    yh = 1.0
    nflg = 0

    if aa <= 1.0 or bb <= 1.0:
        dithresh = 1.0e-6
        rflg = 0
        a = aa
        b = bb
        y0 = yy0
        x = a / (a + b)
        y = incbet(a, b, x)
        return _incbi_ihalve(
            dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0
        )
    else:
        dithresh = 1.0e-4

    # approximation to inverse function
    yp = -ndtri(yy0)

    if yy0 > 0.5:
        rflg = 1
        a = bb
        b = aa
        y0 = 1.0 - yy0
        yp = -yp
    else:
        rflg = 0
        a = aa
        b = bb
        y0 = yy0

    lgm = (yp * yp - 3.0) / 6.0
    x = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0))
    d = yp * sqrt(x + lgm) / x - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (
        lgm + 5.0 / 6.0 - 2.0 / (3.0 * x)
    )
    d *= 2.0
    if d < MINLOG:
        x = 1.0
        return _incbi_under(rflg, x)

    x = a / (a + b * exp(d))
    y = incbet(a, b, x)
    yp = (y - y0) / y0
    if abs(yp) < 0.2:
        return _incbi_newt(
            dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0
        )
    else:
        return _incbi_ihalve(
            dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0
        )


def _incbi_done(rflg, x):
    """Final test in incbi."""
    if rflg:
        if x <= MACHEP:
            x = 1.0 - MACHEP
        else:
            x = 1.0 - x
    return x


def _incbi_under(rflg, x):
    """Underflow handler in incbi."""
    x = 0.0
    return _incbi_done(rflg, x)


class IhalveRepeat(Exception):
    pass

    # Resort to interval halving if not close enough.


def _incbi_ihalve(dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0):
    """Interval halving in incbi."""
    while 1:
        try:
            dir = 0
            di = 0.5
            for i in range(100):
                if i != 0:
                    x = x0 + di * (x1 - x0)
                    if x == 1.0:
                        x = 1.0 - MACHEP
                    if x == 0.0:
                        di = 0.5
                        x = x0 + di * (x1 - x0)
                        if x == 0.0:
                            return _incbi_under(rflg, x)
                    y = incbet(a, b, x)
                    yp = (x1 - x0) / (x1 + x0)
                    if abs(yp) < dithresh:
                        return _incbi_newt(
                            dithresh,
                            rflg,
                            nflg,
                            a,
                            b,
                            x0,
                            yl,
                            x1,
                            yh,
                            y0,
                            x,
                            y,
                            aa,
                            bb,
                            yy0,
                        )
                    yp = (y - y0) / y0
                    if abs(yp) < dithresh:
                        return _incbi_newt(
                            dithresh,
                            rflg,
                            nflg,
                            a,
                            b,
                            x0,
                            yl,
                            x1,
                            yh,
                            y0,
                            x,
                            y,
                            aa,
                            bb,
                            yy0,
                        )
                if y < y0:
                    x0 = x
                    yl = y
                    if dir < 0:
                        dir = 0
                        di = 0.5
                    elif dir > 3:
                        di = 1.0 - (1.0 - di) * (1.0 - di)
                    elif dir > 1:
                        di = 0.5 * di + 0.5
                    else:
                        di = (y0 - y) / (yh - yl)
                    dir += 1
                    if x0 > 0.75:
                        if rflg == 1:
                            rflg = 0
                            a = aa
                            b = bb
                            y0 = yy0
                        else:
                            rflg = 1
                            a = bb
                            b = aa
                            y0 = 1.0 - yy0
                        x = 1.0 - x
                        y = incbet(a, b, x)
                        x0 = 0.0
                        yl = 0.0
                        x1 = 1.0
                        yh = 1.0
                        raise IhalveRepeat
                else:
                    x1 = x
                    if rflg == 1 and x1 < MACHEP:
                        x = 0.0
                        return _incbi_done(rflg, x)
                    yh = y
                    if dir > 0:
                        dir = 0
                        di = 0.5
                    elif dir < -3:
                        di *= di
                    elif dir < -1:
                        di *= 0.5
                    else:
                        di = (y - y0) / (yh - yl)
                    dir -= 1
            if x0 >= 1.0:
                x = 1.0 - MACHEP
                return _incbi_done(rflg, x)
            if x <= 0.0:
                return _incbi_under(rflg, x)
        except IhalveRepeat:
            continue


def _incbi_newt(dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0):
    """Newton's method for incbi."""
    if nflg:
        return _incbi_done(rflg, x)
    nflg = 1
    lgm = lgam(a + b) - lgam(a) - lgam(b)

    for i in range(8):
        # Compute the function at this point.
        if i != 0:
            y = incbet(a, b, x)
        if y < yl:
            x = x0
            y = yl
        elif y > yh:
            x = x1
            y = yh
        elif y < y0:
            x0 = x
            yl = y
        else:
            x1 = x
            yh = y
        if x == 1.0 or x == 0.0:
            break
        # Compute the derivative of the function at this point.
        d = (a - 1.0) * log(x) + (b - 1.0) * log(1.0 - x) + lgm
        if d < MINLOG:
            return _incbi_done(rflg, x)
        if d > MAXLOG:
            break
        d = exp(d)
        # Compute the step to the next approximation of x.
        d = (y - y0) / d
        xt = x - d
        if xt <= x0:
            y = (x - x0) / (x1 - x0)
            xt = x0 + 0.5 * y * (x - x0)
            if xt <= 0.0:
                break
        if xt >= x1:
            y = (x1 - x) / (x1 - x0)
            xt = x1 - 0.5 * y * (x1 - x)
            if xt >= 1.0:
                break
        x = xt
        if abs(d / x) < 128.0 * MACHEP:
            return _incbi_done(rflg, x)
    # Did not converge.
    dithresh = 256.0 * MACHEP
    return _incbi_ihalve(
        dithresh, rflg, nflg, a, b, x0, yl, x1, yh, y0, x, y, aa, bb, yy0
    )
