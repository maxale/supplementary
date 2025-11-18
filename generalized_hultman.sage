'''
Supplementary SageMath code ver. 20251104 for the paper:
* N. Alexeev, A. Pologova, M. A. Alekseyev, Generalized Hultman Numbers and Cycle Structures of Breakpoint Graphs, Journal of Computational Biology 24:2 (2017), 93-105.
  DOI: http://doi.org/10.1089/cmb.2016.0190
  Preprint: arXiv:1503.05285 [q-bio.GN] https://doi.org/10.48550/arXiv.1503.05285
'''

# Base ring in variable u
K.<u> = PolynomialRing(QQ)

def L0(f, n):
    # print("L0:",f,n)
    s = f.parent().gens()
    return sum(sum((i-1)*s[j-1]*s[i-j-1]*derivative(f, s[i-2]) for j in (1..i-1)) for i in (2..n))


def L1(f, n):
    # print("L1:",f,n)
    s = f.parent().gens()
    return sum((i-1)^2 * s[i-1] * derivative(f, s[i-2]) for i in (2..n))


def L2(f, n):
    # print("L2:",f,n)
    s = f.parent().gens()
    return sum(s[i] * sum(j*(i-j)*derivative(f, s[j-1], s[i-j-1]) for j in (1..i-1)) for i in (2..n))


def Ln(f, n):
    # print("L3:",f,n)
    s = f.parent().gens()
    return sum((i - 1)*u*s[i-1]*derivative(f, s[i-2]) for i in (2..n))


def FG(n, orient, multichr):
    s = polygens(K,n+1,'s')
    ff = s[0]
    for k in (1..n-1):
        ff = 1/k*(L0(ff, n) + orient*L1(ff, n) + (1 + orient)*L2(ff, n) + multichr*Ln(ff, n))
    return ff


# function G_n(0;s1,s2,...)     (co-oriented unichromosomal)
def G0(n):
    return FG(n, 0, 0).specialization({u:0})


# function G_n(u;s1,s2,...)     (co-oriented multichromosomal)
def Gu(n):
    return FG(n, 0, 1)


# function F_n(0;s1,s2,...)     (arbitrarily-oriented unichromosomal)
def F0(n):
    return FG(n, 1, 0).specialization({u:0})


# function F_n(u;s1,s2,...)     (arbitrarily-oriented multichromosomal)
def Fu(n):
    return FG(n, 1, 1)


# Hultman numbers
def H(n, m):
    f = G0(n+1)
    x = polygen(QQ)
    return ZZ( f.substitute({v: x for v in f.variables()})[m] )


# signed Hultman numbers H^{\pm}(n,m)
def Hs(n, m):
    f = F0(n+1)
    x = polygen(QQ)
    return ZZ( f.substitute({v: x for v in f.variables()})[m] )


# H^1_2(n,d)
def H12(n, d):
    return Hs(n-1, n-d)


# signed/unsigned H^h_3(n,d)
def Hh3(h, n, d, signed=True):
    f = Fu(n) if signed else Gu(n)
    f = f.map_coefficients(lambda c: c[h-1]).specialization({u:0})    # # coefficient of u^(h-1)
    s = f.parent().gens()
    x = polygen(QQ)
    f = f.substitute({s[i-1]:x^(i%2) for i in (1..n+1)})
    return ZZ( f[n-2*d] )       # coefficient of x^(n-2*d)


# signed/unsigned H^h_4(n,d)
def Hh4(h, n, d, signed=True):
    f = Fu(n) if signed else Gu(n)
    f = f.map_coefficients(lambda c: c[h-1]).specialization({u:0})    # # coefficient of u^(h-1)
    s = f.parent().gens()
    x = polygen(K)
    f = f.substitute({s[i-1]:x^((3-i%3)%3) for i in (1..n+1)}).change_ring(ZZ)
    return sum(sum(f[2*(n+i-3*d)+j] for i in range(3)) for j in range(2))


# A264614 Irregular triangle read by rows: T(n,k) = number of unsigned unichromosonal genomes with n genes at 3-break distance k from a fixed genome, 0 <= k <= floor(n/2).
def a264614_row(n):
    return [Hh3(1, n, d, False) for d in range(n//2 + 1)]

# A264615 Irregular triangle read by rows: T(n,k) = number of signed unichromosonal genomes with n genes at 3-break distance k from a fixed genome, 0 <= k <= floor(n/2).
def a264615_row(n):
    return [Hh3(1, n, d, True) for d in range(n//2 + 1)]

# A264616 Irregular triangle read by rows: T(n,k) = number of unsigned unichromosonal genomes with n genes at 4-break distance k from a fixed genome, 0 <= k <= floor((n+1)/3).
def a264616_row(n):
    return [Hh4(1, n, d, False) for d in range((n+1)//3 + 1)]

# A264617 Irregular triangle read by rows: T(n,k) = number of signed unichromosonal genomes with n genes at 4-break distance k from a fixed genome, 0 <= k <= floor((n+1)/3).
def a264617_row(n):
    return [Hh4(1, n, d, True) for d in range((n+1)//3 + 1)]

# A132803 Number of simple permutations in S_n, i.e., those whose "cycle graph" (or "breakpoint graph") only contains alternating cycles of length at most 3.
def a132803(n):
    #    f = Gu(n).specialization({u:1})        # use this in multichromosomal case
    f = G0(n+1)
    return f( *(ZZ(i<=2) for i in range(f.parent().ngens())) )

