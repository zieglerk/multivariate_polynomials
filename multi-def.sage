'''
TOC
1. abbreviation b_{r,n}
2. auxiliary functions u, v, w and their plots
3. generating functions P, R, I, S, Q, E, A and their coefficients
4. approximation functions \eta, \rho, \eps
5. plots of the relative errors
6. tables of explicit q-series
'''

r = var('r')
q = var('q')
z = var('z')
k = var('k')

# 1. abbreviation b_{r,n}
# =======================

def b(r,n):
    return binomial(n+r,r)

# 2. auxiliary functions u, v, w
# ==============================

def u(r,n,k):
    return b(r,k)+b(r,n-k)-2

def v(r,n,s,k):
    return b(r,k)+b(r,n-s*k)-2

def v_plot(r,n_list,s):
    '''
    graph containing v_{r,n,s} for all n in n_list from 1 to n/s
    '''
    G = plot([])
    for n in n_list:
        G += plot(v(r,n,s,k), 1, n/s, color='black')
        # marking integer points
        for x in srange(1,n//s+1):
            G += point([x,v(r,n,s,x)], pointsize=15, color='black')
        G += text("n = "+n.str(), (n/s+0.25, v(r,n,s,n/s)), color='black', fontsize=12)
    # add a little space in x-direction for the last label
    G.axes_range(G.xmin(),G.xmax()+0.1,G.ymin(),G.ymax())
    return G

def w(r,n,k):
    return k*(b(r,n/k)-1)

def w_plot(r,n_list):
    '''
    graph containing w_{r,n} for all *composite* n in n_list from 1 to n.
    '''
    G = plot([])
    for n in n_list:
        if n.is_prime()==false:
            l = n.divisors()[1]
            G += plot(w(r,n,k), l, n, color='black')
            # mark divisors
            for d in divisors(n)[1:]:
                G += point([d, w(r,n,d)], color='black', pointsize=15)
            G += text("n = "+n.str(), (n+0.75, w(r,n,n)), color='black', fontsize = 12)
    # add a little space in x-direction for the last label
    G.axes_range(G.xmin(),G.xmax()+0.6,G.ymin(),G.ymax())
    return G

# 3. generating functions P, R, I, S, Q, E, A and their coefficients
# ==================================================================

def P(r, n, q):
    '''number of monic polynomials in r variables at degree n over GF(q); denoted \# P_{r,n}(q) in the accompanying paper.
    '''
    return q^(b(r,n)-1)*(1-q^(-b(r-1,n)))/(1-q^(-1))

def PP(r,z,N):
# generating function for P truncated at degree N
    return sum(P(r,i,q)*z^i for i in range(0,N+1)).simplify_full()

def II(r,z,N):
# generating function for irreducible monic polynomials up to and including degree N
    S = sum(moebius(k)/k*log(PP(r,z^k,N)) for k in range(1,N+1)) # check that looking up to degree N is enough
    S = S.taylor(z,0,N)
    return S.simplify_full()

def RR(r,z,N):
# generating function for reducible monic polynomials up to and including degree N
    S = PP(r,z,N) - II(r,z,N)
    S = S.taylor(z,0,N)
    return S.simplify_full()

def SS(r,s,z,N):
    '''generating function for s-powerful monic polynomials up to and including degree N'''
    S = PP(r,z,N)/PP(r,z^s,N) # is it enough to look up to degree N in the denominator
    S = S.taylor(z,0,N)
    return S.simplify_full()

def QQ(r,s,z,N):
# generating function for s-power*less* monic polynomials up to and including degree N
    S = PP(r,z,N) - SS(r,s,z,N)
    S = S.taylor(z,0,N)
    return S.simplify_full()

def A(r, n, q):
    '''number of absolutely irreducible monic polynomials in r variables at degree n over GF(q)'''
    if n == 0:
        return 0
    S = II(r,z,n)
    return sum(sum(moebius(s)*S.coeff(z^(n/k)).substitute(q=q^s) for s in divisors(k))/k for k in divisors(n))

# TODO this function call can be optimized.  A -- and therefore II is called to often.

def AA(r,z,N):
    '''generating function for A truncated at degree N.'''
    return sum(A(r,i,q)*z^i for i in range(0,N+1)).simplify_full()

def EE(r,z,N):
    '''generating function for the number of relatively irreducible monic polynomials truncated at degree N'''
    return II(r,z,N)-AA(r,z,N)

# extract the coefficients by AA.coeff(z^n).simplify_full()

# 4. approximation functions \rho, \eta, \eps
# ===========================================

def rho(r,n,q):
    S = P(r,1,q)*P(r,n-1,q)/(1-q^(-b(r-1,n-1)))
    return S.simplify_full()

def eta(r,n,s,q):
    S = P(r,1,q)*P(r,n-s,q)
    return S.simplify_full()

def eps(r,n,q):
    ell = n.divisors()[1]
    S = P(r,n/ell,q^ell)/(ell*(1-q^(-ell*b(r-1,n/ell))))
    return S.simplify_full()

# 5. plotting the relative error
# ==============================

def RR_vs_rho_plot(r,n_list):
    '''
    graph containing relative errors (R_{r,n}(q)-rho_{r,n}(q))/rho_{r,n}(q) and the upper bound of
    '''
    upper_bound = 1/((1-q^(-1))*(1-q^(-r)))
    G = plot(upper_bound, 2, 20, color='black')
    G += text("$1/((1-\mathbf{q}^{-r})(1-\mathbf{q}^{-1}))$", (10, 1.5), color='black', fontsize=12)
    print 'calling generating function up to N =', max(n_list)
    S = RR(r,z,max(n_list))
    for n in n_list:
        print 'generating graph for n =', n
        numerator = expand(S.coeff(z^n) - rho(r,n,q))
        denominator = expand( rho(r,n,q) * q^(-binomial(n+r-2,r-1)+r*(r+1)/2))
        quotient = (numerator/denominator).simplify_full()
        X = srange(2,20,0.1)
        Y = [quotient(x).n() for x in X]
        G += list_plot(zip(X,Y), plotjoined=True, color='black')
        G += text("$n=%s$"%str(n), (1, quotient(2)), fontsize=12, color='black')
    return G

def QQ_vs_eta_plot(r,n_list,s):
    '''
    graph containing relative errors (Q_{r,n,s}(q)-eta_{r,n,s}(q))/eta_{r,n,s}(q) and the upper bound of
    '''
    G = plot([])
    print 'calling generating function up to N =', max(n_list)
    S = QQ(r,s,z,max(n_list))
    for n in n_list:
        print 'generating graph for n =', n
        numerator = expand(S.coeff(z^n) - eta(r,n,s,q))
        denominator = expand( eta(r,n,s,q) * q^(-binomial(n-s+r,r)+binomial(n-2*s+r,r)+r*(r+1)/2))
        quotient = (numerator/denominator).simplify_full()
        X = srange(2,20,0.1)
        Y = [quotient(x).n() for x in X]
        G += list_plot(zip(X,Y), plotjoined=True, color='black')
        G += text("$n=%s$"%str(n), (1, quotient(2)), fontsize=12, color='black')
    return G

def EE_vs_eps_plot(r,n_list):
    '''
    graph containing relative errors (E_{r,n}(q)-eps_{r,n}(q))/eta_{r,n}(q) and the upper bound of
    '''
    G = plot([])
    print 'calling generating function up to N =', max(n_list)
    S = EE(r,z,max(n_list))
    for n in n_list:
        print 'generating graph for n =', n
        ell = n.divisors()[1]
        numerator = expand(S.coeff(z^n) - eps(r,n,q))
        denominator = expand( eps(r,n,q) * q^(-(ell-1)*(b(r-1,n/ell)-r)-1))
        quotient = (numerator/denominator).simplify_full()
        X = srange(2,10,0.1)
        Y = [quotient(x).n() for x in X]
        G += list_plot(zip(X,Y), plotjoined=True, color='black')
        G += text("$n=%s$"%str(n), (1, quotient(2)), fontsize=12, color='black')
    return G


# 6. tables of explicit q-series
# ==============================

def R_table(r,N):
    S = RR(r,z,N)
    # start of the table
    t  = [r"\begin{tabular}{>{\centering\arraybackslash}p{0.07\textwidth}p{0.86\textwidth}}"]
    t.append(r"\toprule")
    # first row
    t.append(r" $n$  & $\# \Rone{{ {0},n }} (\FF_{{q}})$ \\".format(r))
    t.append(r"\midrule")
    t.append(r" 1  & $0$ \\")
    # further rows
    for n in range(2,N+1):
        t.append(r"{0} & $({1})/{2}$ \\".format(n, latex(n*S.coeff(z^n)), n))
    t.append(r"\bottomrule")
    # add the last line and return
    t.append(r"\end{tabular}")
    return ''.join(t)

def Q_table(r,N,s):
    S = QQ(r,s,z,N)
    # start of the table
    t  = [r"\begin{tabular}{>{\centering\arraybackslash}p{0.07\textwidth}p{0.86\textwidth}}"]
    t.append(r"\toprule")
    # first row
    t.append(r" $n$  & $\# \Qone{{ {0},n,{1} }} (\FF_{{q}})$ \\".format(r,s))
    t.append(r"\midrule")
    t.append(r" {0}  & $0$ \\".format(srange(s).__str__()[1:-1]))
    # further rows
    for n in range(s,N+1):
        t.append(r"{0} & $ {1} $ \\".format(n, latex(S.coeff(z^n))))
    t.append(r"\bottomrule")
    # add the last line and return
    t.append(r"\end{tabular}")
    return ''.join(t)

def E_table(r,N):
    S = EE(r,z,N)
    # start of the table
    t  = [r"\begin{tabular}{>{\centering\arraybackslash}p{0.07\textwidth}p{0.86\textwidth}}"]
    t.append(r"\toprule")
    # first row
    t.append(r" $n$  & $\# \Eone{{ {0},n }} (\FF_{{q}})$ \\".format(r))
    t.append(r"\midrule")
    t.append(r" 1  & $0$ \\")
    # further rows
    for n in range(2,N+1):
        t.append(r"{0} & $({1})/{2}$ \\".format(n, latex(n*S.coeff(z^n)), n))
    t.append(r"\bottomrule")
    # add the last line and return
    t.append(r"\end{tabular}")
    return ''.join(t)

# checking cor 5.24 for 5 > n >= s >= 2
'''
n = 2
s = 2
N = QQ(r,s,z,n)
N = N.coeff(z^n).simplify_full()
D = RR(r,z,n)
D = D.coeff(z^n).simplify_full()
quot = N/D
modquot = quot/q^(-binomial(r+n-1,r)+binomial(r+n-s,r))
modquot = modquot.simplify_full()
left = 1/6
right = 19
p= plot(left,2,20, color='red')
p+= plot(right,2,20, color='red')
for R in srange(2,6):
    q2 = modquot(r=R)
    q2 = q2.simplify_full()
    p += plot(q2,2,20)
p.show()
'''

# checking Cor 6.36 for 5 > n >= 2
'''
n = 4
ell = n.divisors()[1]
N = EE(r,z,n)
N = N.coeff(z^n).simplify_full()
D = II(r,z,n)
D = D.coeff(z^n).simplify_full()
quot = N/D
modquot = quot/q^(-binomial(r+n,r)+ell*binomial(r+n/ell,r)-ell+1)
modquot = modquot.simplify_full()
left = 1/(8*ell)
right = 2/ell
p= plot(left,2,20, color='red')
p+= plot(right,2,20, color='red')
for R in srange(2,6):
    q2 = modquot(r=R)
    q2 = q2.simplify_full()
    p += plot(q2,2,20)
p.show()
'''
