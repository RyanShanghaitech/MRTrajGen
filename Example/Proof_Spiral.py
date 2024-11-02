from sympy import *

def saveLatex(s:Symbol):
    f = open("proof.md", "a")
    f.write("$$\n")
    f.write(latex(s) + "\n")
    f.write("$$\n")

# symbols and functions
t = symbols(r't')
tht = Function(r'\theta')(t)
rho = Function(r'\rho')(tht)

# definition of kx and ky
kx = rho*cos(tht)
ky = rho*sin(tht)

# slew rate
s2 = kx.diff(t, 2)**2 + ky.diff(t, 2)**2

tht_d0 = Symbol(r'\theta_t^{(0)}')
tht_d1 = Symbol(r'\theta_t^{(1)}')
tht_d2 = Symbol(r'\theta_t^{(2)}')
rho_d0 = Symbol(r'\rho_\theta^{(0)}')
rho_d1 = Symbol(r'\rho_\theta^{(1)}')
rho_d2 = Symbol(r'\rho_\theta^{(2)}')
s2 = s2.subs({
    tht.diff(t, 0): tht_d0,
    tht.diff(t, 1): tht_d1,
    tht.diff(t, 2): tht_d2,
    rho.diff(tht, 0): rho_d0,
    rho.diff(tht, 1): rho_d1,
    rho.diff(tht, 2): rho_d2,
})
s2 = s2.expand().collect(tht_d2)

s = Symbol(r"s")
a = s2.coeff(tht_d2,2).simplify()
b = s2.coeff(tht_d2,1).simplify()
c = s2.coeff(tht_d2,0).simplify()
s2 = a*tht_d2**2 + b*tht_d2 + c

dictRep = {
    tht_d2:Symbol(r"\textcolor{cyan}{\theta^{(2)}}"),
}
saveLatex(Eq(Symbol("s^2"), s2.subs(dictRep)))
saveLatex(Eq(Symbol("a"), a.subs(dictRep)))
saveLatex(Eq(Symbol("b"), b.subs(dictRep)))
saveLatex(Eq(Symbol("c"), (c-s**2).subs(dictRep)))
