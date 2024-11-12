import sympy as sp

filename = "result.md"

def save_latex(s:sp.Equality):
    with open(filename, "a") as f:
        f.write("$$\n")
        f.write(sp.latex(s) + "\n")
        f.write("$$\n")

def save_ccode(s:sp.Equality):
    with open(filename, "a") as f:
        f.write("```\n")
        f.write(sp.ccode(s.lhs) + " = ")
        f.write(sp.ccode(s.rhs) + ";\n")
        f.write("```\n")

# symbols and functions
Np = sp.Symbol(r"N_p")
u = sp.Symbol(r"u") # undersamp ratio
tht0 = sp.Symbol(r"\theta_0")

t = sp.Symbol(r't')
tht = sp.Function(r'\theta')(t)
rho = 0.5 - (0.5)/(2*sp.pi)*(u)/(Np/2)*tht

# definition of kx and ky
kx = rho*sp.cos(tht + tht0)
ky = rho*sp.sin(tht + tht0)

# slew rate
s2 = kx.diff(t, 2)**2 + ky.diff(t, 2)**2
 
s2 = s2.expand().collect(tht.diff(t,2))

# clear file
open(filename, "w").close()

# save latex code
dictRep = {
    tht.diff(t,0): sp.Symbol(r'\textcolor{cyan}{\theta_t^{(0)}}'),
    tht.diff(t,1): sp.Symbol(r'\textcolor{cyan}{\theta_t^{(1)}}'),
    tht.diff(t,2): sp.Symbol(r'\textcolor{cyan}{\theta_t^{(2)}}'),
}
save_latex(sp.Eq(sp.Symbol("s")**2, s2.subs(dictRep)))

# save c code
dictRep = {
    Np: sp.Symbol(r'dNp'),
    u: sp.Symbol(r'dU'),
    tht.diff(t,0): sp.Symbol(r'dD0Tht'),
    tht.diff(t,1): sp.Symbol(r'dD1Tht'),
    tht.diff(t,2): sp.Symbol(r'dD2Tht'),
}
save_ccode(sp.Eq(sp.Symbol("double dA"), s2.coeff(tht.diff(t,2),2).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dB"), s2.coeff(tht.diff(t,2),1).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dC"), s2.coeff(tht.diff(t,2),0).subs(dictRep)-sp.Symbol("dS")**2).simplify())
