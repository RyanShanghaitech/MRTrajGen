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
rho = sp.Function(r'\rho')(t)
tht = (2*sp.pi)/(0.5)*(Np/2)/(u)*rho

# definition of kx and ky
kx = rho*sp.cos(tht + tht0)
ky = rho*sp.sin(tht + tht0)

# slew rate
s2 = kx.diff(t, 2)**2 + ky.diff(t, 2)**2

s2 = s2.expand().collect(rho.diff(t,2))

# clear file
open(filename, "w").close()

# save latex code
dictRep = {
    rho.diff(t,0): sp.Symbol(r'\textcolor{cyan}{\rho_t^{(0)}}'),
    rho.diff(t,1): sp.Symbol(r'\textcolor{cyan}{\rho_t^{(1)}}'),
    rho.diff(t,2): sp.Symbol(r'\textcolor{cyan}{\rho_t^{(2)}}'),
}
save_latex(sp.Eq(sp.Symbol("s")**2, s2.subs(dictRep)))

# save c code
dictRep = {
    Np: sp.Symbol(r'dNp'),
    u: sp.Symbol(r'dU'),
    rho.diff(t,0): sp.Symbol(r'dD0Rho'),
    rho.diff(t,1): sp.Symbol(r'dD1Rho'),
    rho.diff(t,2): sp.Symbol(r'dD2Rho'),
}
save_ccode(sp.Eq(sp.Symbol("double dA"), s2.coeff(rho.diff(t,2),2).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dB"), s2.coeff(rho.diff(t,2),1).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dC"), s2.coeff(rho.diff(t,2),0).subs(dictRep)-sp.Symbol("dS")**2).simplify())
