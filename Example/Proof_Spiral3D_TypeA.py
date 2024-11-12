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
u_phi = sp.Symbol(r"u_\phi") # undersamp ratio
u_tht = sp.Symbol(r"u_\theta") # undersamp ratio
phi0 = sp.Symbol(r"\phi_0")
tht0 = sp.Symbol(r"\theta_0")

t = sp.Symbol(r"t")
phi = sp.Function(r"\phi")(t)
tht = sp.sqrt((2*u_phi)/(u_tht))*sp.sqrt(phi)
rho = sp.sqrt(u_tht*u_phi/2)/(sp.pi*Np)*sp.sqrt(phi)

# kx, ky, kz
kx = rho*sp.sin(tht + tht0)*sp.cos(phi + phi0)
ky = rho*sp.sin(tht + tht0)*sp.sin(phi + phi0)
kz = rho*sp.cos(tht + tht0)

# slew rate
s2 = kx.diff(t,2)**2 + ky.diff(t,2)**2 + kz.diff(t,2)**2

s2 = s2.subs({phi.diff(t,2): 0})
s2 = s2.expand().collect(phi.diff(t,1))

# clear file
open(filename, "w").close()

# save latex code
dictRep = {
    phi.diff(t,0): sp.Symbol(r'\textcolor{cyan}{\phi_t^{(0)}}'),
    phi.diff(t,1): sp.Symbol(r'\textcolor{cyan}{\phi_t^{(1)}}'),
}
save_latex(sp.Eq(sp.Symbol("s")**2, s2.subs(dictRep)))

# save c code
dictRep = {
    Np: sp.Symbol(r'dNp'),
    u_phi: sp.Symbol(r'dUPhi'),
    u_tht: sp.Symbol(r'dUTht'),
    phi0: sp.Symbol(r'dPhi0'),
    tht0: sp.Symbol(r'dTht0'),
    phi.diff(t,0): sp.Symbol(r'dD0Phi'),
    phi.diff(t,1): sp.Symbol(r'dD1Phi'),
}
save_ccode(sp.Eq(sp.Symbol("double dA"), s2.coeff(phi.diff(t,1),4).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dB"), s2.coeff(phi.diff(t,1),2).subs(dictRep)).simplify())
save_ccode(sp.Eq(sp.Symbol("double dC"), s2.coeff(phi.diff(t,1),0).subs(dictRep)-sp.Symbol("dS")**2).simplify())
