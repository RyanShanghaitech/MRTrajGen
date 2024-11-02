import sympy as sp

useLaTeX = 1
useCxx = 0

def saveLatex(s):
    f = open("latex.md", "a")
    f.write("$$\n")
    f.write(sp.latex(s) + "\n")
    f.write("$$\n")

if useLaTeX:
    Np = sp.Symbol(r"N_p")
    pi = sp.Symbol(r"\pi")
    t = sp.Symbol(r"t")
    dt = sp.Symbol(r"\Delta t")
    u_phi = sp.Symbol(r"u_\phi") # undersamp ratio
    u_tht= sp.Symbol(r"u_\theta") # undersamp ratio
    phi0 = sp.Symbol(r"\phi_0")
    tht0 = sp.Symbol(r"\theta_0")
else:
    Np = sp.Symbol(r"Np")
    pi = sp.Symbol(r"pi")
    t = sp.Symbol(r"t")
    dt = sp.Symbol(r"dt")
    u_phi = sp.Symbol(r"u_phi") # undersamp ratio
    u_tht= sp.Symbol(r"u_tht") # undersamp ratio
    phi0 = sp.Symbol(r"phi0")
    tht0 = sp.Symbol(r"tht0")

phi = sp.Function(r"\phi")(t)
tht = sp.sqrt((2*u_phi)/(u_tht))*sp.sqrt(phi)
rho = sp.sqrt(u_tht*u_phi/2)/(pi*Np)*sp.sqrt(phi)

kx = rho*sp.sin(tht + tht0)*sp.cos(phi + phi0)
ky = rho*sp.sin(tht + tht0)*sp.sin(phi + phi0)
kz = rho*sp.cos(tht + tht0)

s2 = kx.diff(t,2)**2 + ky.diff(t,2)**2 + kz.diff(t,2)**2

if useLaTeX:
    phi_d0 = sp.Symbol(r"\phi^{(0)}")
    phi_d1 = sp.Symbol(r"\phi^{(1)}")
    phi_d2 = sp.Symbol(r"\phi^{(2)}")
else:
    phi_d0 = sp.Symbol(r"phi_d0")
    phi_d1 = sp.Symbol(r"phi_d1")
    phi_d2 = sp.Symbol(r"phi_d2")
s2 = s2.subs(
    {
        phi:phi_d0,
        phi.diff(t,1):phi_d1,
        phi.diff(t,2):0, # phi_d2,
    },
)
s2 = s2.expand().collect(phi_d1)

s = sp.Symbol(r"s")
a = s2.coeff(phi_d1,4).simplify()
b = s2.coeff(phi_d1,2).simplify()
c = (s2.coeff(phi_d1,0) - s**2).simplify()

if useLaTeX:
    dictRep = {
        phi_d0:sp.Symbol(r"\textcolor{cyan}{\phi^{(0)}}"),
        phi_d1:sp.Symbol(r"\textcolor{cyan}{\phi^{(1)}}"),
        phi_d2:sp.Symbol(r"\textcolor{cyan}{\phi^{(2)}}"),
    }
    saveLatex(sp.Eq(sp.Symbol("s^2"), s2.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("a"), a.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("b"), b.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("c"), c.subs(dictRep)))
    saveLatex("")
elif useCxx:
    print("double a = ", sp.cxxcode(a), ";\n", sep="")
    print("double b = ", sp.cxxcode(b), ";\n", sep="")
    print("double c = ", sp.cxxcode(c), ";\n", sep="")
else:
    print("a = ", sp.pycode(a), sep="")
    print("b = ", sp.pycode(b), sep="")
    print("c = ", sp.pycode(c), sep="")
    