import sympy as sp

useLaTeX = 1
useCxx = 0

def saveLatex(s):
    f = open("latex.md", "a")
    f.write("$$\n")
    f.write(sp.latex(s) + "\n")
    f.write("$$\n")

t = sp.Symbol(r"t")
if useLaTeX:
    Np = sp.Symbol(r"N_p")
    pi = sp.Symbol(r"\pi")
    u_phi = sp.Symbol(r"u_\phi") # undersamp ratio
    u_tht= sp.Symbol(r"u_\theta") # undersamp ratio
    phi0 = sp.Symbol(r"\phi_0")
    tht0 = sp.Symbol(r"\theta_0")
elif useCxx:
    Np = sp.Symbol(r"dNp")
    pi = sp.Symbol(r"m_dPi")
    u_phi = sp.Symbol(r"dUPhi") # undersamp ratio
    u_tht= sp.Symbol(r"dUTht") # undersamp ratio
    phi0 = sp.Symbol(r"dPhi0")
    tht0 = sp.Symbol(r"dTht0")
else:
    Np = sp.Symbol(r"Np")
    pi = sp.Symbol(r"pi")
    u_phi = sp.Symbol(r"u_phi") # undersamp ratio
    u_tht= sp.Symbol(r"u_tht") # undersamp ratio
    phi0 = sp.Symbol(r"phi0")
    tht0 = sp.Symbol(r"tht0")

tht = sp.Function(r"\theta")(t)
phi = sp.sqrt(2*u_tht/u_phi)*sp.sqrt(tht)
rho = u_phi/(2*pi*Np)*phi

kx = rho*sp.sin(tht + tht0)*sp.cos(phi + phi0)
ky = rho*sp.sin(tht + tht0)*sp.sin(phi + phi0)
kz = rho*sp.cos(tht + tht0)

s2 = kx.diff(t,2)**2 + ky.diff(t,2)**2 + kz.diff(t,2)**2

if useLaTeX:
    tht_d0 = sp.Symbol(r"\theta^{(0)}")
    tht_d1 = sp.Symbol(r"\theta^{(1)}")
    tht_d2 = sp.Symbol(r"\theta^{(2)}")
if useCxx:
    tht_d0 = sp.Symbol(r"dD0Tht")
    tht_d1 = sp.Symbol(r"dD1Tht")
    tht_d2 = sp.Symbol(r"dD2Tht")
else:
    tht_d0 = sp.Symbol(r"tht_d0")
    tht_d1 = sp.Symbol(r"tht_d1")
    tht_d2 = sp.Symbol(r"tht_d2")
s2 = s2.subs(
    {
        tht:tht_d0,
        tht.diff(t,1):tht_d1,
        tht.diff(t,2):0, # tht_d2,
    },
)
s2 = s2.expand().collect(tht_d1)

if useLaTeX:
    s = sp.Symbol(r"s")
elif useCxx:
    s = sp.Symbol(r"dS")
else:
    s = sp.Symbol(r"s")
a = s2.coeff(tht_d1,4).simplify()
b = s2.coeff(tht_d1,2).simplify()
c = (s2.coeff(tht_d1,0) - s**2).simplify()

if useLaTeX:
    dictRep = {
        tht_d0:sp.Symbol(r"\textcolor{cyan}{\theta^{(0)}}"),
        tht_d1:sp.Symbol(r"\textcolor{cyan}{\theta^{(1)}}"),
        tht_d2:sp.Symbol(r"\textcolor{cyan}{\theta^{(2)}}"),
    }
    saveLatex(sp.Eq(sp.Symbol("s^2"), s2.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("a"), a.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("b"), b.subs(dictRep)))
    saveLatex(sp.Eq(sp.Symbol("c"), c.subs(dictRep)))
    saveLatex("")
elif useCxx:
    print("double dA = ", sp.cxxcode(a), ";", sep="")
    print("double dB = ", sp.cxxcode(b), ";", sep="")
    print("double dC = ", sp.cxxcode(c), ";", sep="")
else:
    print("a = ", sp.pycode(a), sep="")
    print("b = ", sp.pycode(b), sep="")
    print("c = ", sp.pycode(c), sep="")
