from sympy import *

def prExpr(expr):
    print("$$")
    print_latex(expr)
    print("$$")

# symbols
symT = symbols('t')
symGamma = symbols(r'\gamma')
symSRLim = symbols('s')

# functions
# funA = Function('A')
funTht = Function(r'\theta')(symT)
funRho = Function(r'\rho')(funTht)

# definition of kx and ky
kx = funRho*cos(funTht)
ky = funRho*sin(funTht)

# slew rate
sr = sqrt((1/symGamma*kx.diff(symT, 2))**2 + (1/symGamma*ky.diff(symT, 2))**2)

symD0RhoTht = Symbol(r'\rho_\theta^(0)')
symD1RhoTht = Symbol(r'\rho_\theta^(1)')
symD2RhoTht = Symbol(r'\rho_\theta^(2)')
sr = sr.subs({
    funRho.diff(funTht, 0): symD0RhoTht,
    funRho.diff(funTht, 1): symD1RhoTht,
    funRho.diff(funTht, 2): symD2RhoTht,
})

symD0ThtTime = Symbol(r'\theta_t^(0)')
symD1ThtTime = Symbol(r'\theta_t^(1)')
symD2ThtTime = Symbol(r'\theta_t^(2)')
sr = sr.subs({
    funTht.diff(symT, 0): symD0ThtTime,
    funTht.diff(symT, 1): symD1ThtTime,
    funTht.diff(symT, 2): symD2ThtTime,
})

sr = sr.expand().simplify()
sr = sr.collect(symD2ThtTime)
prExpr(Eq(symSRLim, sr))
