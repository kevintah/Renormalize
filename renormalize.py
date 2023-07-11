import sympy as sp
import numpy as np
from sympy import diff
import math as math
from builtins import chr

x,t,m = sp.symbols('x t m')
lbda = sp.Symbol(chr(955))
phi = sp.Function(chr(0x03C6))(x,t)


L = 1/2 * ( phi.diff(t,1))**2 +  1/2 *( phi.diff(x,1) )**2  - 1/2 * m**2 * phi -  (1/math.factorial(4)) * lbda * phi**4

phiDot =  phi.diff(t,1)
phix =  phi.diff(x,1)
dL_dphiDot = L.diff( phiDot , 1 )
dL_dphix = L.diff( phix , 1 )
euler_lagrange_eq =- dL_dphiDot.diff(t,1) - dL_dphix.diff(x,1) + L.diff(phi)
pde = sp.Eq(euler_lagrange_eq, 0)

print(pde)






