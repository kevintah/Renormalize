import sympy as sp
import numpy as np
from sympy import diff
import math as math
from builtins import chr
from sympy import init_printing

init_printing(order='none')

x,t,m,k,y, eye = sp.symbols('x t m k y I')
kx,kt = sp.symbols('kx kt')
lbda = sp.Symbol(chr(955))
x_y = sp.Symbol('x-y')
phi = sp.Function(chr(0x03C6))(x,t)
dalembertian_symbol = sp.Symbol(chr(0x25A1))
integral_symbol = sp.Symbol(chr(8747))
kronecker_delta_symbol = sp.Symbol(chr(0x03B4))

xshift = x-y
tshift = t-y

D_x_y = sp.Symbol('D(x-y)')
D_k = sp.Symbol('D(k)')
D_x = sp.Symbol('D(x)')



L = 1/2 * ( phi.diff(t,1))**2 +  1/2 *( phi.diff(x,1) )**2  - 1/2 * m**2 * phi**2 -  (1/math.factorial(4)) * lbda * phi**4

phiDot =  phi.diff(t,1)
phix =  phi.diff(x,1)
dL_dphiDot = L.diff( phiDot , 1 )
dL_dphix = L.diff( phix , 1 )
euler_lagrange_eq =- dL_dphiDot.diff(t,1) - dL_dphix.diff(x,1) + L.diff(phi)
pde = sp.Eq(euler_lagrange_eq, 0)


lbda_to_zero = euler_lagrange_eq.subs(lbda, 0)



# operator method and fourier method



operatorEquation = sp.Eq(dalembertian_symbol - m**2, 1)

#Take fourier tranform on both sides

fourierTransformedOperatorEquation = str(integral_symbol)+ '1/'+str(math.pi**2) +'*'+ str((dalembertian_symbol - m**2) *D_k)+'*' + str( sp.exp(sp.I *k*x))+ '=' + str(integral_symbol)+'1/'+str(math.pi**2)  + str(kronecker_delta_symbol)+'(x)'
fourierTransformedOperatorEquationCompled = str(integral_symbol)+'1/'+str(math.pi**2)+'*'  + str((k**2 - m**2) *D_k)+'*' + str( sp.exp(sp.I *k*x))+ '=' + str(integral_symbol)+'1/'+str(math.pi**2)+'*'  + str(kronecker_delta_symbol)+'(x)'

#Read off D_k

read_Off_D_k =  1/ (math.sqrt(math.pi))**2  * 1/(-k**2 + m**2)
fourier_transformed_read_Off_D_k = str(integral_symbol)+ 'dk'+str(1/math.pi**2)+ '*'+ str( sp.exp(- sp.I *k*x)) + str(1/ (math.sqrt(math.pi))**2  * 1/(-k**2 + m**2))
D_x = fourier_transformed_read_Off_D_k
print(D_x)
















