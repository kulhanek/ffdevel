from sympy import *
import math

r, pb1, pb2 = symbols('r pb1 pb2')

# DOI: 10.1021/acs.jctc.6b00969, SI
# overlap integral (eq 36)
# NOTE: we use beta = 1/alpha

pb   = (pb1 + pb2)/2
pr   = pb*r
pfe  = exp(-pr)
s    = pb**3 / (192 * math.pi) *( 3.0 + 3.0*pr +pr**2 ) * pfe

ds   =   diff(s,r)
tt   = - diff(log(s),r) * r
dtt  =   diff(tt,r)

[replacements,  reduced_exprs] = cse( (s,ds,tt,dtt) )

# print variables
print('    real(DEVDP)             :: ', end='')
first = True
for rep in replacements:
    if not first:
        print(',',end='')
    print("{a:}".format(a=str(rep[0])),end='')
    first = False
print()

# print replacements
print()
for rep in replacements:
    print("        {a:3} = {b:}".format(a=str(rep[0]),b=str(rep[1])))

# print results
print()
print("        {a:3} = {b:}".format(a='s',  b=str(reduced_exprs[0])))
print("        {a:3} = {b:}".format(a='ds', b=str(reduced_exprs[1])))
print("        {a:3} = {b:}".format(a='tt', b=str(reduced_exprs[2])))
print("        {a:3} = {b:}".format(a='dtt',b=str(reduced_exprs[3])))
