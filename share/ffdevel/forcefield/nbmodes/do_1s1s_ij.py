from sympy import *
import math

r, pb1, pb2 = symbols('r pb1 pb2')

# DOI: 10.1021/acs.jctc.6b00969, SI
# overlap integral (eq 32,33,34)
# NOTE: we use beta = 1/alpha

pfe1 = exp(-pb1*r)
pfe2 = exp(-pb2*r)
pa1  = 1.0/pb1
pa2  = 1.0/pb2
g0ab = - 4*pa1**2*pa2**2/(pa1**2-pa2**2)**3
g0ba = - 4*pa2**2*pa1**2/(pa2**2-pa1**2)**3
g1ab = pa1/(pa1**2-pa2**2)**2
g1ba = pa2/(pa2**2-pa1**2)**2
s    = 1.0/(8.0*math.pi*r)*( (g0ab+g1ab*r)*pfe1 + (g0ba+g1ba*r)*pfe2)

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
