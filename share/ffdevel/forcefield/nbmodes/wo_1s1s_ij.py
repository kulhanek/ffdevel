from sympy import *
import math

r, pb1, pb2 = symbols('r pb1 pb2')

# DOI: 10.1063/1.5081060
# WF overlap integral (eq 27)

X    = (pb1/2)**2 - (pb2/2)**2
pr1  = pb1*r/2.0
pfe1 = exp(-pr1)
pr2  = pb2*r/2.0
pfe2 = exp(-pr2)
sce  = pb1*(r*X-2*pb2)*pfe2 + pb2*(r*X+2*pb1)*pfe1
s    = sqrt(pb1**3*pb2**3/r)*(1.0/(2*X**3*sqrt(r)))*sce

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
