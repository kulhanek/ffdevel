from sympy import *

r, pb1, pb2 = symbols('r pb1 pb2')

pr1     = pb1*r
pr2     = pb2*r
pfe1    = exp(-pr1)
pfe2    = exp(-pr2)
a2      = pb1 ** 2
b2      = pb2 ** 2
inva2b2 = 1.0/(b2 - a2)
ci      =  + pfe1 * (b2*inva2b2)**2 * (1.0 - 2.0*a2*inva2b2 + 0.5*pr1) + pfe2 * (a2*inva2b2)**2 * (1.0 + 2.0*b2*inva2b2 + 0.5*pr2)

sci     = ci
sdci    = diff(ci,r)

[replacements,  reduced_exprs] = cse( (sci,sdci) )

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
print("        {a:3} = {b:}".format(a='ci', b=str(reduced_exprs[0])))
print("        {a:3} = {b:}".format(a='dci',b=str(reduced_exprs[1])))
