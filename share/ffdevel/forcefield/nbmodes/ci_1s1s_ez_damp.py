from sympy import *

r, pb1 = symbols('r pb1')

# DOI: 10.1002/jcc.20520
# eq 20a

ci   = exp(-pb1*r)*(1 + pb1*r/2)

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
