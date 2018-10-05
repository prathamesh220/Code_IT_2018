n=10
s=2.0
range(1,n+1)
F=[i**(-s) for i in range(1,11)]
P=[i/sum(F) for i in F ]

for i in range(n):
    print i+1, P[i]
    print '\n'
    
print ';'    