import matplotlib.pyplot as plt
import numpy as np
import model






mu = 40.0
Tiarr = [0.025, 0.5, .75, 1.0, 1.25, 1.5, 1.75, 2.0]

r = np.linspace(0, 1800.0, 50000.0)

L = 149.867 / .004
d = .8770405
F = 22.035

w0 = 487.98634

#fig, ax = plt.subplots()
#for Ti in Tiarr:
#    m = model.forward3(r, L, d, F, Ti, 40.0,w0)
#    m /= m.max()
#    ax.plot(r, m, label="{0:4.3f}".format(Ti))
#ax.legend()
#ax.set_xlim(500,750)
#plt.show()


Varr = [0.0, 1.0, 2.0, 3.0, 4.0]
Ti = 0.6
fig, ax = plt.subplots()
for idx, V in enumerate(Varr):
    m = model.forward3(r, L, d, F, Ti, 40.0,w0, V=V, nlambda=512)
    m /= m.max()
    if idx == 0:
        i = np.where(np.logical_and(r>500.0, r<750))[0]
        print i[0]
        ii = np.argmax(m[i])
        iii = i[0]+ii
        ax.plot(r**2, m, label="{0}".format(V))
    else:
        i = np.where(np.logical_and(r>500.0, r<750))[0]
        ii = np.argmax(m[i])
        print r[ii+i[0]]
        #dr = r[ii+i[0]] - r[iii]
        dr2 = r[ii+i[0]]**2 - r[iii]**2
        ax.plot(r**2-dr2, m, label="{0}".format(V))

ax.legend()
plt.show()
