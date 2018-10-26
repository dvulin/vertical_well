# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 01:28:47 2018

@author: d
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import fluids

# http://www.coolprop.org/coolprop/HighLevelAPI.html

"""------------------- input vrijednosti ------------------------------------------------------"""
fluid = 'H2O'           # radni fluid, H2O, CO2
# http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids

e = 0.0004              # hrapavost, m
d = 0.10                # promjer busotine, m
q = 9640.49             # protok, m3/dan
proizvodna = True       # flag za proizvodnu ili utisnu
T = 175                 # temperatura, C (175, 30)
p_i=200                  # tlak na dnu busotine, bar (80, 300)
dz=100                   # korak proracunavanja po dubinama
dubina=2500             # dubina dna busotine (2500, 3975)


eD=e/d                                                              # relativna hrapavost
z = np.linspace(dubina, 0, int(dubina / dz + 1))                    # dubine, m
ro_sc = CP.PropsSI('D','T',(273.15+15.6),'P',101325,fluid)          # gustoca pri standardnim uvjetima, kg/m3


def funcRe(ro, v, d, mi):
    """
    ro - gustoca [kg/m3]
    v  - srednja brzina [m/s]
    d  - unutrasnji promjer cijevi [m]
    mi - dinamicka viskoznost fluida [Pas]
    """
    return (ro * v * d / mi)


def funcReg(Re, e):
    if Re < (3600 / e):
        if Re < 2000: fr = "strujno"
        if Re < 100: fr = "laminarno jak utjecaj Re"
        if Re < 1: fr = "puzece"
        if Re >= 2000 and Re < 4000: fr = "kriticno"
        if Re >= 4000: fr = "tranzicijsko/turbulentno"
    else:
        fr = "turbulentno (potpuno)"
    return (fr)


def funcf(Re, e, d):
    # racuna (Darcyev) friction factor (Haaland, 1983)
    return ((-1.8 * np.log10(6.9 / Re + (e / (3.7 * d) ** 1.11))) ** (-2))


def funcdpf(f, dz, d, ro, v):
    # Darcy–Weisbachova jednadzba gubitaka radi trenja
    d_p = f * (dz / d) * ro * (v ** 2) * 0.5            # gubitci radi trenja
    return (d_p)


g = 9.8066  # akceleracija zbog gravitacije m2/s
q = q / 86400.  # m3/s

p = p_i
pg = []
Reg = []
flow_type = []
dens=[]
mu=[]

for zi in z:
    pg.append(p)
    ro = CP.PropsSI('D','T',(273.15+T),'P',(p*1e5),fluid)   # gustoca pri trenutnom tlaku
    mi = CP.PropsSI('V','T',(273.15+T),'P',(p*1e5),fluid)   # viskoznost pri trenutnom tlaku
    dens.append(ro)
    mu.append(mi)
    print(p, ro, mi)
    q_rc = q * ro_sc / ro               # volumni protok u busotini pri zadanoj gustoci i protoku na povrsini
    v = q_rc / (np.pi * 0.25 * d ** 2)                    # brzina protjecanja m/s
    Re = funcRe(ro, v, d, mi)                             # Reynoldsov broj
    Reg.append(Re)
    flow_type.append(funcReg(Re, e))

    # Clamond, Didier, 2009. “Efficient Resolution of the Colebrook Equation.”
    f = fluids.friction_factor(Re=Re, eD=eD)                # Clamondova jednadzba - bolja od Haalandove
    dpf = funcdpf(f, dz, d, ro, v)

    if proizvodna:
        dp = ro * g * dz - ro * v * dz + dpf                # - ro * v * dz ?
    else:
        dp = ro * g * dz - ro * v * dz - dpf

    dp = dp * 1e-5
    p = p - dp
    if p < 0: break

THP='{0:.4g}'.format(pg[-1])

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title('Proizvodna busotina: %s, p_BHP= %s bar \n broj koraka: %s \n q= %s m3/s'
             % (proizvodna, p_i, len(pg), '{0:.2g}'.format(q)))
ax.set_xlabel('dubina, m')
ax.set_ylabel('tlak, bar')
ax.legend(loc='best', title=fluid)
plt.ylim(1, p_i + 20)
ax.annotate((str(THP)+ ' bar'), xy=(z[-1]+500, (int(pg[-1] -5 ))))
plt.plot(z[:len(pg)], pg, linestyle='-', marker='o')
plt.show()

print("rjesenje za tlak na uscu: p_wh = %s bar \n" % THP)


"""------------------------------ plotaj gustoce i viskoznosti -------------------------------"""

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(221)
ax.set_title('viskoznost (T=%s C)' % (T))
ax.set_xlabel('p, bar')
ax.set_ylabel('viskoznost, Pa s')
ax.legend(loc='right', title=fluid)
ax.plot(pg, mu, linestyle='-', marker='o')

ax2 = fig.add_subplot(222)
ax2.set_title('gustoca (T=%s C)' % (T))
ax2.set_xlabel('p, bar')
ax2.set_ylabel('gustoca, kg/m^3')
ax2.legend(loc='best', title=fluid)
ax2.plot(pg, dens, linestyle='-', marker='o')
plt.show()




