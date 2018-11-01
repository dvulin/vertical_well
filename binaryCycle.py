import sys
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import fluids
import scipy.optimize as sci
import pandas as pd

class binaryCycle(object):
    def __init__(self, fluid='H2O', q=9640.49, d=0.219, T=175., p=80., h=2500., prod=True):
        """
        :param fluid: naziv fluida (prema definicijama u CoolPropu
        :param q: volumni protok fluida, m3/dan
        :param d: promjer busotine, m
        :param T: temperatura, C (175, 30)
        :param p: tlak na dnu busotine, bar
        :param h: dubina dna busotine, m
        :param prod: flag za proizvodnu ili utisnu (True/False)
        """
        self.e = 0.00004             # hrapavost, m
        self.d = d                   # promjer busotine, m
        self.q = q/86400.            # protok, m3/s
        self.proizvodna = prod       # flag za proizvodnu ili utisnu
        self.T = T                   # temperatura, C (175, 30)
        self.p_i=p                   # tlak na dnu busotine, bar (80, 300)
        self.Tc=CP.PropsSI('Tcrit', fluid)      # kriticna temperatura zadanog fluida, K
        self.dubina=h                # dubina dna busotine (2500, 3975)
        self.fluid=fluid             # tip fluida (CO2, H2O)
        self.region=None             # podrucje protjecanja
        self.ff=None                 # faktor trenja (friction factor) za zadane uvjete
        self.Ppump=0                 # potrebna snaga pumpe, kW
        self.Pcomp=0                 # potrebna snaga kompresora, kW
        self.info='No info'          # informacije za ispis

        self.ro=[]          # prazna lista gustoca koje se popne u kalkulaciji protoka u busotini
        self.mi=[]          # prazna lista viskoznosti
        self.Re_i=[]        # prazna lista Reynoldsovih brojeva
        self.fr=[]          # prazna lista faktora trenja
        self.flowType=[]    # zapis tipova protjecanja
        self.p=[]           # tlak na dubini z

        if prod:
            self.dz = 50  # korak proracunavanja po dubinama
        else:
            self.dz = 50

        self.eD = self.e / self.d  # relativna hrapavost
        self.z = np.linspace(self.dubina, 0, int(self.dubina / self.dz + 1))    #  dubine, m
        self.ro_sc = CP.PropsSI('D', 'T', (273.15 + 15.6), 'P', 101325, self.fluid)  # gustoca pri standardnim uvjetima, kg/m3

    def funcRe(self, ro, v, d, mi):
        """
        Racuna Reynoldsov broj
        ro - gustoca [kg/m3]
        v  - srednja brzina [m/s]
        d  - unutrasnji promjer cijevi [m]
        mi - dinamicka viskoznost fluida [Pas]
        :return: Reynoldsov broj
        """
        self.Re=ro * v * d / mi
        return (ro * v * d / mi)

    def funcReg(self, Re, e=None):
        """
        Racuna podrucje protjecanja fluida
        :param e: hrapavost, m
        :return: podrucje protjecanja (broj, opis)
        """
        region=['strujno', 'laminarno jak utjecaj Re',
                 'puzece', 'kriticno',
                 'tranzicijsko/turbulentno', 'turbulentno (potpuno)']
        if e==None: e = self.e
        if Re < (3600 / e):
            if Re < 2000: fr = 0
            if Re < 100: fr = 1
            if Re < 1: fr = 2
            if Re >= 2000 and Re < 4000: fr = 3
            if Re >= 4000: fr = 4
        else:
            fr = 5
        self.region=fr
        return (fr, region[fr])

    def funcf(self, Re, e, d):
        # racuna (Darcyev) friction factor (Haaland, 1983)
        self.ff=(-1.8 * np.log10(6.9 / Re + (e / (3.7 * d) ** 1.11))) ** (-2)
        return (self.ff)

    def funcdpf(self, f, dz, d, ro, v):
        """
        vraca pad tlaka, tj. gubitke radi trenja za zadani interval dz i uvjete protjecanja
        :param f:   Darcyev friction factor
        :param dz:  interval, m
        :param d:   promjer cijevi, m
        :param ro:  gustoca fluida (na zadanom intervalu)
        :param v:
        :return:
        """
        # Darcy–Weisbachova jednadzba gubitaka radi trenja
        d_p = f * (dz / d) * ro * (v ** 2) * 0.5            # gubitci radi trenja
        return (d_p)

    def optWellPressures(self, WHP=1, phaseCheck=False):
        """
        :param WHP : trazeni tlak na uscu
        :return: tlak fluida na uscu (bar) i snaga pumpe potrebna za podizanje tlaka na tu razinu (W)
        """
        if WHP<1:
            WHP=1e5
        else:
            WHP=WHP*1e5
        g = 9.8066              # akceleracija zbog gravitacije m2/s
        p = self.p_i*1e5        # pocetni tlak (dno busotine, Pa)
        if not self.proizvodna:
            p=p-WHP             # oduzeti tlak na wellheadu ukoliko postoji

        e, eD=self.e, self.eD
        pg, Reg, flow_type, dens, mu, dpf, ro = [], [], [], [], [], 0, 0
        pump=0
        for zi in self.z:
            self.p.append(p)
            ro = CP.PropsSI('D', 'T', (273.15 + self.T), 'P', p , self.fluid)         # gustoca pri trenutnom tlaku
            mi = CP.PropsSI('V', 'T', (273.15 + self.T), 'P', p , self.fluid)         # viskoznost pri trenutnom tlaku
            self.ro.append(ro)
            self.mi.append(mi)
            q_rc = self.q * self.ro_sc / ro      # volumni protok u busotini pri zadanoj gustoci i protoku na povrsini
            v = q_rc / (np.pi * 0.25 * self.d ** 2)          # brzina protjecanja m/s
            Re = self.funcRe(ro, v, self.d, mi)              # Reynoldsov broj
            self.Re_i.append(Re)
            self.flowType.append(self.funcReg(Re, e))
            f = fluids.friction_factor(Re=Re, eD=eD)         # Clamondova jednadzba - bolja od Haalandove
                                        # Clamond, Didier, 2009. “Efficient Resolution of the Colebrook Equation.”
            self.fr.append(f)
            dpf = self.funcdpf(f, self.dz, self.d, ro, v)
            if self.proizvodna:
                dp = ro * g * self.dz - ro * v * self.dz + dpf
            else:
                dp = ro * g * self.dz + ro * v * self.dz - dpf
            p = p - dp
            if p<101325:
                self.info='izlaz na:' + str(zi) + ' m \n'
                ro_avg=(ro+CP.PropsSI('D', 'T', (273.15 + self.T), 'P', WHP, self.fluid))*0.5
                if self.proizvodna:
                    pump=zi*ro_avg*g+WHP
                else:
                    pump = zi * ro_avg * g - WHP
                self.info += 'potrebno pumpati: ' + str(pump) + ' Pa \n'
                self.info += 'ro_avg = ' + str(ro_avg) + ' kg/m3 \n'
                break

        self.Ppump=self.q*pump*0.001    # snaga pumpe, kW
        self.info += 'Potrebna snaga pumpe = ' + str(self.Ppump) + ' kW \n'