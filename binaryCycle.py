import sys
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import fluids
from scipy.optimize import minimize
import pandas as pd

class wellBore(object):
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
        fr=0
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
            if (p<101325 and self.proizvodna):
                self.info='izlaz na:' + str(zi) + ' m \n'
                ro_avg=(ro+CP.PropsSI('D', 'T', (273.15 + self.T), 'P', WHP, self.fluid))*0.5
                pump=zi*ro_avg*g+WHP
                break
        if (self.proizvodna==False):
            self.info = ''
            ro_avg = (ro + CP.PropsSI('D', 'T', (273.15 + self.T), 'P', WHP, self.fluid)) * 0.5
            pump = p - WHP


        self.info += 'potrebno pumpati: ' + str(pump) + ' Pa \n'
        self.info += 'ro_avg = ' + str(ro_avg) + ' kg/m3 \n'

        self.Ppump=self.q*pump*0.001    # snaga pumpe, kW
        self.info += 'Potrebna snaga pumpe = ' + str(self.Ppump) + ' kW \n'

class ORC(object):
    """
    assumptions and boundary conditions:
    - input heat source (wellbore) temperature and pressure are constant
    - all processes are steady-state and adiabatic.
    - friction and heat losses are neglected.
    - evaporator and condenser efficiency is 100 %
    - isentropic efficiency of pump and turbine are 0.75.
    """

    def __init__(self, fluid='Isopentane', geofluid='H2O', T_prod=175, T_inj=50, p_prod=10, T_out=35, mass_flow_geo=55.9636):
        self.fluid=fluid
        self.geofluid=geofluid
        self.m_gf=mass_flow_geo                 # mass flow rate of geofluid (kg/s)
        self.m_wf=10                            # mass flow rate of working fluid (kg/s)
        self.Tc = CP.PropsSI('Tcrit', fluid)    # critical temperature, K
        self.pc = CP.PropsSI('Pcrit', fluid)    # critical pressure, Pa
        self.p_prod=p_prod*1e5                  # pressure of geofluid, Pa
        self.T_inj=T_inj+273.15                 # (assumed) injection temperature of geofluid, K
        self.T_prod=T_prod+273.15               # produced geofluid temperature, K
        self.T_in=self.T_prod                   # = maximum temperature after superheating (evaporator), K
        self.T_out=T_out+273.15                 # lowest temperature, after cooling (K)
                                                # = maximum cooling fluid temperature
        self.info=''                            # info warning messages

        self.dT_pinch=5          # temperature difference between pinch point and workin fluid temperature, K
        self.dp=1.5e5            # pressure difference at the pump, Pa

    def pump(self, dp=None, T3=None):
        """
        isentropic pump (adiabatic and reversible)
        :param dp: desired pressure after pump
        :param T3: (calculated) temperature after condenser
        :return: sets states (h, s, T, p, dens, Q) at point 1 and 2 (before and after pump) in ORC system
        """

        if dp==None:
            dp=self.dp
        else:
            self.dp=dp

        if T3 == None: self.T3 = self.T_out
        self.T3c = self.T3                          # temperature at the inlet of condenser, K

        self.p3=CP.PropsSI('P', 'T', self.T3, 'Q', 0, self.fluid)

        self.p3c = CP.PropsSI('P', 'T', self.T3c, 'Q', 1, self.fluid)
        self.s3c = CP.PropsSI('S', 'T', self.T3c, 'Q', 1, self.fluid)
        self.h3c = CP.PropsSI('H', 'T', self.T3c, 'Q', 1, self.fluid)

        self.h3=CP.PropsSI('H', 'T', self.T3, 'Q', 0, self.fluid)       # J/kg
        self.s3=CP.PropsSI('S', 'T', self.T3, 'Q', 0, self.fluid)       # J/kg/K
        self.dens3 = CP.PropsSI('D', 'T', self.T3, 'Q', 0, self.fluid)

        self.p4=self.p3+dp
        self.s4 =self.s3
        self.T4 = CP.PropsSI('T', 'S', self.s4, 'P', self.p4, self.fluid)
        self.h4=CP.PropsSI('H', 'S', self.s4, 'P', self.p4, self.fluid)
        self.dens4 = CP.PropsSI('D', 'S', self.s4, 'P', self.p4, self.fluid)
        self.W_pump=(self.h4-self.h3)*self.m_wf
        if self.W_pump<0:
            self.W_pump=0

    def evaporator(self, dTsh=1):
        """
        calculates preheater, boiler and superheater, takung into account pinch point analysis
        :param dTsh: temperature difference between boiler and superheater, K
        :return:
        """
        self.p1=self.p4

        # boiler, x=1 (preheater)
        self.T1b = CP.PropsSI('T', 'P', self.p1, 'Q', 0, self.fluid)        # preheater=boiler temperature, K
        self.h1x = CP.PropsSI('H', 'P', self.p1, 'Q', 0, self.fluid)        # J/kg
        self.Qwfph=(self.h1x-self.h4)*self.m_wf                             # working fluid heat flux in preheater, J/s
        self.dens1x = CP.PropsSI('D', 'P', self.p1, 'Q', 0, self.fluid)
        self.s1x = CP.PropsSI('S', 'P', self.p1, 'Q', 0, self.fluid)        # J/kg/K

        # boiler, y=1 (boiler)
        self.h1y = CP.PropsSI('H', 'P', self.p1, 'Q', 1, self.fluid)
        self.dens1y = CP.PropsSI('D', 'P', self.p1, 'Q', 1, self.fluid)
        self.s1y = CP.PropsSI('S', 'P', self.p1, 'Q', 1, self.fluid)
        self.Qwfb=(self.h1y-self.h1x)*self.m_wf                             # working fluid heat flux in boiler, J/s

        # superheater, y=1
        self.T1=self.T1b+dTsh
        self.s1 = CP.PropsSI('S', 'P', self.p1, 'T', self.T1, self.fluid)
        if self.s1<self.s3c:        # check for expander to work outside two-phase region
            self.s1=self.s3c
            self.T1 = CP.PropsSI('T', 'P', self.p1, 'S', self.s1, self.fluid)
        self.h1=CP.PropsSI('H', 'P', self.p1, 'T', self.T1, self.fluid)
        self.Q1=(self.h1-self.h1y)*self.m_wf                                # wf heat flux in superheater, J/s

        self.Q_heating=self.Qwfph+self.Qwfb+self.Q1

    def expander(self):
        self.s2=self.s1
        self.p2=self.p3c
        self.T2=CP.PropsSI('T', 'P', self.p2, 'S', self.s2, self.fluid)
        self.h2=CP.PropsSI('H', 'P', self.p2, 'S', self.s2, self.fluid)
        self.W_expander = (self.h1-self.h2)*self.m_wf
        if self.W_expander<0:
            self.W_expander=0

    def condenser(self, coolingFluid='H2O', Tin_cooling=None):
        """

        :return:
        """
        # precooling
        pin_cooling=1.5*1e5               # inlet pressure of cooling fluid
        if Tin_cooling==None:
            Tin_cooling=12+273.15       # inlet temperature of cooling fluid
        hin_cooling=CP.PropsSI('H', 'P', pin_cooling, 'T', Tin_cooling, coolingFluid)
        Tout_cooling=self.T3c
        hout_cooling = CP.PropsSI('H', 'P', pin_cooling, 'T', Tout_cooling, coolingFluid)
        self.Q_cooler=(self.h2-self.h3c)*self.m_wf                # heat flux in cooler, J/s
        self.m_cf = self.Q_cooler/(hout_cooling-hin_cooling)      # mass flow rate of cooling fluid

        # condenser
        self.Q_condenser=(self.h3c-self.h3)*self.m_wf

    def pinchAnalysis(self):
        """
        sets geofluid enthalpies to the same relative enthalpy point and calculates temperatures at each
        point of interest (to be comparable with enthalpies and temperatures of working fluid)
        :return:
        """
        hgf0=CP.PropsSI('H', 'P', self.p_prod, 'T', self.T_inj, self.geofluid)
        Tgf = []
        dhgf=self.Qwf/self.m_gf                  # changes in entalpies for geofluid
        h_gf=np.array([hgf0, hgf0+dhgf[0], hgf0+dhgf[0]+dhgf[1], hgf0+dhgf[1]+dhgf[2]])
        for i, h in enumerate(h_gf):
            Tgf.append(CP.PropsSI('T', 'P', self.p_prod, 'H',h, self.geofluid))
        Tgf=np.array(Tgf)
        h_gf=h_gf-hgf0
        return (Tgf, h_gf)

    def plotORC(self, hgf):

        plt.plot(self.Qwf, self.Tgf[1:], marker='o', linestyle='-', color='r')
        plt.plot(self.Qwf, self.Twf[1:], marker='o', linestyle='-', color='b')
        plt.xlabel('Q, J/s')
        plt.ylabel('T, K')
        plt.grid(which='both')
        plt.show()

        plt.plot(hgf, self.Tgf, marker='o', linestyle='-', color='r')
        plt.plot(hgf, self.Twf, marker='o', linestyle='-', color='b')
        plt.xlabel('h, J/kg')
        plt.ylabel('T, K')
        plt.grid(which='both')
        plt.show()

        # # plot T-s
        # s=[self.s1, self.s2, self.s3c, self.s3, self.s4, self.s1x, self.s1y]
        # T=[self.T1, self.T2, self.T3c, self.T3, self.T4, self.T1b, self.T1b]
        # sQ=np.arange(self.s3-10, self.s2+10, 50).tolist()
        #
        # s_crit = (CP.PropsSI('S', 'P', self.pc, 'T', self.Tc, self.fluid))
        # TQ=[]
        # for s in sQ:
        #     if s<s_crit:
        #         TQ.append(CP.PropsSI('T', 'S', s, 'Q', 1, self.fluid))
        #     else:
        #         TQ.append(CP.PropsSI('T', 'S', s, 'Q', 0, self.fluid))
        #
        # plt.plot(sQ, TQ, linestyle='-', color='black')
        # plt.plot(s, T, marker='o', linestyle='-', color='b')
        # plt.xlabel('s, J/kg/K')
        # plt.ylabel('T, K')
        # plt.grid(which='both')
        # plt.show()

    def runORC(self, plot=False):
        self.pump(self.dp)
        self.evaporator()
        self.expander()
        self.condenser()
        pinched=False
        self.Qwf = np.array([self.Qwfph, self.Qwfb, self.Q1]).cumsum()
        self.Tgf, hgf = self.pinchAnalysis()
        self.Twf=[self.T4, self.T1b, self.T1b, self.T1]
        if plot: self.plotORC(hgf)
        while not pinched:
            if (self.Tgf[1]-self.dT_pinch)<self.T1b:
                self.dp-=0.25*1E5
                self.pump(dp=self.dp)
                self.evaporator()
                self.expander()
                self.condenser()
                #print('evaporator: pressure reduced by 0.25 bar (to %s)' %(self.dp/1e5))
            elif (self.Tgf[1]-self.dT_pinch-self.T1b)>5:
                self.dp+=0.25*1E5
                self.pump(dp=self.dp)
                self.evaporator()
                self.expander()
                self.condenser()
                #print('evaporator: pressure increased by 0.25 bar (to %s)' %(self.dp/1e5))
            else:
                pinched=True

        self.eta_ORC=(self.W_expander-self.W_pump)/self.Q_heating       # thermal efficiency

        self.Qwf = np.array([self.Qwfph, self.Qwfb, self.Q1]).cumsum()
        self.Tgf, hgf = self.pinchAnalysis()
        self.Twf=[self.T4, self.T1b, self.T1b, self.T1]
        if plot: self.plotORC(hgf)

    def summarizeORC(self):
        column_names=['stage', 'T (°C)', 'p (bar)', 's (J/kg/K)', 'h (J/kg)', 'Q (W)']
        stages = ['2 - expander', '3 - condenser y', '3 - condenser x', '4 - pump',
                  '1 - preheater x', '1 - boiler', '1 - superheater']
        t=np.array([self.T2, self.T3c, self.T3, self.T4, self.T1b, self.T1b, self.T1])-273.15
        p=np.array([self.p2, self.p3c, self.p3, self.p4, self.p4, self.p4, self.p1])/1e5
        s = np.array([self.s2, self.s3c, self.s3, self.s4, self.s1x, self.s1y, self.s1])
        h = np.array([self.h2, self.h3c, self.h3, self.h4, self.h1x, self.h1y, self.h1])
        Q = np.array([self.W_expander, self.Q_cooler, self.Q_condenser, self.W_pump, self.Qwfph, self.Qwfb, self.Q1])
        self.tableORC=pd.DataFrame(columns=column_names)
        self.tableORC['stage']=stages
        self.tableORC['T (°C)']=t
        self.tableORC['p (bar)']=p
        self.tableORC['s (J/kg/K)']=s
        self.tableORC['h (J/kg)']=h
        self.tableORC['Q (W)']=Q

        self.summary = pd.DataFrame(columns=['parameter', 'value', 'unit'])
        self.summary.loc['T_inj']=['geofluid injection temperature', (self.T_inj-273.15), 'C']
        self.summary.loc['T3c']=['cooling fluid temperature', (self.T3c - 273.15), 'C']
        self.summary.loc['m_wf'] = ['mass flow rate (working fluid)', self.m_wf, 'kg/s']
        self.summary.loc['m_gf'] = ['mass flow rate (geofluid)', self.m_gf, 'kg/s']
        self.summary.loc['m_cf'] = ['mass flow rate (cooling fluid))', self.m_cf, 'kg/s']
        self.summary.loc['eta_ORC'] = ['thermal efficiency', (self.W_expander - self.W_pump) / self.Q_heating, 'fraction']
        self.summary.loc['P_net'] = ['net power', (self.W_expander - self.W_pump)/1E6, 'MW']

        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 300)


    def optimizeORC(self):
        def objectiveORC(x):
            self.T_inj = x[0]
            self.T_out=x[1]
            self.m_wf=x[2]

            if self.T_inj<298.15:
                self.T_inj=298.15
            if self.T_inj>(self.T_prod-5):
                self.T_inj = (self.T_prod - 5)
            if self.T_out>self.T_inj:
                self.T_out=self.T_inj-3
            if self.T_out<288.15:
                self.T_out=288.15
            if self.m_wf<1:
                self.m_wf=1
            if self.m_wf>200:
                self.m_wf=200
            print (self.T_inj, self.T_out, self.m_wf)
            self.runORC()
            return(1/self.eta_ORC)

        x0=np.array([self.T_inj, self.T_out, self.m_wf])
        opt_r=minimize(objectiveORC,x0,
                       method='nelder-mead',
                       options={'xtol': 1e-8, 'disp': False})

        self.T_inj=opt_r.x[0]
        self.T_out=opt_r.x[1]
        self.m_wf=opt_r.x[2]





