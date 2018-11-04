from binaryCycle import wellBore as b

fluid='H2O'
d = 0.219               # promjer busotine, m
q=9640.49               # protok, m3/dan
proizvodna = True       # flag za proizvodnu ili utisnu
T = 175.                # temperatura, C (175, 30)
p_i=80.                 # tlak na dnu busotine, bar (80, 300)
dz=50.                  # korak proracunavanja po dubinama
dubina=2500.            # dubina dna busotine (2500, 3975)


Wprod=b(fluid=fluid, q=q, d=d, T=T, p=p_i, h=dubina, prod=proizvodna)
Wprod.optWellPressures(20)      # postavlja se ciljni WHP 20 bar

Winj=b(fluid=fluid, q=q, d=d, T=30, p=300, h=3975, prod=False)
Winj.optWellPressures(20)        #

print ('zbirni podatci za proizvodnu busotinu, Wprod \n --------------------------------------------------------------')
print (Wprod.info)

print ('zbirni podatci za utisnu busotinu, Winj \n ---------------------------------------------------------------')
print (Winj.info)

'----------------- analiza optimalnog tlaka WHP -----------------------------------------------------------'
p, wp, wi, wtot =[], [], [], []
for WHP in range(5, 100, 5):
    Wprod.optWellPressures(WHP)
    Winj.optWellPressures(WHP)
    p.append(WHP)
    wp.append(Wprod.Ppump)
    wi.append(Winj.Ppump)
    wtot.append(Wprod.Ppump+Winj.Ppump)

import pandas as pd
import matplotlib.pyplot as plt

column_names=['WHP', 'Wp', 'Wi', 'total']
optWHP=pd.DataFrame(columns=column_names)
optWHP.WHP=p
optWHP.Wp=wp
optWHP.Wi=wi
optWHP.total=wtot

title='snaga pumpe na proizvodnoj (Wp), utisnoj (Wi) i ukupna snaga (total)'
ax1=optWHP.plot(x='WHP', y='Wp', title=title)
ax2=optWHP.plot(x='WHP', y='Wi', color='green', ax=ax1)
ax3=optWHP.plot(x='WHP', y='total', color='red', ax=ax1)
ax1.set_ylabel("potrebna snaga, kW")
ax1.set_xlabel("tlak na uscu busotina, bar")
plt.show()

