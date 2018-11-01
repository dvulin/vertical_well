from  binaryCycle import binaryCycle as b

fluid='H2O'
d = 0.219               # promjer busotine, m
q=9640.49               # protok, m3/dan
proizvodna = False      # flag za proizvodnu ili utisnu
T = 30.                 # temperatura, C (175, 30)
p_i=300.                # tlak na dnu busotine, bar (80, 300)
dz=50.                  # korak proracunavanja po dubinama
dubina=3975.            # dubina dna busotine (2500, 3975)

c=b(fluid=fluid, q=q, d=d, T=T, p=p_i, h=dubina, prod=proizvodna)
c.optWellPressures(1)
print (c.info)