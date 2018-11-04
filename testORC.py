from binaryCycle import ORC as a
ORC=a()

ORC.runORC()
print (ORC.eta_ORC)
ORC.optimizeORC()
print (ORC.eta_ORC)

ORC.Tgf, hgf = ORC.pinchAnalysis()
ORC.plotORC(hgf)

print (ORC.m_cf)
ORC.summarizeORC()
print (ORC.tableORC)

import matplotlib.pyplot as plt
ORC.tableORC.plot(x='s (J/kg/K)', y='T (°C)')
plt.show()