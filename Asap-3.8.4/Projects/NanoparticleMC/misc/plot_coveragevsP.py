import os
if not os.path.exists('AdsorptionParameters.py'):
	os.symlink('../AdsorptionParameters.py', 'AdsorptionParameters.py')
import matplotlib.pyplot as plt
import asap3.nanoparticle_mc.langmuirExpression as le
import numpy as np

#Test the CO Plot now

#Get one coverage for each CN for each pressure at 300 K

x = np.linspace(1,1E4,100)
covs = np.zeros((len(x),16))
for i in range(len(x)):
	covs[i] = le.getCoverages(T=300,P=x[i],species="AuCO")


CN_4 = [c[4] for c in covs]
CN_6 = [c[6] for c in covs]
CN_8 = [c[8] for c in covs]
CN_9 = [c[9] for c in covs]
CN_12 = [c[12] for c in covs]
plt.plot(x,CN_4,label="CN 4")
plt.plot(x,CN_6,label="CN 6")
plt.plot(x,CN_8,label="CN 8")
plt.plot(x,CN_9,label="CN 9")
plt.plot(x,CN_12,label="CN 12")
plt.legend(loc=3)
plt.ylim([0,1.1])
plt.xlabel("Pressure[Pa]",fontsize=20)
plt.ylabel("Coverage",fontsize=20)
plt.show()
