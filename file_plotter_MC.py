import numpy as np
import matplotlib.pyplot as plt

with open('MC.txt','r') as f:
	mc = f.read()
	
mc = mc.split('\n')
ep = [row.split(' ')[0] for row in mc]
tab = [row.split(' ')[1] for row in mc]
value = [row.split(' ')[2] for row in mc]
val = [float(count) for count in value]

plt.figure(1)
plt.subplot(211)
#plt.xscale('log')
plt.yscale('log')
plt.title('Monte Carlo CS Simulation')
plt.hist(val, bins = np.arange(0,1.0333,.0333), histtype = 'step')

plt.subplot(212)
plt.title('Cross-section Integration for omega 2 = 79,000 and another value')
plt.plot(ep,tab,'o')

plt.show()