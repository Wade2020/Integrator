import numpy as np
import matplotlib.pyplot as plt

with open('Data.txt','r') as f:
	data = f.read()
	
data = data.split('\n')
ep = [row.split(' ')[0] for row in data]
tab = [row.split(' ')[1] for row in data]
value = [row.split(' ')[2] for row in data]
val = [float(count) for count in value]
tab2 = [row.split(' ')[3] for row in data]
value2 = [row.split(' ')[4] for row in data]
val2 = [float(count) for count in value2]

plt.figure(1)
plt.subplot(211)
plt.hist(val, bins = np.arange(0,1.0333,.0333), alpha=0.8, histtype = 'step')
plt.hist(val2, bins = np.arange(0,1.0333,.0333), alpha=0.8, histtype = 'step')
#plt.xscale('log')
plt.yscale('log')
plt.title('Normalized Monte Carlo Cross-section')



plt.subplot(212)
plt.title('Cross-section Integration for omega 2 = 79,000 and another value')
plt.plot(ep,tab,'o')
plt.plot(ep,tab2,'o')

plt.show()

#plt.hist(val, bins = np.arange(25,99985,3332))
#plt.hist(val2, bins = np.arange(25,1580,51))