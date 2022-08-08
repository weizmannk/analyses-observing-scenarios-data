
import numpy as np

import matplotlib.pyplot as plt

from astropy.io import ascii


 lc = ascii.read('ZTF21abbzjeq.csv', format='csv')
    
plt.figure()

for filter in ["g", "r"] :
        tab = lc[(lc['filters'] == filter)]
        
        tab.sort('time')
        idx = np.where(tab['mag'] > 50)[0]
        tab['mag'][idx] = np.nan


        t = tab['time']
        y = tab['mag']
        dy = tab["magerr"]
        
        if filter == 'g':
            plt.errorbar(t, y, dy, fmt='d', label='ZTF g', c='darkgreen')
        else:
            plt.errorbar(t, y, dy, fmt='o', label='ZTF r', c='darkred')

plt.gca().invert_yaxis()
plt.legend(loc='best')
plt.grid(True)
plt.xlabel(r'Time since first ZTF detection [days]')
plt.ylabel(r'absolute magnitude (AB)')
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.title('Kilonova follow-up ')
plt.show()