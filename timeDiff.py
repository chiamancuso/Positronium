import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas

def gaussian_func(x, A, mu, sigma):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))

def chi2red(sel,are1or2gaussians,*p):
    if are1or2gaussians==1:
        chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-5
    else:
        chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-8
    #print(chi2, dof)
    return round(chi2/dof,2)


ch = ['12','23','13']
Title =['Channel 1: BERNA', 'Channel 2: DEMO', 'Channel 3: BLACK']
energy = ['E1','E2','E3']
time = ['t1','t2','t3']
coeff = [495.3662907422378, 499.20938900287547, 457.96722725099886]
data=[]
en = [258984.8434375515, 258388.2395968352, 238664.36333105486]
devEn = [115.89422582052933, 140.38598828118543, 112.05306994484995]


data = pd.read_csv('../Ortho_120_120.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
print(data['E1'])
for n,j in enumerate(energy):
    data[j] = data[j]/coeff[n]


t12 = data['t1'] - data['t2']
t23 = data['t2'] - data['t3']
t13 = data['t1'] - data['t3']
dataT = [t12, t23, t13]

for n,j in enumerate(energy):
    # hist,bins = np.histogram(data[j], 1000)
    # m = np.where(hist == np.max(hist))
    #
    data = data.loc[(data[j]<(2000)),:]

plt.figure(figsize=(8,5))
plt.hist2d(data['E1'], data['E2'],bins=20, label='Title')
plt.colorbar()
plt.show()





for n,j in enumerate(energy):
    sel = (dataT[n]>(-30))&(data[n]<(30))
    #dataT[n][sel]
    histT,binsT = np.histogram(dataT[n], 1000)

    bincT = 0.5*(binsT[:-1]+binsT[1:])

    plt.figure(figsize=(8,5))
    plt.hist(dataT[n],bins=100,label='Title')
    plt.show()
