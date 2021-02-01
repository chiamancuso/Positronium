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
tDiff = ['t12', 't23', 't13']
coeff = [495.3662907422378, 499.20938900287547, 457.96722725099886]
data=[]
en = [258984.8434375515, 258388.2395968352, 238664.36333105486]
devEn = [115.89422582052933, 140.38598828118543, 112.05306994484995]


data = pd.read_csv('../Ortho_120_120.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
# print(data['E1'])
for n,j in enumerate(energy):
    data[j] = data[j]/coeff[n]


t12 = data['t1'] - data['t2']
t23 = data['t2'] - data['t3']
t13 = data['t1'] - data['t3']
dataT = {'t12':t12, 't23':t23, 't13':t13}
df = pd.DataFrame(data=dataT)

for n,j in enumerate(energy):
    # hist,bins = np.histogram(data[j], 1000)
    # m = np.where(hist == np.max(hist))
    data = data.loc[(data[j]<(2000)),:]
    data = data.loc[((data['t1']-data['t2'])<(0)),:]
    data = data.loc[((data['t1']-data['t2'])>(-70)),:]
    data = data.loc[((data['t2']-data['t3'])<(70)),:]
    data = data.loc[((data['t2']-data['t3'])>(0)),:]

for i,en in enumerate(energy):
    hist,bins = np.histogram(data[en], 60)
    binc = 0.5*(bins[:-1]+bins[1:])
    m = np.where(hist == np.max(hist))
    defin = (binc>(binc[m][0]-np.abs(1.1*binc[m][0])))*(binc<(binc[m][0]+np.abs(1.1*binc[m][0])))

    histCut = np.delete(hist[defin], np.where(hist[defin] == 0))
    bincCut = np.delete(binc[defin], np.where(hist[defin] == 0))
    err_hist=np.sqrt(histCut)


    par, par_var=curve_fit(gaussian_func, bincCut, histCut, p0=[np.max(histCut), binc[m][0], binc[m][0]*0.1],sigma=err_hist,absolute_sigma=True)
    sel_off = (bincCut>(par[1]-3*par[2]))*(bincCut<(par[1]+3*par[2]))
#    mean = np.mean(binc[sel_off])
    par_off, par_var_off = curve_fit(gaussian_func, bincCut[sel_off], histCut[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)
    print(par_off)
    plt.figure(figsize=(8,5))
    plt.scatter(bincCut,histCut,label='1st peak',color='lightblue')
    plt.errorbar(bincCut, histCut, err_hist,fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)

    plt.plot(np.linspace(par[1]-3*par[2],par[1]+3*par[2],1000), gaussian_func(np.linspace(par[1]-3*par[2],par[1]+3*par[2],1000), *par),label='peak fit',color='lightblue')

    plt.show()


# plt.figure(figsize=(8,5))
# plt.hist2d(data['E1'], data['E2'],bins=20, label='Title')
# plt.colorbar()
# plt.show()




#
# for n,j in enumerate(tDiff):
#     df=df.loc[(df[j]>(-35))&(df[j]<(35)),:]
# #print(df.info())
#
#
#
#     # histT,binsT = np.histogram(dataT[n], 1000)
#
#     # bincT = 0.5*(binsT[:-1]+binsT[1:])
#
# plt.figure(figsize=(8,5))
# plt.hist(data['E1'],bins=100,label='Title')
# plt.hist(data['E2'],bins=100,label='Title')
# plt.hist(data['E3'],bins=100,label='Title')
# plt.show()
