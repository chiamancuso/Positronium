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

time=['14','24']
energy = ['E1','E2']
ch = ['14','24']
data=[]

for k,el in enumerate(time):
    data.append(pd.read_csv('../Coinc'+el+'noDel.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4']))


for n,j in enumerate(energy):
    hist,bins = np.histogram(data[n][j], 2000)
    m = np.where(hist == np.max(hist))

    data[n]=data[n].loc[(data[n][j]>(bins[m][0]-0.1*bins[m][0]))&(data[n][j]<(bins[m][0]+0.1*bins[m][0])),:]

t14 = data[0]['t1'] - data[0]['t4']
t24 = data[1]['t2'] - data[1]['t4']

dataT = [t14, t24]
#
# plt.figure(figsize=(8,5))
# plt.hist(t14,bins=100,label='Title')
# plt.show()
    # plt.hist(data[n][j],bins=100,label='Title')

for n,j in enumerate(energy):
    histT,binsT = np.histogram(dataT[n], 1000)

    bincT = 0.5*(binsT[:-1]+binsT[1:])

    # plt.figure(figsize=(8,5))
    # plt.hist(dataT[n],bins=100,label='Title')
    # plt.show()
    # plt.hist(data[n][j],bins=100,label='Title')

    ma = np.where(histT == np.max(histT))
    # print(bincT[ma])
    defin = (bincT>(bincT[ma][0]-np.abs(1.1*bincT[ma][0])))*(bincT<(bincT[ma][0]+np.abs(1.1*bincT[ma][0])))
    #defin = (bincT>(bincT[ma][0]-0.3*np.abs(bincT[ma][0])))*(bincT<(bincT[ma][0]+0.3*np.abs(bincT[ma][0])))


    histTcut = np.delete(histT[defin], np.where(histT[defin] == 0))
    bincTcut = np.delete(bincT[defin], np.where(histT[defin] == 0))

    err_histT=np.sqrt(histTcut)

    par, par_var = curve_fit(gaussian_func, bincTcut, histTcut,p0=[np.max(histTcut),bincT[ma][0],bincT[ma][0]],sigma=err_histT,absolute_sigma=True)
        # print(par)
    sel_off = (bincT>(par[1]-3*np.abs(par[2])))*(bincT<(par[1]+3*np.abs(par[2])))
        #print(bincT[sel_off])
    par_new, par_var_new = curve_fit(gaussian_func, bincTcut, histTcut, p0=par, sigma=err_histT,absolute_sigma=True)

    mu_err = np.sqrt(np.diag(par_var_new)[1])
    sigma_err = np.sqrt(np.diag(par_var_new)[2])

    print(n+1, par_new[0],par_new[1], mu_err, par_new[2], sigma_err)

    plt.figure(10*n+1, figsize=(8,5))
    g1=plt.scatter(bincTcut,histTcut,label='1st peak',color='lightblue')
    plt.errorbar(bincTcut, histTcut, err_histT,fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)

    g2=plt.plot(np.linspace(bincT[ma][0]-np.abs(bincT[ma][0]),bincT[ma][0]+np.abs(bincT[ma][0]),1000), gaussian_func(np.linspace(bincT[ma][0]-np.abs(bincT[ma][0]),bincT[ma][0]+np.abs(bincT[ma][0]),1000), *par),label='peak fit',color='lightblue')

    plt.title('Coincidence {}'.format(j),fontsize=17)
    plt.xlabel('Time [s]',fontsize=17)
    plt.ylabel('Count of photons',fontsize=17)
    plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
    plt.savefig('Coincidence_'+ch[n]+'.png')
