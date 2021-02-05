import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv
import pandas
plt.rcParams['axes.facecolor'] = 'white'

def Cal(E,a):
    i=a*E
    return(i)

def Calibration(ADC,err_ADC,E_exp,Channel_number):  ##Channel number is a word; ex: "channel 2"
    """Calibration fit C=f(E)"""

    par_Cal, par_Cal_var=curve_fit(Cal,E_exp,ADC,sigma=err_ADC,absolute_sigma=True)

    a=round(par_Cal[0],3)                   ### a coef to go like C=a*E
    err_a=round(np.sqrt(np.diag(par_Cal_var)[0]),3)


    #plt.figure(Channel_number)
    plt.errorbar(E_exp,ADC,yerr=err_ADC,fmt='b.',capsize=4)
    plt.plot(E_exp,Cal(E_exp,par_Cal[0]),c='k',label="a=%s"%a+"$\pm$%s"%err_a)
    plt.title('Calibration detector{} ADC=f(E)'.format(Channel_number),fontsize=20)
    plt.xlabel('E (keV)')
    plt.ylabel('ADC')
    # plt.annotate(r'$\a$',xy=(1000,np.mean(ADC)), xycoords='data',xytext=(+10, +30), textcoords='offset points', fontsize=16)
    # plt.annotate(f'={round(a,3)}+-{round(err_a,3)}',xy=(1200,np.mean(ADC)), xycoords='data',xytext=(+10, +30), textcoords='offset points', fontsize=16)
    plt.legend()
    plt.savefig('Calibration_channel_{}'.format(Channel_number))
    plt.show()

    return ([a,err_a])


def gaussian_func(x, A, mu, sigma, a, b):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))+a*x+b

def _2_gaussian_func(x, A1, mu1, sigma1, A2, mu2, sigma2, a, b):
    return A1 * np.exp( - (x - mu1)**2 / (2 * sigma1**2))+a*x+b+A2 * np.exp( - (x - mu2)**2 / (2 * sigma2**2))


def chi2red(sel,are1or2gaussians,*p):
    if are1or2gaussians==1:
        chi2=sum(((hist[sel] - gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-5
    else:
        chi2=sum(((hist[sel] - _2_gaussian_func(binc[sel],*p) )/ err_hist[sel]) ** 2)
        dof=len(binc[sel])-8
    #print(chi2, dof)
    return round(chi2/dof,2)

Elements=['Co','Na']#,'Cs']
Channel=['ch4']
Title =['Channel 4: TERESA']

for k,el in enumerate(Elements):
    data= pd.read_csv('../Calib'+el+'4new2.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])

    ##
    channel_4=data.loc[(data.iloc[:,6]>=1),:]
#     channel_2=data.loc[(data.iloc[:,2]>=1),:]
#     channel_3=data.loc[(data.iloc[:,4]>=1),:]

    plt.figure(10*k,figsize=(16,10))
    plt.hist(channel_4['E4'],bins=400,label=Title[0])
#     plt.hist(channel_3['E3'],bins=1000,label=Title[1])
#     plt.hist(channel_2['E2'],bins=1000,label=Title[2])
    plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
    plt.title('Histograms for channel 4'+el,fontsize=17)
    plt.xlabel('ADC',fontsize=17)
    plt.ylabel('Count of photons',fontsize=17)
    for tickLabel in plt.gca().get_xticklabels()+plt.gca().get_yticklabels():
        tickLabel.set_fontsize(15)
    #plt.show()
    plt.savefig('Spectrach4_'+el)

    E4 = np.array(channel_4['E4'])
    E4 = np.delete(E4, np.where(E4 == 0))
#     E2 = np.array(channel_2['E2'])
#     E3 = np.array(channel_3['E3'])

    E = [E4]
    print(E)
    parameters1 = [['Channel', 'Amplitude', 'Average', 'AverageError', 'Sigma', 'SigmaError', 'a', 'b','integral']]
    paramCo1 = [['Channel', 'Amplitude', 'Average', 'AverageError', 'Sigma', 'SigmaError', 'a', 'b']]
    #parameters2 = [['Channel', 'Amplitude', 'Average', 'Sigma', 'a', 'b']]
    parameters2 = []
    paramCo2=[]

    for i,ch in enumerate(Channel):

        hist,bins = np.histogram(E[i], 400)
        binc = 0.5*(bins[:-1]+bins[1:])
        err_hist=np.sqrt(hist)
        conr = (binc>2*10**5)
        # plt.figure(figsize=(8,5))
        # plt.scatter(binc,hist,label='try',color='lightblue')
        # #plt.hist(E[i],bins=100,label='here')
        # plt.show()
        m = np.where(hist == np.max(hist))
        if (el == 'Co'):
            sel_2max = (binc>(binc[m]*1.07))
            maxim = np.where(hist == np.max(hist[sel_2max]))
            defin = (binc>(binc[m]*0.93))*(binc<binc[maxim]*1.07)
            mean_2max = np.mean(binc[sel_2max])
            #print(binc[m], binc[maxim])
            par_2max, par_va_2max = curve_fit(_2_gaussian_func, binc[defin], hist[defin], p0=[np.max(hist), binc[m][0], binc[m][0]*0.01, np.max(hist), binc[maxim][0], binc[maxim][0]*0.01, 1, 1 ],sigma=err_hist[defin],absolute_sigma=True)
            #print(*par_2max)
            sel_off_2max = (binc>(par_2max[1]-3*par_2max[2]))*(binc<(par_2max[4]+3*par_2max[5]))
            par_2max_off, par_va_2max_off = curve_fit(_2_gaussian_func, binc[sel_off_2max], hist[sel_off_2max], p0=par_2max, sigma=err_hist[sel_off_2max],absolute_sigma=True)

            chi2=chi2red(defin,2,*par_2max)
            a=chi2**0.5
         #   print(a)

            par_2max_off, par_va_2max_off = curve_fit(_2_gaussian_func, binc[sel_off_2max], hist[sel_off_2max], p0=par_2max, sigma=a*err_hist[sel_off_2max],absolute_sigma=True)

            mu_err = np.sqrt(np.diag(par_va_2max_off)[1])
            sigma_err = np.sqrt(np.diag(par_va_2max_off)[2])

            mu_err_2 = np.sqrt(np.diag(par_va_2max_off)[4])
            sigma_err_2 = np.sqrt(np.diag(par_va_2max_off)[5])

         #   print(i+1, par_2max_off[0],par_2max_off[1], mu_err, par_2max_off[2], sigma_err, par_2max_off[6],par_2max_off[7])

            paramCo1.append([i+1, par_2max_off[0],par_2max_off[1], mu_err, par_2max_off[2], sigma_err, par_2max_off[6],par_2max_off[7]])
            paramCo2.append([i+1, par_2max_off[3],par_2max_off[4], mu_err_2, par_2max_off[5], sigma_err_2, par_2max_off[6],par_2max_off[7]])


            plt.figure(10*k+i+1, figsize=(16,10))
            g1=plt.scatter(binc[defin],hist[defin],label='1st peak',color='lightblue')
            plt.errorbar(binc[defin], hist[defin], err_hist[defin],fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)
            g2=plt.plot(binc[defin], _2_gaussian_func(binc[defin], *par_2max),label='peak fit'+' ($\chi^2$/dof = %s'%chi2+')',color='lightblue')

        else:
            sel = (binc>(binc[m]*0.93))*(binc<(binc[m]*1.07))
            mean = np.mean(binc[sel])

            par, par_var=curve_fit(gaussian_func, binc[sel], hist[sel], p0=[np.max(hist), mean, mean*0.01, 1, 1],sigma=err_hist[sel],absolute_sigma=True)
            sel_off = (binc>(par[1]-3*par[2]))*(binc<(par[1]+3*par[2]))
        #    mean = np.mean(binc[sel_off])
            par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=err_hist[sel_off],absolute_sigma=True)

            chi2=chi2red(sel_off,1,*par_off)
            a=chi2**0.5
#            print(a)

            par_off_, par_var_off = curve_fit(gaussian_func, binc[sel_off], hist[sel_off], p0=par, sigma=a*err_hist[sel_off],absolute_sigma=True)

            mu_err = np.sqrt(np.diag(par_var_off)[1])
            sigma_err = np.sqrt(np.diag(par_var_off)[2])
            parameters1.append([i+1, par_off[0],par_off[1], mu_err, par_off[2], sigma_err, par_off[3],par_off[4],np.sqrt(2*np.pi)*par_off[0]*par_off[2]])

            plt.figure(10*k+i+1, figsize=(16,10))
            g1=plt.scatter(binc[sel_off],hist[sel_off],label='1st peak',color='lightblue')
            plt.errorbar(binc[sel_off], hist[sel_off], err_hist[sel_off],fmt='.', color='lightblue',ecolor='lightgray', elinewidth=2, capsize=0)
            g2=plt.plot(binc[sel_off], gaussian_func(binc[sel_off], *par_off),label='1st peak fit'+' ($\chi^2$/dof = %s'%chi2+')',color='lightblue')

            if (el=='tatat'): # or 'Na'

                sel_off2 = (binc>(par[1]+3*par_off[2]))
                ma = np.where(hist == np.max(hist[sel_off2]))
                c=binc[ma]
                d=np.max(c)
                # print(binc[sel_off2])
                # print(c)
                # print(type(c))
                sel_new = (binc>(d*0.90))*(binc<(d*1.10))
                mean_new = np.mean(binc[sel_new])


                par_new, par_var_new=curve_fit(gaussian_func, binc[sel_new], hist[sel_new], p0=[np.max(hist[sel_off2]), mean_new, mean_new*0.01, 1, 1],sigma=err_hist[sel_new],absolute_sigma=True)
                sel_off_new = (binc>(par_new[1]-3*par_new[2]))*(binc<(par_new[1]+3*par_new[2]))

                par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off_new], hist[sel_off_new], p0=par_new, sigma=err_hist[sel_off_new],absolute_sigma=True)

                chi2=chi2red(sel_off_new,1,*par_off)
                a=chi2**0.5
                #print(a)
                par_off, par_var_off = curve_fit(gaussian_func, binc[sel_off_new], hist[sel_off_new], p0=par_new, sigma=a*err_hist[sel_off_new],absolute_sigma=True)

                mu_err = np.sqrt(np.diag(par_var_off)[1])
                sigma_err = np.sqrt(np.diag(par_var_off)[2])
                parameters2.append([i+1, par_off[0], par_off[1], mu_err, par_off[2], sigma_err, par_off[3],par_off[4],np.sqrt(2*np.pi)*par_off[0]*par_off[2]])


                g3=plt.scatter(binc[sel_off_new],hist[sel_off_new],label='2nd peak',color='orange')
                #chi2a=round(chi2/a**2,2)
                plt.errorbar(binc[sel_off_new], hist[sel_off_new], err_hist[sel_off_new],fmt='.',color='orange',ecolor='lightgray', elinewidth=2, capsize=0)
                g4=plt.plot(binc[sel_off_new], gaussian_func(binc[sel_off_new], *par_off),label='2nd peak fit'+' ($\chi^2$/dof = %s'%chi2+')',color='orange')

        plt.title(Title[i]+' '+el,fontsize=17)
        plt.xlabel('ADC',fontsize=17)
        plt.ylabel('Count of photons',fontsize=17)
        for tickLabel in plt.gca().get_xticklabels()+plt.gca().get_yticklabels():
            tickLabel.set_fontsize(15)
        plt.legend(bbox_to_anchor=(0.987, 0.982), loc=1, borderaxespad=0.)
        plt.savefig('FitPeak2_channel_'+el+ch.format(i))
        #plt.show()
    paramTot = parameters1 + parameters2
    if (el == 'Co'):
        paramTotCo = paramCo1 + paramCo2
        np.savetxt("ParametersCo1.csv", paramTotCo, delimiter="\t", fmt='%s')
    #print(paramTotCo, paramCo1, paramCo2)

    np.savetxt("Parameters"+el+".csv", paramTot, delimiter="\t", fmt='%s')



arr = [497387.27955199406,563055.9877577076,537355.6533663797]
arr2 = [18378.91000969369, 19719.546785745424, 20018.95363230451]
E_exp=np.array([1173.,1332.,1275.]) ##expected energy of the pics from the spectrums Co Na Cs in keV

print(Calibration(arr,arr2,E_exp,4))

    #Calibration(ADC,err_ADC,E_exp,Channel_number):
