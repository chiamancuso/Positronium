import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = False
from prompt_toolkit import HTML
from scipy.odr import ODR, Model, Data, RealData
import numpy as np
from scipy.stats import norm
import pandas as pd
from scipy.optimize import curve_fit
from terminaltables import AsciiTable
import csv


def simplegaussianl(par, x):
    return par[0] * np.exp( - (x - par[1])**2 / (2 * par[2]**2))+par[3]

def simplegaussian(x, A, mu, sigma, B):
    return A * np.exp( - (x - mu)**2 / (2 * sigma**2))+B



Ld=[138,142,146,150,154,158,162,166,170,172,174,176,178,180,182,184,186,188,190,194,198,202,206,210]
#Ldstr=['138','142','146','150','154','158','162','166','170','172','174','176','178','180','182','184','186','188','190','194','198','202','206','210']

res=[]

for k, angle in enumerate(Ld):
    data= pd.read_csv(str(angle)+'deg.txt',delimiter="\t",names = ['t1','E1','t2','E2','t3','E3','t4','E4'])
    channel_1=data.loc[(data.iloc[:,0]>=1),:]
    channel_2=data.loc[(data.iloc[:,2]>=1),:]
    E1 = np.array(channel_1['E1'])
    E2 = np.array(channel_2['E2'])
    #print (str(angle)+'deg.txt')
    hist,bins = np.histogram(E1, 1000)
    binc = (0.5*(bins[:-1]+bins[1:]))
#    binc = [int(x) for x in binc]
#     res.append(np.max(hist))
    #plt.figure(k)
  #  plt.plot(binc,hist)
    
    conr = (binc>2*10**5)
    m = np.where(hist == np.max(hist[conr]))
#     print(binc[m][0],k)
    sel = (binc>(binc[m][0]*0.93))*(binc<(binc[m][0]*1.07))
    res.append(sum(hist[sel]))

# print(res)

err_res=np.sqrt(res)
errx=1                                            #CHANGE ERRORS ON X
err_resx=len(Ld)*[errx]

Lc=np.linspace(140,220,200)

def graph(a,toc): #a:odr with xerrors or curvefit without? toc: number of the figure
    plt.figure(toc,figsize=(15,10))
    for tickLabel in plt.gca().get_xticklabels()+plt.gca().get_yticklabels():
        tickLabel.set_fontsize(15)
    plt.title("Number of coincidence counts in the peak of $Na^{22}$ wrt angle",fontsize=17)
    plt.xlabel(r'Angle $\theta$ (°)',fontsize=17)
    plt.ylabel('# coincidence events',fontsize=17)
    plt.scatter(Ld,res,label='\n Data with parameters'+\
               ' \n $16$ $cm$ source-detectors'+\
               ' \n $150$ $s/points$')
    if a==False:
        par, par_va = curve_fit(simplegaussian, Ld, res, p0=[70000, 180, 10,1000],sigma=err_res,absolute_sigma=True)
        chi2=round(sum(((res - simplegaussian(Ld,*par) )/ err_res) ** 2)/(len(Ld)-3),2)
        plt.plot(Lc,simplegaussian(Lc, *par),color='gold',label='Fit with '+r'$A\exp\{\frac{-(\theta-\mu)^2}{2\sigma^2}\}+Cst$'+\
                  ' \n $A =$%s'%int(par[0])+' $ \pm$ %s'%int(np.sqrt(np.diag(par_va)[0]))+' #'+\
                  ' \n $\mu =$ %s'%round(par[1],1)+' $\pm$ %s'%round(np.sqrt(np.diag(par_va)[1]),1)+'°'+\
                  ' \n $\sigma =$ %s'%round(par[2],1)+ '$\pm$ %s'%round(np.sqrt(np.diag(par_va)[2]),1)+'°'+\
                  ' \n $Cst=$ %s'%int(par[3])+' $\pm$ %s'%int(np.sqrt(np.diag(par_va)[3]))+' #'+\
                  ' \n $\chi^2/dof = $ %s'%chi2)
        plt.errorbar(Ld, res, err_res,fmt='.',label=r'$y=\sqrt{counts}$ '+\
                     '\n $x=0°$', color='black',ecolor='lightgray', elinewidth=3, capsize=0)

    else:
        data = RealData(Ld,res,err_resx,err_res)
        model = Model(simplegaussianl)

        odr = ODR(data, model, [73021, 183, 11,1208])
        odr.set_job(fit_type=2)
        output = odr.run()
        
        xn = Lc
        yn = simplegaussianl(output.beta, xn)
        
        #pl.hold(True)
        #plot(Ld,res,'ro')
        #print(x,y)
        plt.plot(xn,yn,'k-',label=' ODR leastsq fit'+\
            ' \n $\chi^2/dof = $ %s'%round(output.sum_square/(len(Ld)-3),2)+'\n')
        
        odr.set_job(fit_type=0)
        output = odr.run()
        par,par_va=output.beta,output.cov_beta
        yn = simplegaussianl(output.beta, xn)
        plt.plot(xn,yn,color='gold',label='ODR fit '+r'$A\exp\{\frac{-(\theta-\mu)^2}{2\sigma^2}\}+Cst$'+\
          ' \n $A =$%s'%int(par[0])+' $ \pm$ %s'%int(np.sqrt(np.diag(par_va)[0]))+' #'+\
          ' \n $\mu =$ %s'%round(par[1],1)+' $\pm$ %s'%round(np.sqrt(np.diag(par_va)[1]),1)+'°'+\
          ' \n $\sigma =$ %s'%round(par[2],1)+ '$\pm$ %s'%round(np.sqrt(np.diag(par_va)[2]),1)+'°'+\
          ' \n $Cst=$ %s'%int(par[3])+' $\pm$ %s'%int(np.sqrt(np.diag(par_va)[3]))+' #'+\
          ' \n Sum of squares/dof $= $ %s'%round(output.sum_square/(len(Ld)-3),2))
        plt.legend(loc=0)
        plt.errorbar(Ld, res, err_res,err_resx,label=r'$y=\sqrt{counts}$ '+\
                     '\n $x=$%s'%errx+'°',fmt='.', color='black',ecolor='lightgray', elinewidth=3, capsize=0)


    plt.gca().set_xlim(140,220)
    plt.legend(bbox_to_anchor=(0.68, 0.58), loc=1, borderaxespad=0.,prop={'size':14})
    plt.yscale('log')
    plt.ylim(1e2,1e5)



graph(False,1)
graph(True,2)
plt.savefig('Angle analysis')
