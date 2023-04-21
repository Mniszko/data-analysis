import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import config

l=config.l
d=config.d
h=config.h
I=config.I

#odwrotność tensora
def sig_xx_value(I,Vc,VH):
    return np.abs(Vc/(Vc**2*l/d+VH**2*d/l)*I)
def sig_xy_value(I,Vc,VH):
    return np.abs(VH/(Vc**2*l**2/d**2+VH**2)*I)

#podobno xx i xy są błędnie nazwane, trzeba się temu lepiej przyjżeć

#funkcje dopasowywane
def sigma2_xx(B, n, mu):
    return (n*mu/(1+mu**2*B**2))*(-1.6e-19)
def sigma2_xy(B, n, mu):
    return (n*mu*mu*B/(1+mu**2*B**2))*(-1.6e-19)


def sigma4_xx(B, n1, mu1,n2,mu2):
    return (n1*mu1/(1+mu1**2*B**2) + n2*mu2/(1+mu2**2*B**2))*1.6e-19
def sigma4_xy(B, n1, mu1,n2,mu2):
    return (n1*mu1*mu1*B/(1+mu1**2*B**2)+n2*mu2*mu2*B/(1+mu2**2*B**2))*1.6e-19

#funkcje dopasowywane ustandaryzowane
def double_fit_2_model(B,n,mu):
    return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)
def double_fit_4_model(B,n1,mu1,n2,mu2):
    return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)



def dopasowanie(B,sdh,h,sig_xx,sig_xy,p0,bounds=([-np.inf, -np.inf], [np.inf, np.inf])):
    #dopasowujemy do 0-1.5 T
    B_lim=np.max(B)
    connected_B = np.concatenate((B,B+B_lim))
    connected_rho = np.concatenate((sdh,h))
    connected_sigma = np.concatenate((sig_xx,sig_xy))
    #modele definiowane wewnątrz, żeby skorzystać z zewnętrznego B_lim, dla wygody w programowaniu

    def double_fit_2_model(B,n,mu):
        return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)
    def double_fit_4_model(B,n1,mu1,n2,mu2):
        return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)

    #ograniczenia lepiej zdefiniować dokładnie jeśli chce się dokładnie dopasować krzywe

    #bounds = ([-np.inf, 0.5e1], [np.inf, 5e2])
    bounds = ([-np.inf, -np.inf], [np.inf, np.inf])


    fit_B = B[(B<1.5)]
    connected_fit_B = np.concatenate((connected_B[:len(fit_B)],connected_B[len(B):len(B)+len(fit_B)]))
    connected_fit_sigma = np.concatenate((connected_sigma[:len(fit_B)],connected_sigma[len(B):len(B)+len(fit_B)]))
    par,cov = curve_fit(double_fit_2_model,connected_fit_B,connected_fit_sigma,p0=p0,bounds=bounds)

    print('parametry modelu 1 nośnikowego (jednostki SI):')
    names = ['n [1/m^2]','mu [m^2/Vs]']
    for i,j in enumerate(par):
        print(f'parametr {names[i]} : {"{:e}".format(j)} +- {"{:e}".format(np.sqrt(cov[i][i]))}')

    plt.plot(connected_B,connected_sigma,c='b',label='dane')
    plt.title('model 1-nośnikowy')
    plt.plot(connected_B,double_fit_2_model(connected_B,*par),c='r',label='dopasowanie')
    plt.legend()
    plt.show()

def parametry_dopasowania(B,sdh,h,sig_xx,sig_xy,p01,p02,bounds=([-np.inf, -np.inf], [np.inf, np.inf])):
    """
    funkcja łączy listy danych B i tensora sigma, zostawia tylko pierwsze 0-1.5 T
    dla każdego sigma mapując dane tak, żeby zgadzały się z poprzednimi, następnie
    fitując krzywe klasyczne tensora sigma do tak wyciętych danych
    """
    #dopasowujemy do 0-1.5 T
    B_lim=np.max(B)
    connected_B = np.concatenate((B,B+B_lim))
    connected_rho = np.concatenate((sdh,h))
    connected_sigma = np.concatenate((sig_xx,sig_xy))
    #modele definiowane wewnątrz, żeby skorzystać z zewnętrznego B_lim, dla wygody w programowaniu

    #model 1-nośnikowy
    def double_fit_2_model(B,n,mu):
        return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)
    #model 2-nośnikowy
    def double_fit_4_model(B,n1,mu1,n2,mu2):
        return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)

    #ograniczenia lepiej zdefiniować dokładnie jeśli chce się dokładnie dopasować krzywe

    #wycinanie danych w przybliżeniu klasycznych
    fit_B = B[(B<1.5)]
    connected_fit_B = np.concatenate((connected_B[:len(fit_B)],connected_B[len(B):len(B)+len(fit_B)]))
    connected_fit_sigma = np.concatenate((connected_sigma[:len(fit_B)],connected_sigma[len(B):len(B)+len(fit_B)]))

    #model 1 nośnikowy
    par1,cov1 = curve_fit(double_fit_2_model,connected_fit_B,connected_fit_sigma,p0=p01,bounds=bounds)
    #model 2 nośnikowy
    par2,cov2 = curve_fit(double_fit_4_model,connected_fit_B,connected_fit_sigma,p0=p02,bounds=bounds)
    return par1, cov1, par2, cov2

def wykres_dopasowania(B,sdh,h,sig_xx,sig_xy,par1,cov1,par2,cov2):
    B_lim=np.max(B)
    connected_B = np.concatenate((B,B+B_lim))
    connected_rho = np.concatenate((sdh,h))
    connected_sigma = np.concatenate((sig_xx,sig_xy))
    def double_fit_2_model(B,n,mu):
        return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)
    def double_fit_4_model(B,n1,mu1,n2,mu2):
        return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)
    
    print('parametry modelu 1 nośnikowego (jednostki SI):')
    names = ['n [1/m^2]','mu [m^2/Vs]']
    for i,j in enumerate(par1):
        print(f'parametr {names[i]} : {"{:e}".format(j)} +- {"{:e}".format(np.sqrt(cov1[i][i]))}')
    print('parametry modelu 1 nośnikowego (jednostki SI):')
    names = ['n [1/m^2]','mu [m^2/Vs]','n [1/m^2]','mu [m^2/Vs]']
    for i,j in enumerate(par2):
        print(f'parametr {names[i]} : {"{:e}".format(j)} +- {"{:e}".format(np.sqrt(cov2[i][i]))}')

    plt.plot(connected_B,connected_sigma,c='b',label='dane')
    plt.title('oba modele')
    plt.plot(connected_B,double_fit_2_model(connected_B,*par1),c='r',label='1-nośnikowy')
    plt.plot(connected_B,double_fit_4_model(connected_B,*par2),c='g',label='2-nośnikowy')
    plt.legend()
    plt.show()
