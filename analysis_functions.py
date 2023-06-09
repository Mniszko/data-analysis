import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import least_squares

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

def parametry_dopasowania(B,sdh,h,sig_xx,sig_xy,p01,p02,bounds1=([-np.inf, -np.inf], [np.inf, np.inf]),bounds2=([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),low_cut=0,high_cut=1.5):
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

    fit_B = B[(B >= low_cut) & (B <= high_cut)]
    connected_fit_B = np.concatenate((fit_B,fit_B+B_lim))
    connected_fit_sigma = np.concatenate(
        (sig_xx[(B >= low_cut) & (B <= high_cut)],
        sig_xy[(B >= low_cut) & (B <= high_cut)])
    )

    #model 1 nośnikowy
    par1,cov1 = curve_fit(double_fit_2_model,connected_fit_B,connected_fit_sigma,p0=p01,bounds=bounds1)
    #model 2 nośnikowy
    par2,cov2 = curve_fit(double_fit_4_model,connected_fit_B,connected_fit_sigma,p0=p02,bounds=bounds2)
    return par1, cov1, par2, cov2

def parametry_dopasowania_ls(B,sdh,h,sig_xx,sig_xy,p01,p02,bounds1=([-np.inf, -np.inf], [np.inf, np.inf]),bounds2=([-np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf]),low_cut=0,high_cut=1.5):

    B_lim=np.max(B)
    connected_B = np.concatenate((B,B+B_lim))
    connected_rho = np.concatenate((sdh,h))
    connected_sigma = np.concatenate((sig_xx,sig_xy))

    def double_fit_2_model(B,n,mu):
        return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)

    def double_fit_4_model(B,n1,mu1,n2,mu2):
        return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)

    fit_B = B[(B >= low_cut) & (B <= high_cut)]
    connected_fit_B = np.concatenate((fit_B,fit_B+B_lim))
    connected_fit_sigma = np.concatenate(
        (sig_xx[(B >= low_cut) & (B <= high_cut)],
        sig_xy[(B >= low_cut) & (B <= high_cut)])
    )

    # Residuals for the two models
    def residuals_1(p, B, sigma):
        return sigma - double_fit_2_model(B, *p)
    
    def residuals_2(p, B, sigma):
        return sigma - double_fit_4_model(B, *p)

    # Fitting with least_squares
    result_1 = least_squares(residuals_1, p01, bounds=bounds1, args=(connected_fit_B, connected_fit_sigma))
    result_2 = least_squares(residuals_2, p02, bounds=bounds2, args=(connected_fit_B, connected_fit_sigma))

    return result_1, result_2

def wykres_dopasowania_ls(B, sig_xx, sig_xy, result_1, result_2):
    plt.figure()

    B_lim = np.max(B)
    connected_B = np.concatenate((B,B+B_lim))

    def double_fit_2_model(B,n,mu):
        return sigma2_xx(B,n,mu)*np.heaviside(-B+B_lim,0.5) + sigma2_xy(B-B_lim,n,mu)*np.heaviside(B-B_lim,0.5)

    def double_fit_4_model(B,n1,mu1,n2,mu2):
        return sigma4_xx(B,n1,mu1,n2,mu2)*np.heaviside(-B+B_lim,0.5) + sigma4_xy(B-B_lim,n1,mu1,n2,mu2)*np.heaviside(B-B_lim,0.5)

    # Generate fitted data
    fit_1 = double_fit_2_model(connected_B, *result_1.x)
    fit_2 = double_fit_4_model(connected_B, *result_2.x)

    # Plot initial data
    plt.plot(B, sig_xx, c='b', label='Initial data - sig_xx')
    plt.plot(B+B_lim, sig_xy, c='b', label='Initial data - sig_xy')

    # Plot fitted data
    plt.plot(connected_B, fit_1, 'r', label='model 1 nośnikowy')
    plt.plot(connected_B, fit_2, 'g', label='model 2 nośnikowy')

    plt.xlabel('B')
    plt.ylabel('sigma')
    plt.legend()

    print('parametry modelu 1 nośnikowego (jednostki SI):')
    names = ['n [1/m^2]','mu [m^2/Vs]']
    for i,j in enumerate(result_1.x):
        print(f'parametr {names[i]} : {"{:e}".format(j)}')
    print('parametry modelu 1 nośnikowego (jednostki SI):')
    names = ['n [1/m^2]','mu [m^2/Vs]','n [1/m^2]','mu [m^2/Vs]']
    for i,j in enumerate(result_2.x):
        print(f'parametr {names[i]} : {"{:e}".format(j)}') #macierz kowariancji trzeba by było inaczej policzyć

    plt.show()
    


def wykres_dopasowania(B,sdh,h,sig_xx,sig_xy,par1,cov1,par2,cov2,fig=None,ax=None):
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

    if fig==None and ax==None:
        plt.plot(connected_B,double_fit_2_model(connected_B,*par1),c='r',label='1-nośnikowy')
        plt.plot(connected_B,double_fit_4_model(connected_B,*par2),c='g',label='2-nośnikowy')
        plt.plot(connected_B,connected_sigma,c='b',label='dane')
        plt.legend()
    elif fig!=None and ax!=None:
        ax.plot(connected_B,double_fit_2_model(connected_B,*par1),c='r',label='1-nośnikowy')
        ax.plot(connected_B,double_fit_4_model(connected_B,*par2),c='g',label='2-nośnikowy')
        ax.plot(connected_B,connected_sigma,c='b',label='dane')
        ax.legend()
    else:
        raise Exception("funkcja wykres_dopasowania ma źle zdefiniowane zmienne fig i ax")
    
def transformata_fouriera(B,sdh,cut_amplitude,cut_high_amp,do_plot=True):
    #skrajne częstotliwości wycinanego przedziału
    CUT_AMPLITUDE = cut_amplitude
    CUT_HIGH_AMP = cut_high_amp
    def zero_cut(x,y,cut_amp):
        return (np.heaviside(x-cut_amp,0.5) + np.heaviside(-x - cut_amp,0.5))*y
    def high_cut(x,y,cut_amp):
        return (np.heaviside(-x+cut_amp,0.5) + np.heaviside(x + cut_amp,0.5)-1)*y

    #sygnał i interpolacja
    x_space = 1/B
    signal = sdh

    signal_x_space = np.linspace(np.min(x_space),20,9000)
    signal = np.interp(signal_x_space,np.flip(x_space),np.flip(signal))

    #transformata fouriera
    signal_frequencies = np.fft.fft(signal)
    original_sig_freq = signal_frequencies.copy()
    freq_points = np.fft.fftfreq(len(signal_x_space),signal_x_space[1]-signal_x_space[0])

    #wycinanie przedziału częstotliwości
    signal_frequencies=zero_cut(freq_points,np.abs(signal_frequencies),CUT_AMPLITUDE)
    signal_frequencies=high_cut(freq_points,np.abs(signal_frequencies),CUT_HIGH_AMP)

    #odwrotna transformata
    retrived_signal = np.fft.ifft(signal_frequencies)

    return signal_x_space,signal,signal_frequencies,original_sig_freq,freq_points,retrived_signal

def wykres_fouriera(signal_x_space,signal,signal_frequencies,original_sig_freq,freq_points,retrived_signal):

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(10, 12))
    # przestrzeń sygnału
    ax1.plot(signal_x_space, signal, label='signal (sdh)',c='r')
    ax1.set_xlabel('1/$B$ [T$^{-1}$]')
    ax1.set_ylabel('Amplitude')
    #ax1.set_xlim(0,0.8)
    #ax1.set_xticks([0,0.12,0.178])
    ax1.grid()
    ax1.legend()
    # przestrzeń częstotliwości
    #ax2.plot(freq_points, np.abs(signal_frequencies), label='cropped fourier transform (sdh)',c='b')
    ax2.plot(freq_points, np.abs(original_sig_freq), label='original fourier transform (sdh)',c='salmon')
    ax2.set_xlabel('Frequency ($B$)')
    ax2.set_ylabel('Magnitude')
    #ax2.set_xlim(10,25)
    #ax2.set_ylim(0.0001,0.004)
    #ax2.set_xticks([10,11,12,13,14,15,16,17,18,19,20])
    ax2.legend()
    # odwrotna transformacja po wydzięciu częstotliwości
    #ax3.plot(signal_x_space,retrived_signal,label='inverse fourier transform (sdh)',c='b')
    ax3.plot(signal_x_space, signal, label='signal (sdh)',c='salmon')
    ax3.set_xlabel('1/$B$ [T$^{-1}$]')
    ax3.set_ylabel('Amplitude')
    ax3.set_xlim(0.1,0.5)
    ax3.legend()
    fig.tight_layout()
    return fig, (ax1, ax2, ax3)

def koncentracja_z_fouriera(Bperiod,Bfrequency):
    #print('na podstawie odległości maksimów (0.058 1/T) częstotliwość wynosi ~17.2 T')
    planck = 6.626e-34 #Js
    elementary_charge = 1.6e-19 #C
    n = 2*elementary_charge/Bperiod/planck
    print(f'koncentracja nośników n={"{:e}".format(n)}')
    print('na podstawie odczytów z transformacji fouriera:')
    n = 2*elementary_charge*Bfrequency/planck
    print(f'koncentracja nośników n={"{:e}".format(n)}')

def SdH_without_straight_line(x,y):
    def line_model(x,a,b):
        return a*x+b
    par,cov = curve_fit(line_model,x,y)
    y2 = y - line_model(x,*par)
    return y2

#sig_xx[(B >= low_cut) & (B <= high_cut)]
#fit_B = B[(B >= low_cut) & (B <= high_cut)]