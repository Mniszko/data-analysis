"""
---------------------------------------------------
definiowanie funkcji i preambuła
---------------------------------------------------
"""

#predefiniowane funkcje:
from prereq import *
from scipy.optimize import curve_fit

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

"""
---------------------------------------------------
dane próbki i procesowanie danych pomiarowych
---------------------------------------------------
"""

l=4.5e-3 #+-0.5
d=3e-3 #+-0.5
h=20e-9*10
I=1e-6 #zależnie od pomiaru

name = './20230330-002.lvm' #nazwa pliku z pomiarami
values = read_data(name)
numarr = turning_points(values[0],True) #tu powinny pojawić się wszystkie 3 punkty na których trzeba przeciąć dane pomiarowe, czasem trzeba dopisać jeden czy dwa ręcznie
print(numarr)
num1,num2,num3=numarr
B, USdH, UH = main_processing(*values,l,d,h,I,numarr)

"""
---------------------------------------------------
rysowanie danych
---------------------------------------------------
"""

#oporności
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 12))

ax1.set_title('SdH')
ax1.plot(B,np.abs(USdH),c='r',label='SdH')
ax1.set_xlim(0,9)
ax1.legend()
ax2.set_title('Hall')
ax2.plot(B,np.abs(UH),c='r',label='Hall')
ax2.set_xlim(0,9)
ax2.legend()

fig.show()

#przewodności
plt.plot(B,sig_xx_value(1e-6,USdH,UH),label='sig_xx',c='r')
plt.xlim(0,2)
plt.ylim(0,0.003)
plt.title('\sigma_{xx}')
plt.legend()
plt.show()

plt.plot(B,sig_xy_value(1e-6,USdH,UH)/10,label='sig_xy',c='g')
plt.ylim(0,0.0017)
plt.title('\sigma_{xy}')
plt.legend()
plt.show()


"""
---------------------------------------------------
funkcja dopasowująca i rysująca wykresy
---------------------------------------------------
"""

def dopasowanie(B,sdh,h,sig_xx,sig_xy,p0):
    #dopasowujemy do 0-1.5 T
    B_lim=np.max(B)
    connected_B = np.concatenate((B,B+B_lim))
    connected_rho = np.concatenate((sdh,h))
    connected_sigma = np.concatenate((sig_xx,sig_xy))
    #modele definiowane wewnątrz, żeby skożystać z zewnętrznego B_lim, dla wygody w programowaniu

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

"""
---------------------------------------------------
przykładowe dopasowanie
---------------------------------------------------
"""


p0=[-1.097695e+16,1.242392e+00] #pierwsze zgadnięcia parametrów
#trzeba ograniczyć B do 9 tesli max
B = B[(B <= 9)]
USdH = USdH[:len(B)]
UH = UH[:len(B)]
dopasowanie(B,UH,USdH,sig_xx_value(1e-6,USdH,UH),sig_xy_value(1e-6,USdH,UH),p0) 
