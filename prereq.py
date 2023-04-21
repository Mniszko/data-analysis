import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import locale


plt.rc('text', usetex=False)
plt.rcParams['figure.figsize']=[12,8]
plt.rc('font', size=10*4/2)          # controls default text sizes
plt.rc('axes', titlesize=10*4/2)     # fontsize of the axes title
plt.rc('axes', labelsize=10*4/2)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=8*4/2)    # fontsize of the tick labels
plt.rc('ytick', labelsize=8*4/2)    # fontsize of the tick labels
plt.rc('legend', fontsize=8*4/2)    # legend fontsize
plt.rc('figure', titlesize=16*4/2)
locale.setlocale(locale.LC_NUMERIC, "pl_PL.UTF-8")
plt.rcParams['axes.formatter.use_locale']=True




#definicje klas i funkcji

def simplifier(arr,cut_last=True,cut_first=True):
    chng = arr[0]
    if cut_first:
        finarr = []
    else:
        finarr = []
    for itr in arr[1:]:
        if abs(chng-itr) > 10:
            finarr.append(itr)
        chng = itr

    if cut_last:
        finarr = finarr[:-1]
    
    return finarr

def turning_points_homebrew(B,auto_simplify=False):
    znak = 1
    arri=[]
    for i in range(len(B)-1):
        ZNAK_TEMP = 1
        x1 = B[i]
        x2 = B[i+1]
        chng = x2-x1
        if chng == 0:
            ZNAK_TEMP = znak
        elif chng == abs(chng):
            ZNAK_TEMP = 1
        elif chng != abs(chng):
            ZNAK_TEMP = -1
        
        if (ZNAK_TEMP != znak):
            arri.append(i)
        
        znak = ZNAK_TEMP

    if auto_simplify:
        arri = simplifier(arri)
    return arri

def turning_points(array,auto_simplify=False):
    ''' turning_points(array) -> min_indices, max_indices
    Finds the turning points within an 1D array and returns the indices of the minimum and 
    maximum turning points in two separate lists.
    '''
    idx_max = []
    if (len(array) < 3): 
        return idx_max

    NEUTRAL, RISING, FALLING = range(3)
    def get_state(a, b):
        if a < b: return RISING
        if a > b: return FALLING
        return NEUTRAL

    ps = get_state(array[0], array[1])
    begin = 1
    for i in range(2, len(array)):
        s = get_state(array[i - 1], array[i])
        if s != NEUTRAL:
            if ps != NEUTRAL and ps != s:
                if s == FALLING: 
                    idx_max.append((begin + i - 1) // 2)
                else:
                    idx_max.append((begin + i - 1) // 2)
            begin = i
            ps = s
    if auto_simplify:
        idx_max = simplifier(idx_max)
    return idx_max

def extract_and_interpolate(xfield,variable,num1,num2,num3,Low=0,High=10,Xnumber=1000):

    xfield1 = xfield[:num1]
    variable1 = variable[:num1]
    xfield2 = -xfield[num1:num2]
    variable2 = variable[num1:num2]
    xfield3 = xfield[num2:num3]
    variable3 = variable[num2:num3]
    xfield4 = -xfield[num3:]
    variable4 = variable[num3:]

    x1 = np.linspace(Low,High,Xnumber)
    y1 = np.interp(x1,xfield1,variable1)

    x2 = -np.linspace(Low,High,Xnumber)
    y2 = np.interp(x2,xfield2,variable2)

    x3 = np.linspace(Low,High,Xnumber)
    y3 = np.interp(x3,xfield3,variable3)

    x4 = -np.linspace(Low,High,Xnumber)
    y4 = np.interp(x4,xfield4,variable4)
    return [y1,y2,y3,y4]

def read_data(filename):
    dane = pd.read_csv(filename,delimiter='\t',skiprows=21)
    H = dane['DC Voltage 27 [mV]']/1e3
    SdH = dane['DC Voltage [mV]']/1e3
    B = dane['B field [T]']
    return B, SdH, H


def main_processing(B,SdH,H,l,d,h,I,numarr):
    num1, num2, num3 = numarr
    

    B1 = B[:num1]
    SdH1 = SdH[:num1]
    H1 = H[:num1]

    B2 = B[num1:num2]
    SdH2 = SdH[num1:num2]
    H2 = H[num1:num2]

    #pole ujemne:
    B3 = B[num2:num3]
    SdH3 = SdH[num2:num3]
    H3 = H[num2:num3]

    B4 = B[num3:]
    SdH4 = SdH[num3:]
    H4 = H[num3:]



    #interpolacja:
    B1_int = np.linspace(0.01,9.6,1000)
    H1_int = np.interp(B1_int,B1,H1)
    SdH1_int = np.interp(B1_int,B1,SdH1)

    B2_int = np.linspace(0.01,9.6,1000)
    H2_int = np.interp(B2_int,np.flip(B2),np.flip(H2))
    SdH2_int = np.interp(B2_int,np.flip(B2),np.flip(SdH2))

    B3_int = np.linspace(0.01,9.6,1000)
    H3_int = np.interp(B3_int,B3,H3)
    SdH3_int = np.interp(B3_int,B3,SdH3)

    B4_int = np.linspace(0.01,9.6,1000)
    H4_int = np.interp(B4_int,np.flip(B4),np.flip(H4))
    SdH4_int = np.interp(B4_int,np.flip(B4),np.flip(SdH4))

    #pełne dane
    USdH = (SdH1_int - SdH2_int - SdH3_int + SdH4_int)/4
    UH = (H1_int-H2_int+H3_int-H4_int)/4

    return B1_int, USdH, UH