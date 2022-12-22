#SergioLondono_AG1 Exoplanetas y Astrobiologia.

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

#Define constants 
M_s = 0.45          # Star mass in solar mass  
R_s = 0.46          # Star radius in solar radius
T_s = 3350           # Star temperature in Kelvin
Period = 2.6439     # Planet Period in days
exent = 0.16        # Exentricity 

## Read the data.
RV_Hires = pd.read_csv('Datos/RV_HIRES.csv', sep=";", decimal=",")              #Radial speed in m/s
Transito_Sp = pd.read_csv('Datos/transito_Spitzer.csv',sep=";", decimal=",")    #Flux vs time

#linea de tendenci and Deap transit
Transito_Sp['averange'] = Transito_Sp['Flux'].rolling(10).mean()
TD_max= max (Transito_Sp.iloc[10:,3])
TD_min= min (Transito_Sp.iloc[10:,3])
TD=(TD_max-TD_min)/TD_max
print('transit deep: ',TD)

# Planet radius
R_ps=(TD*R_s**2)**(0.5)                                                         #Planet radius in solars radius
R_p=R_ps*695700                                                                 #Planet radius in km
print('Planer radius: ',R_p)

# planet orbit radius 
Period_sid=Period/365.25636                                                     #Period in sideral year    
a_p = (M_s*Period_sid**2)**(1/3)                                                    #Planer orbit radius in U.A.
print ('Orbit radius: ',a_p)

# curve_fit to RV
def fun(x,a,b,c):
    return a*np.cos(b*x+c)
RV_ajus, cov = curve_fit(fun,RV_Hires.iloc[:,3],RV_Hires.iloc[:,1])
x_fit=np.linspace(0,1,1000)
RV_sin= fun(x_fit,RV_ajus[0],RV_ajus[1],RV_ajus[2])

# Estimation of Period and amplitude
P=abs(2*np.pi/RV_ajus[1])                                                       #Period in phase.
K=abs(RV_ajus[0])                                                               #Amplitude in [m/S]
print('K: ',K,'P: ', P)

#Minimun planet mass
mass_minP_J = 4.93E-3*(1-exent)**(0.5)*K*Period**(1/3)*M_s**(2/3)               #Planet mass in terms od jupyter mass
mass_minP=mass_minP_J*1.898E27                                                  #planet mas in kg 
print('planet mass[Mj]',mass_minP_J,'planet mass[kg]',mass_minP)

#uncertainty of 10%
dmass=4.93E-3*(1-exent)**(0.5)*K*Period**(1/3)*M_s**(-1/3)*(2/3)*0.1
print('uncertainy',dmass)

#densidad del planeta
rho_P=mass_minP*1000/(4/3*np.pi*(R_p*100000)**3)
print('dencity [g/cm^3]: ',rho_P)

#Equilibrium temperature
T_eq=T_s*(1-0)**(1/4)*((R_s*695700)/(a_p*1.496e+8))**(0.5)
print('Teq:',T_eq)

plt.figure('radial velocity vs time')
plt.errorbar(RV_Hires.iloc[:,0],RV_Hires.iloc[:,1],RV_Hires.iloc[:,2],fmt='o',ms=3,capsize=2)
plt.xlabel('time [Days]')
plt.ylabel('Radial velocity [m/s]')

plt.figure('radial velocity vs phase')
plt.subplots
plt.errorbar(RV_Hires.iloc[:,3],RV_Hires.iloc[:,1],RV_Hires.iloc[:,2],fmt='o',ms=3,capsize=2)
plt.plot(x_fit,RV_sin)
plt.xlabel('Phase')
plt.ylabel('Radial velocity [m/s]')

plt.figure('Transite data with MA steps 10')
plt.errorbar(Transito_Sp.iloc[:,0],Transito_Sp.iloc[:,1],Transito_Sp.iloc[:,2], fmt='o',ms=3,capsize=2)
plt.plot(Transito_Sp['HJD'],Transito_Sp['averange'])
plt.xlabel('Time in days')
plt.ylabel('Flux in %')
# plt.show()