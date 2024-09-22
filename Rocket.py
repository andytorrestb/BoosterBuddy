import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def interpolate(Ma, x2, x1, alpha):
    dxdh = (x2[alpha]- x1[alpha]) / (x2['x'] - x1['x'])

    x = x1[alpha] + (Ma - x1['x']) * dxdh

    return float(x)

class Rocket:

  def __init__(self, params):

    # Set supplied data
    self.m_i = params['m_i']
    self.m_p = params['m_p']  
    self.I_sp = params['I_sp']
    self.theta_i = params['theta_i']
    # self.t_b = params['t_b']
    self.m_dot = params['m_dot']


    # Set burn time
    self.t_b = self.m_p / self.m_dot
    self.n_t = params['n_t']
    # Set geometry calculations.
    self.d_ref = params['d_ref']
    self.A = 0.25*math.pi*self.d_ref**2 # assuming circle CSA


    # Set thrust value (eq X.XX)
    g_o = 9.81
    self.F = self.m_dot*self.I_sp*g_o

    return

  def set_CD_data(self):
    CD = pd.read_csv('CD_vs_Ma.csv')

    Ma = CD['x']
    plt.plot(Ma, CD['alpha10'], label = 'alpha = 10')
    plt.plot(Ma, CD['alpha08'], label = 'alpha = 08')
    plt.plot(Ma, CD['alpha06'], label = 'alpha = 06')
    plt.plot(Ma, CD['alpha04'], label = 'alpha = 04')
    plt.plot(Ma, CD['alpha00'], label = 'alpha = 00')

    plt.xlabel('Mach Number')
    plt.ylabel('Drag Coefficient')

    plt.savefig('CD_vs_Ma.png')
    plt.clf()

    self.CD = CD
    return

  def set_CL_data(self):
    CL = pd.read_csv('CL_vs_Ma.csv')

    Ma = CL['x']

    plt.plot(Ma, CL['alpha10'], label = 'alpha = 10')
    plt.plot(Ma, CL['alpha08'], label = 'alpha = 08')
    plt.plot(Ma, CL['alpha06'], label = 'alpha = 06')
    plt.plot(Ma, CL['alpha04'], label = 'alpha = 04')
    plt.plot(Ma, CL['alpha00'], label = 'alpha = 00')

    plt.xlabel('Mach Number')
    plt.ylabel('Lift Coefficient')

    plt.savefig('CL_vs_Ma.png')
    plt.clf()

    self.CL = CL
    return

  def get_CL(self, Ma):
    # if Ma == 0:
    #   print("object not moving (Mach not set)")
    #   return
    # print(self.CL)
    CL = self.CL
    # Find rows of data used for interpolation.
    for i in range(len(CL)):
        # Save current row for reference.
        row  = CL.iloc[i]

        # If h is found, return data for that row. 
        if row['x'] == Ma:
            return row.to_dict()

        # If we pass h, save data for interpolation.
        if row['x'] > Ma:
            x2 = CL.iloc[i]
            x1 = CL.iloc[i-1]
            break
        
        # Rocket has left the atmopsphere.
        if i == len(CL) - 1:
            print("Mach number not supported")

    return interpolate(Ma, x2, x1, 'alpha00')
        


  def get_CD(self, Ma):
    CD = self.CD
    # Find rows of data used for interpolation.
    for i in range(len(CD)):
        # Save current row for reference.
        row  = CD.iloc[i]

        # If h is found, return data for that row. 
        if row['x'] == Ma:
            return row.to_dict()['alpha00']

        # If we pass h, save data for interpolation.
        if row['x'] > Ma:
            x2 = CD.iloc[i]
            x1 = CD.iloc[i-1]
            break
        
        # Rocket has left the atmopsphere.
        if i == len(CD) - 1:
            print("Mach number not supported")

    return interpolate(Ma, x2, x1, 'alpha00')



