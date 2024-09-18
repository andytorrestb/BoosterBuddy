import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
class Rocket:

  def __init__(self, params):

    # Set supplied data
    self.m_i = params['m_i']
    self.m_p = params['m_p']  
    self.I_sp = params['I_sp']
    self.d_ref = params['d_ref']
    self.theta_i = params['theta_i']
    self.t_b = params['t_b']
    self.A = 0.25*math.pi*self.d_ref**2 # assuming circle CSA

    # Set mass flow rate
    self.m_dot = self.m_p / self.t_b

    # Set thrust value (eq X.XX)
    g_o = 9.81
    self.F = self.m_dot*self.I_sp*g_o

    return

  def get_CD_data(self):
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

    return

  def get_CL_data(self):
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

    return


