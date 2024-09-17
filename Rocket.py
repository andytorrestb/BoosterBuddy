import numpy as np

class Rocket:

  def __init__(self, params):

    self.m_i = params['m_i']
    self.m_p = params['m_p']  
    self.I_sp = params['I_sp']
    self.m_dot = params['m_dot']
    self.d_ref = params['d_ref']
    self.theta_i = params['theta_i']


    return

  def set_rocket_thrust(self):
    """
        F=m_dot*I_sp*g_o (eq X.XX)
    """
    g_o = 9.81
    self.F = self.m_dot*self.I_sp*g_o

