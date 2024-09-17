import numpy as np

class Rocket:

  def __init__(self, params):

    # Set supplied data
    self.m_i = params['m_i']
    self.m_p = params['m_p']  
    self.I_sp = params['I_sp']
    self.d_ref = params['d_ref']
    self.theta_i = params['theta_i']
    self.t_b = params['t_b']

    # Set mass flow rate
    self.m_dot = self.m_p / self.t_b

    # Set thrust value (eq X.XX)
    g_o = 9.81
    self.F = self.m_dot*self.I_sp*g_o

    return



