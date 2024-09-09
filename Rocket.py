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


  def sim_flight(self):
    # Constant variables
    F = self.F
    m=200.0
    psi = 0.0
    theta = 0.0
    Cd = 0.10
    rho = 1.0
    A = 1.0
    g = 10.0

    # time variables
    t_b = 6 # sec
    n_t = 100
    dt = t_b/n_t
    t_vals = np.linspace(0, t_b, n_t)

    # velocity variables
    u_vals = [0]

    #time loop
    for t in t_vals:
      u_k = u_vals[-1]

      # dudt = F/m - (Cd/(2*m))*rho*A*u_k**2 - g
      print(u_k)
      dudt = F/m - (Cd/(2*m))*rho*A*pow(u_k, 2) - g

      u_vals.append(u_k + dudt*dt)



    u_vals.pop()
    return t_vals, u_vals
