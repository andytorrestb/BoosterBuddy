import numpy as np

class FlightSimulator:

  def __init__(self, rocket, atmosphere):
    self.rocket = rocket
    self.atmosphere = atmosphere
    return

  def set_dt(self, t_b, n_t):
    return


  def run(self):
    # Constant variables.
    F = self.rocket.F
    m= self.rocket.m_i
    m_dot = self.rocket.m_dot
    psi = 0.0
    theta = 0.0
    Cd = 0.10
    rho = 1.0
    A = 1.0
    g = 10.0

    # Time variables.
    t_b = 60 # sec
    n_t = 500
    dt = t_b/n_t
    t_vals = np.linspace(0, t_b, n_t)

    # Velocity values.
    u_vals = [0]

    # Altitude values.
    h_vals = [0]

    # mass values
    m_vals = [self.rocket.m_i]

    #time loop
    for t in t_vals:
      # Save data from previous iteration.
      u_k = u_vals[-1] # velocity
      h_k = h_vals[-1]
      print(u_k)

      # Gather atmospheric data.
      # Get P, T, rho, for a given altitude



      # Calculate rocket dynamics
      dudt = F/m - (Cd/(2*m))*rho*A*pow(u_k, 2) - g*dt
      
      u = u_k + dudt*dt
      h = h_k + u_k*t_b + 0.5*dudt*pow(dt, 2)
      
      # Update rocket properties.
      m = m_vals[-1] - m_dot*dt

      # Save to results array.
      u_vals.append(u)
      h_vals.append(h)
      m_vals.append(m)

    


    u_vals.pop()
    h_vals.pop()
    return t_vals, u_vals, h_vals
