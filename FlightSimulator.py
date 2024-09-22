import numpy as np
import math
class FlightSimulator:

  def __init__(self, rocket, atmosphere):
    self.rocket = rocket
    self.atmosphere = atmosphere
    return

  def set_dt(self, t_b, n_t):
    return

  def run_1D(self):
    # Constant variables.
    F = self.rocket.F
    m_dot = self.rocket.m_dot
    psi = 0.0
    theta = 0.0
    A = self.rocket.A
    g = 10.0

    # Time variables.
    t_b = self.rocket.t_b
    n_t = 500 
    dt = t_b/n_t
    t_vals = np.linspace(0, t_b, n_t)

    # Velocity values  (wrt to flight).
    u_vals = []
    u_o = 0

    # Altitude values.
    h_vals = []
    h_o = 0

    # Range values.
    x_vals =  []
    x_o = 0

    # Mass values.
    m_vals = []

    # Atmospheric data
    rho_vals = []
    P_vals = []
    T_vals = []

    # Aerodynamic data 
    Ma_vals = []
    CD_vals = []
    CL_vals = []

    #time loop
    for t in t_vals:
      # Save data from previous iteration. 
      # Use initial conditions of needed.
      if len(u_vals) == 0:
        u_k = u_o
        h_k = h_o
      else:
        u_k = u_vals[-1] # velocity
        h_k = h_vals[-1] 

      print(u_k, h_k)
      print('========')

      # Gather atmospheric data.
      # Get P, T, rho, for a given altitude
      atm_sample = self.atmosphere.get_atm_props(h_k)
      P = atm_sample['PPSL']*101325
      rho = atm_sample['rho']
      T = atm_sample['T']

      # Calculate Mach Number 
      a_k = np.sqrt(1.4*287*T)
      Ma_k = u_k / a_k

      # Calculate aerodynamic coeeficients
      CD_k = 0

      # Update rocket mass.
      if len(m_vals) == 0:
        m_k = self.rocket.m_i - m_dot*dt
      else:
        m_k = m_vals[-1] - m_dot*dt

      # Calculate rocket dynamics
      dudt = F/m_k - (CD_k/(2*m_k))*rho*A*pow(u_k, 2) - g

      u = u_k + dudt*dt
      h = h_k + u_k*dt + 0.5*dudt*pow(dt, 2)

      # Save rocket dynamics to an arrays.

      u_vals.append(u)
      h_vals.append(h)
      m_vals.append(m_k)

      # Save aerodynamic data to arrays.
      Ma_vals.append(Ma_k)
      CD_vals.append(CD_k)
      # CL_vals.append(CL_k)

      # Save atmopspheric data to arrays.
      rho_vals.append(rho)
      P_vals.append(P)
      T_vals.append(T)

    results = {
      't': t_vals,
      'u': u_vals,
      'h': h_vals,
      'rho': rho_vals,
      'P': P_vals,
      'T': T_vals,
      'Ma': Ma_vals
    }

    return results

  def run_2D(self):
      # Constant variables.
      F = self.rocket.F
      m_dot = self.rocket.m_dot
      psi = 0.0
      theta = 0.0
      A = self.rocket.A
      g = 10.0

      # Time variables.
      t_b = self.rocket.t_b
      n_t = 500 
      dt = t_b/n_t
      t_vals = np.linspace(0, t_b, n_t)

      # Velocity values (wrt to flight path).
      u_vals = [] # parallel to path
      u_o = 0

      v_vals = [] # perp to path
      v_o = 0

      # Flight velocities mapped to global reference.
      ux_vals = []
      uy_vals = []
      ux_o = 0
      uy_o = 0

      vx_vals = []
      vy_vals = []
      vx_o = 0
      vy_o = 0

      # Velocities as defined wrt to global reference.
      V_vals = []
      V_o = 0

      Vx_vals = []
      Vx_o = 0

      Vy_vals = []
      Vy_o = 0

      # Altitude values.
      h_vals = []
      h_o = 0

      # Range values.
      r_vals =  []
      r_o = 0

      # Mass values.
      m_vals = []

      # Atmospheric data.
      rho_vals = []
      P_vals = []
      T_vals = []

      # Aerodynamic data.
      Ma_vals = []
      CD_vals = []
      CL_vals = []

      # Flight angle values.
      theta_vals = []
      theta_o = 70  * (math.pi/180)

      #time loop
      for t in t_vals:
        # Save data from previous iteration. 
        # Use initial conditions if needed.
        if len(u_vals) == 0:
          u_k = u_o
          ux_k = ux_o
          uy_k = uy_o
          h_k = h_o
          r_k = r_o
          v_k = v_o
          V_x = Vx_o
          V_y = Vy_o
          theta_k = theta_o

        else:
          u_k = ux_vals[-1] # velocity
          h_k = h_vals[-1] 
          ux_k = ux_vals[-1]
          uy_k = uy_vals[-1]
          h_k = h_vals[-1]
          r_k = r_vals[-1]
          v_k = v_vals[-1]
          Vx_k = Vx_vals[-1]
          Vy_k = Vy_vals[-1]
          theta_k = theta_vals[-1]

        print(u_k, h_k)
        print('========')

        # Gather atmospheric data.
        # Get P, T, rho, for a given altitude
        atm_sample = self.atmosphere.get_atm_props(h_k)
        P = atm_sample['PPSL']*101325
        rho = atm_sample['rho']
        T = atm_sample['T']

        # Calculate Mach Number 
        a_k = np.sqrt(1.4*287*T)
        Ma_k = u_k / a_k

        # Calculate aerodynamic coefficients
        CD_k = self.rocket.get_CD(Ma_k)
        CL_k = self.rocket.get_CL(Ma_k)
        # CD_k = 0
        # CL_k = 0

        # Update rocket mass.
        if len(m_vals) == 0:
          m_k = self.rocket.m_i - m_dot*dt
        else:
          m_k = m_vals[-1] - m_dot*dt

        # Calculate flight derivatives.
        dudt = F/m_k - (CD_k/(2*m_k))*rho*A*pow(u_k, 2) - g*math.sin(theta)
        dvdt = (CD_k/(2*m_k))*rho*A*pow(u_k,2)- g*math.cos(theta)
        
        # Force no rotation if on the launch pad
        if h_k < 100:
          dvdt = 0

        # For no rotation if previous velocity is zero.
        if u_k == 0:
          dthetadt = 0
        else:
          dthetadt = dvdt/u_k

        # Update flight velocity and variables. 
        v = v_k + dvdt*dt
        u = u_k + dudt*dt

        # Handle launch pad logic.
        if h_k < 100:
          theta = theta_o

          # Map velocity components to global frame of reference.
          ux = u*math.cos(theta)
          uy = u*math.sin(theta)
          
          Vx = ux
          Vy = uy
        else:
          # Update theta for current time step.
          theta = theta_k + dthetadt*dt

          # Map velocity components to global frame of reference (parallel to flight).
          ux = u*math.cos(theta)
          uy = u*math.sin(theta)

          # Map velocity components to global frame of reference (perpendicular to flight).
          theta_u = theta
          theta_v = theta_u + 0.5*math.pi
          vx = v*math.cos(theta_v)
          vy = v*math.sin(theta_v)

          # Velocity wrt to glopbal reference frame.
          Vx = ux + vx
          Vy = uy + vy

        # Calculate changes in altitude and range.
        # if len(u_vals) == 0:
        #   dvydt = Vy / dt
        # else:
        #   dvydt = (Vy - Vy_vals[-1]) / dt

        dvydt = dudt*math.sin(theta)

        h = h_k + u_k*dt + 0.5*dvydt*pow(dt, 2)
        r = r_k + Vx*dt

        # Save rocket dynamics to an arrays.
        theta_vals.append(theta)
        v_vals.append(v)
        u_vals.append(u)
        ux_vals.append(ux)
        uy_vals.append(uy)
        Vx_vals.append(Vx)
        Vy_vals.append(Vy)
        h_vals.append(h)
        r_vals.append(r)
        m_vals.append(m_k)

        # Save aerodynamic data to arrays.
        Ma_vals.append(Ma_k)
        CD_vals.append(CD_k)
        CL_vals.append(CL_k)

        # Save atmopspheric data to arrays.
        rho_vals.append(rho)
        P_vals.append(P)
        T_vals.append(T)

      results = {
        't': t_vals,
        'ux': ux_vals,
        'uy': uy_vals,
        'u': u_vals,
        'Vy': Vy_vals,
        'Vx': Vx_vals,
        'h': h_vals,
        'r': r_vals,
        'rho': rho_vals,
        'P': P_vals,
        'T': T_vals,
        'Ma': Ma_vals
      }

      return results