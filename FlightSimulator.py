
import numpy as np
import math

class FlightSimulator:
    # This class defines a flight simulator for a rocket considering atmospheric effects.

    def __init__(self, rocket, atmosphere):
        # Initialize the simulator with a rocket object and an atmosphere object.
        self.rocket = rocket
        self.atmosphere = atmosphere
        self.g = 10.0  # Set the acceleration due to gravity (m/s^2), approximated as 10 m/s^2 for simplicity.  # gravity constant
    
    def set_dt(self):
        # Set the time step based on rocket's burn time and number of time steps.
        self.dt = self.rocket.t_b / self.rocket.n_t  # Calculate the time step (dt) by dividing the total burn time by the number of time steps.
        self.t_vals = np.linspace(0, self.rocket.t_b, self.rocket.n_t + 1)  # Generate an array of time values from 0 to burn time, equally spaced over the number of steps.

    def initialize_values(self):
        # Initialize all required lists for storing simulation results and set initial conditions.
        # Start all lists empty and set initial conditions within the simulation loops
        self.u_vals = []  # Initialize the list to store velocity values at each time step.
        self.h_vals = []  # Initialize the list to store altitude values at each time step.
        self.r_vals = []
        self.x_vals = []
        self.m_vals = []
        self.rho_vals, self.P_vals, self.T_vals = [], [], []
        self.Ma_vals, self.CD_vals, self.CL_vals = [], [], []
        self.theta_vals = []  # for 2D simulation
        self.ux_vals = []  # for 2D simulation
        self.uy_vals = []  # for 2D simulation

    def update_atmosphere(self, altitude):
        # Update atmospheric properties based on current altitude and store the results.
        atm_props = self.atmosphere.get_atm_props(altitude)
        P = atm_props['PPSL'] * 101325
        rho = atm_props['rho']
        T = atm_props['T']
        return P, rho, T  # Return the pressure (P), density (rho), and temperature (T) at the given altitude.

    def run_1D(self):
        # Run a 1D simulation of the rocket's flight using vertical dynamics only.
        self.initialize_values()
        A = self.rocket.A
        u_o = 0
        h_o = 0
        m_o = self.rocket.m_i

        for t in self.t_vals:
            if len(self.u_vals) == 0:
                u_k, h_k, m_k = u_o, h_o, m_o
            else:
                u_k = self.u_vals[-1]
                h_k = self.h_vals[-1]
                m_k = self.m_vals[-1]

            atm_props = self.update_atmosphere(h_k)
            P, rho, T = atm_props

            self.P_vals.append(P)
            self.rho_vals.append(rho)
            self.T_vals.append(T)

            a_k = np.sqrt(1.4 * 287 * T)  # Calculate the speed of sound (a_k) at the given temperature using the ideal gas law.
            Ma_k = u_k / a_k
            CD_k = 0  # placeholder for drag coefficient
            self.CD_vals.append(CD_k)
            self.Ma_vals.append(Ma_k)

            dudt = self.rocket.F / m_k - (CD_k / (2 * m_k)) * rho * A * pow(u_k, 2) - self.g  # Compute the rate of change of velocity (dudt) considering thrust, drag, and gravity.
            u = u_k + dudt * self.dt
            h = h_k + u_k * self.dt + 0.5 * dudt * pow(self.dt, 2)
            m = m_k - self.rocket.m_dot * self.dt

            self.u_vals.append(u)
            self.h_vals.append(h)
            self.m_vals.append(m)

        return {
            'time': self.t_vals,
            'velocity': self.u_vals,
            'altitude': self.h_vals,
            'mass': self.m_vals,
            'pressure': self.P_vals,
            'density': self.rho_vals,
            'temperature': self.T_vals,
            'Mach': self.Ma_vals
        }

    def run_2D(self):
        # Run a 2D simulation considering both vertical and horizontal dynamics.
        self.initialize_values()
        A = self.rocket.A
        theta_o = self.rocket.theta_i * (math.pi / 180)  # Convert initial theta from degrees to radians

        # Simulate rocket physics up until the point of motor burn out.
        for t in self.t_vals:
            # Set conditions from previous time step. Use initial conditions if needed.
            if len(self.u_vals) == 0:
                u_k, h_k, r_k, m_k, theta_k = 0, 0, 0, self.rocket.m_i, theta_o
                ux_k, uy_k = 0, 0  # initial velocity components
            else:
                u_k = self.u_vals[-1]
                h_k = self.h_vals[-1]
                r_k = self.r_vals[-1]
                m_k = self.m_vals[-1]
                theta_k = self.theta_vals[-1]
                ux_k = self.ux_vals[-1]
                uy_k = self.uy_vals[-1]

            # Update atmospheric properties.
            atm_props = self.update_atmosphere(h_k)
            P, rho, T = atm_props
            self.P_vals.append(P)
            self.rho_vals.append(rho)
            self.T_vals.append(T)

            # Update aerodynamic properties.
            a_k = np.sqrt(1.4 * 287 * T)
            Ma_k = u_k / a_k
            CD_k = self.rocket.get_CD(Ma_k)
            # CD_k = 0
            self.CD_vals.append(CD_k)
            self.Ma_vals.append(Ma_k)

            # Update rocket dynamics.
            dudt = self.rocket.F / m_k - (CD_k / (2 * m_k)) * rho * A * pow(u_k, 2) - self.g  # Compute the rate of change of velocity (dudt) considering thrust, drag, and gravity. * math.cos(theta_k)
            dvdt = -self.g * math.sin(theta_k)

            if u_k == 0:  # Avoid division by zero in the calculation of the change in flight angle.
                dthetadt = 0
            else:
                dthetadt = dvdt / u_k

            u = u_k + dudt * self.dt
            v = uy_k + dvdt * self.dt
            theta = theta_k + dthetadt * self.dt  # Update the flight angle (theta) based on its rate of change over the time step.

            ux = u * math.cos(theta)
            uy = u * math.sin(theta)

            # Update rocket location.
            h = h_k + uy * self.dt
            r = r_k + ux * self.dt

            m = m_k - self.rocket.m_dot * self.dt

            self.u_vals.append(u)
            self.ux_vals.append(ux)
            self.uy_vals.append(uy)
            self.h_vals.append(h)
            self.r_vals.append(r)
            self.m_vals.append(m)
            self.theta_vals.append(theta)

        while self.h_vals[-1] > 0:
            # Set initial conditions.
            u_k = self.u_vals[-1]
            h_k = self.h_vals[-1]
            r_k = self.r_vals[-1]
            m_k = self.m_vals[-1]
            theta_k = self.theta_vals[-1]
            ux_k = self.ux_vals[-1]
            uy_k = self.uy_vals[-1]

            # Update atmospheric properties.
            atm_props = self.update_atmosphere(h_k)
            P, rho, T = atm_props
            self.P_vals.append(P)
            self.rho_vals.append(rho)
            self.T_vals.append(T)

            # Update aerodynamic properties.
            a_k = np.sqrt(1.4 * 287 * T)
            Ma_k = u_k / a_k
            CD_k = self.rocket.get_CD(Ma_k)
            # CD_k = 0
            CL_k = self.rocket.get_CL(Ma_k)
            # CL_k = 0
            self.CD_vals.append(CD_k)
            self.Ma_vals.append(Ma_k)

            # Update rocket dynamics.
            dudt =  -1* (CD_k / (2 * m_k)) * rho * A * pow(u_k, 2) - math.sin(theta_k)*self.g  # Compute the rate of change of velocity (dudt) considering thrust, drag, and gravity. * math.cos(theta_k)
            dvdt = (CL_k / (2 * m_k)) * rho * A * pow(u_k, 2)-self.g * math.cos(theta_k)

            if u_k == 0:  # Avoid division by zero in the calculation of the change in flight angle.
                dthetadt = 0
            else:
                dthetadt = dvdt / u_k

            u = u_k + dudt * self.dt
            v = uy_k + dvdt * self.dt
            theta = theta_k + dthetadt * self.dt  # Update the flight angle (theta) based on its rate of change over the time step.

            ux = u * math.cos(theta)
            uy = u * math.sin(theta)

            # Update rocket location.
            h = h_k + uy * self.dt
            r = r_k + ux * self.dt

            self.t_vals = np.append(self.t_vals, self.t_vals[-1] + self.dt)
            self.u_vals.append(u)
            self.ux_vals.append(ux)
            self.uy_vals.append(uy)
            self.h_vals.append(h)
            self.r_vals.append(r)
            self.m_vals.append(m)
            self.theta_vals.append(theta)
            print(h)


        # Added print statement for final flight angle and velocity
        p=2
        print('====================== Burn Out Time Conditions ===========================')
        print(f'Flight path angle: {round(self.theta_vals[-1],p)} radians, Velocity: {round(self.u_vals[-1],p)} m/s, Altitude: {round(self.h_vals[-1],p)} m')

        return {
            't': self.t_vals,
            'u': self.u_vals,
            'ux': self.ux_vals,
            'uy': self.uy_vals,
            'h': self.h_vals,
            'r': self.r_vals,
            'mass': self.m_vals,
            'theta': self.theta_vals,
            'P': self.P_vals,
            'rho': self.rho_vals,
            'T': self.T_vals,
            'Ma': self.Ma_vals
        }

