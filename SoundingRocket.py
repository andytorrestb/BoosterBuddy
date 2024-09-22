import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from Rocket import Rocket

class SoundingRocket(Rocket):

    def plot_altitude(self):

        t_b = self.t_b
        n_t = self.n_t

        t_vals = np.linspace(0, t_b, n_t)

        h_vals = []

        for t in t_vals:
            print(t)
            h_vals.append(self.get_hb_altitude(t))

        # h_vals.reverse()

        plt.plot(t_vals, h_vals)
        plt.savefig('sounding_h_vs_t.png')
        plt.clf()

        return h_vals, t_vals

    def plot_velocity(self):

        t_b = self.t_b
        n_t = self.n_t

        t_vals = np.linspace(0, t_b, n_t)

        u_vals = []

        for t in t_vals:
            print(t)
            u_vals.append(self.get_velocity(t))

        # u_vals.reverse()

        plt.plot(t_vals, u_vals)
        plt.savefig('sounding_u_vs_t.png')
        plt.clf()

        return u_vals, t_vals


    def get_altitude(self, t):
        g = 9.81
        u_eq = self.I_sp*g
        R = 1/((self.m_i - self.m_p) / self.m_i )
        t_b = self.t_b

        k = (u_eq*t_b*R)/(R-1)
        h = (1 - (1 - 1/R)*(t/t_b))
        h = h*math.log(1 - (1 - 1/R)*(t/t_b))
        h = h - (1 - (1 - 1/R)*(t/t_b))
        h = k*h - 0.5*g*t**2
        return h

    def get_hb_altitude(self, t):
        g = 9.81
        u_eq = self.I_sp*g
        R = 1/((self.m_i - self.m_p) / self.m_i )
        t_b = self.t_b
        

        # h_b = -Veq*tb*(log(R)./(R-1))+(Veq.*tb)-(0.5*g0*(tb^2));
 

        h = -1*u_eq*t*(math.log(R) / (R-1)) + u_eq*t - 0.5*g*t**2

        return h



    def get_velocity(self, t):
        g = 9.81
        u_eq = self.I_sp*g
        R = 1/((self.m_i - self.m_p) / self.m_i )
        t_b = self.t_b

        u = -1*u_eq*math.log(1 - (1 - 1/R)*(t/t_b)) - g*t
        return u

    # def
