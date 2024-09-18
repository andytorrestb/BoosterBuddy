import matplotlib.pyplot as plt
import Rocket, FlightSimulator, Atmosphere

# ==============================
# || Vehicle Performance Data ||
# ==============================
rocket_params = {

  't_b': 40, # motor burn time (s)
  'm_i': 4000, # initial mass (kg)
  'm_p': 2500, # propellant mass (kg)
  'I_sp': 260,  # specific impulse (s)

  'd_ref': 1.25, # Reference diameter (m) 

  # Initial flight path angle.
  'theta_i': 70 # degrees

}

# Set rocket data.
r = Rocket.Rocket(rocket_params)
r.get_CD_data()
r.get_CL_data()

# Set atmosphere data.
atm = Atmosphere.Atmosphere()
atm.set_atm_props()
atm.get_atm_props(500000)
input()

# Instantiate flight simulator.
flight_sim = FlightSimulator.FlightSimulator(r, atm)

t_vals, u_vals, h_vals, rho_vals, P_vals, T_vals = flight_sim.run()

plt.plot(t_vals, u_vals)
plt.savefig('u_vs_t.png')
plt.clf()
plt.plot(t_vals, h_vals)
plt.savefig('h_vs_t.png')

plt.clf()
plt.plot(h_vals, P_vals)
plt.savefig('P_vs_h.png')

plt.clf()
plt.plot(h_vals, T_vals)
plt.savefig('T_vs_h.png')