import matplotlib.pyplot as plt
import Rocket, FlightSimulator, Atmosphere

# ==============================
# || Vehicle Performance Data ||
# ==============================
rocket_params = {

  'm_i': 4000, # kg
  'm_p': 2500, # kg
  'I_sp': 260,  # s
  'm_dot': 50, # kg/s

  # Reference diamter for lift and drag.
  'd_ref': 1.25, # m 

  # Initial flight path angle.
  'theta_i': 70 # degrees
}

# Set rocket data.
r = Rocket.Rocket(rocket_params)
r.set_rocket_thrust()

# Set atmosphere data.
atm = Atmosphere.Atmosphere()

# Instantiate flight simulator.
flight_sim = FlightSimulator.FlightSimulator(r, atm)

t_vals, u_vals, h_vals = flight_sim.run()

plt.plot(t_vals, u_vals)
plt.savefig('u_vs_t.png')
plt.clf()
plt.plot(t_vals, h_vals)
plt.savefig('h_vs_t.png')