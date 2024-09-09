import matplotlib.pyplot as plt
import Rocket

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

  # initial flight path angle.
  'theta_i': 70 # degrees
}

r = Rocket.Rocket(rocket_params)
r.set_rocket_thrust()

t_vals, u_vals = r.sim_flight()

plt.plot(t_vals, u_vals)
plt.savefig('test.png')