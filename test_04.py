import matplotlib.pyplot as plt
import Rocket, Atmosphere
import FlightSimulator

# ==============================
# || Vehicle Performance Data ||
# ==============================
rocket_params = {

  # 't_b': 10, # motor burn time (s)
  'm_i': 4000, # initial mass (kg)
  'm_p': 2500, # propellant mass (kg)
  'm_dot': 50, # propellant burn rate (kg/s)
  'I_sp': 260,  # specific impulse (s)

  'd_ref': 1.25, # Reference diameter (m) 

  # Initial flight path angle.
  'theta_i': 70, # degrees
  'n_t': 100
}

# Set rocket data.
r = Rocket.Rocket(rocket_params)
r.set_CD_data()
r.set_CL_data()

# Set atmosphere data.
atm = Atmosphere.Atmosphere()
atm.set_atm_props()
# atm.get_atm_props(50000)
# input()

# Instantiate flight simulator.
flight_sim = FlightSimulator.FlightSimulator(r, atm)
flight_sim.set_dt()

results = flight_sim.run_2D()

plt.plot(results['t'], results['u'])
plt.savefig('u_vs_t.png')
plt.clf()
plt.plot(results['t'], results['h'])
plt.savefig('h_vs_t.png')

plt.clf()
plt.plot(results['h'], results['P'])
plt.savefig('P_vs_h.png')

plt.clf()
plt.plot(results['h'], results['T'])
plt.savefig('T_vs_h.png')

plt.clf()
plt.plot(results['h'], results['rho'])
plt.savefig('rho_vs_h.png')

plt.clf()
plt.plot(results['t'], results['Ma'])
plt.savefig('Ma_vs_t.png')

plt.clf()
plt.plot(results['r'], results['h'])
plt.savefig('h_vs_r.png')

plt.clf()
plt.plot(results['t'], results['theta'])
plt.savefig('theta_vs_t.png')