import matplotlib.pyplot as plt
import Rocket, FlightSimulator, Atmosphere
import SoundingRocket

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

# Plot sounding rocket data
sr = SoundingRocket.SoundingRocket(rocket_params)
h_sounding, t_sounding = sr.plot_altitude()
u_sounding, t_sounding = sr.plot_velocity()

# Plot rocket solution.
r = Rocket.Rocket(rocket_params)
r.set_CD_data()
r.set_CL_data()

# Set atmosphere data.
atm = Atmosphere.Atmosphere()
atm.set_atm_props()


# # Instantiate flight simulator.
flight_sim = FlightSimulator.FlightSimulator(r, atm)

results = flight_sim.run_1D()

plt.plot(results['t'], results['u'], label='rocket-solution')
plt.plot(t_sounding, u_sounding, label='sounding-rocket')
plt.title('Burnout Velocity Comparison')
plt.ylabel('Velocity (U)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('u_compare.png')
plt.clf()

plt.plot(results['t'], results['h'], label='rocket-solution')
plt.plot(t_sounding, h_sounding, label='sounding-rocket')
plt.title('Burnout Altitude Comparison')
plt.ylabel('Altitude (h)')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('h_compare.png')
plt.clf()
