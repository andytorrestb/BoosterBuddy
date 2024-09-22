import math
import pandas as pd

def interpolate(h, x2, x1, prop):
    dxdh = (x2[prop]- x1[prop]) / (x2['h'] - x1['h'])

    x = x1[prop] + (h - x1['h']) * dxdh

    return float(x)

class Atmosphere:
    def __init__(self):
        return

    def get_atm_props(self, h):
 
        atm = self.atm

        h_vals = atm['h']
        
        # for row in atm.iterrows():
        #     print(row[0])
        #     print(row[1])
        #     print()
        #     if row[1]['h'] < h:
        #         pass
        
        # Find rows of data used for interpolation.
        for i in range(len(atm)):
            # Save current row for reference.
            row  = atm.iloc[i]

            # If h is found, return data for that row. 
            if row['h'] == h:
                return row.to_dict()

            # If we pass h, save data for interpolation.
            if row['h'] > h:
                x2 = atm.iloc[i]
                x1 = atm.iloc[i-1]
                break
            
            # Rocket has left the atmopsphere.
            if i == 15:
                atm_sample = {
                    'T': 1000,
                    'rho': 0,
                    'PPSL': 0,
                    'h': h
                }
                return atm_sample

        # Interpolate T, rho, and P
        props = ['T', 'rho', 'PPSL']
        atm_sample = {}
        for prop in props:
            atm_sample[prop] = interpolate(h, x2, x1, prop)

        # Save the h value used.
        atm_sample['h'] = h

        return atm_sample


    def set_atm_props(self):

        # Properties of Earth's Standard Atmosphere. Appendex 2.
        h = [0, 1000, 3e3, 5e3, 10e3, 25e3, 50e3, 75e3, 100e3, 130e3, 160e3, 200e3, 300e3, 400e3, 600e3, 1000e3]
        T = [288, 282, 269, 255, 223, 222, 271, 206, 195, 469, 696, 846, 976, 996, 999, 1000]
        PPSL = [1.0, 0.887, 0.669, 0.5331, 0.26151, 0.02516, 7.8735e-4, 2.0408e-5, 3.1593e-7, 1.2341e-8, 2.9997e-9, 8.3628e-10, 8.6557e-11, 1.4328e-11, 8.1056e-13, 7.4155e-14]
        rho = [1.225, 1.1117, 9.0912e-1, 7.6312e-1, 4.1351e-1, 4.0084e-2, 1.0269e-3, 3.4861e-5, 5.604e-7, 8.152e-9, 1.233e-9, 2.541e-10, 1.916e-11, 2.803e-12, 2.137e-13, 3.561e-15]

        atm = {'h': h, 'T':T, 'PPSL': PPSL, 'rho': rho}
        atm = pd.DataFrame(atm)

        self.atm = atm

        return
