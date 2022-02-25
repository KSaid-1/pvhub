from pvhub import *
import pandas as pd

# Example of data arrays
inp = pd.read_csv("./inputs/example.csv")
SNID = inp["SNID"]
ra = inp["RA_host"]
dec = inp["Dec_host"]
zcmb = inp["zcmb"]
pv = calculate_pv(ra, dec, zcmb)
for sn, p in zip(SNID, pv):
    print(f"PV of object {sn} = {p}")

# Example of querying for a single object
test_RA = 334.6
test_Dec = 40.6
test_zcmb = 0.0029
pv = calculate_pv(test_RA, test_Dec, test_zcmb)
print(f"PV of object at (RA, Dec, zcmb) = ({test_RA}, {test_Dec}, {test_zcmb}): {pv}")

# Turning off extrapolation beyond 2M++, z>0.067
test_zs = [0.05, 0.06, 0.07]
test_RA = [334.6]
test_Dec = [40.6]
pv = calculate_pv(test_RA * 3, test_Dec * 3, test_zs, extrapolation=False)
print("No extrapolation:")
for z, p in zip(test_zs, pv):
    print(f"z={z}, vpec={p}")
