import os
import numpy as np
import pandas as pd
from pvhub import *

# Print out all the models one can use
print("Available model classes:")
print(get_models())
print("")

# Example of querying the 2M++_SDSS model for a single object, with object 2013dy from example.csv
print("Querying a single model for a single object:")

model = TwoMPP_SDSS(verbose=True)
test_RA = 334.6
test_Dec = 40.6
test_zcmb = 0.0029
pv = model.calculate_pv(test_RA, test_Dec, test_zcmb)
print(f"PV of object at (RA, Dec, zcmb) = ({test_RA}, {test_Dec}, {test_zcmb}): {pv}\n")

# Turning off extrapolation beyond 2M++, z>0.067
print("Querying a different model for objects at different redshift, this time with no extrapolation:")

model = TwoMPP_SDSS_6dF(verbose=True)
test_zs = [0.05, 0.06, 0.07]
pv = model.calculate_pv([334.6] * 3, [40.6] * 3, test_zs, extrapolation=False)
for z, p in zip(test_zs, pv):
    print(f"z={z}, vpec={p}")
print("")

# Example of data arrays
print("Querying the same model for an array of other object:")

# First line just works out the relative location of example.py, so we can always find example.csv
current_file = os.path.dirname(inspect.stack()[0][1])
inp = pd.read_csv(os.path.normpath(current_file) + "/example.csv")
pv = model.calculate_pv(inp["RA_host"], inp["Dec_host"], inp["zcmb"])
for sn, p in zip(inp["SNID"], pv):
    print(f"PV of object {sn} = {p}")
print("")

# Finally, loop over all the models and return the predicted peculiar velocities for the example file.
# Store them in a dataframe with the model name as columns, and each SNe as a different row.
print("Building a table of results for all objects and all models:")
pvs = {}
for model_class in get_models():
    model = model_class()
    pvs[model.name] = model.calculate_pv(inp["RA_host"], inp["Dec_host"], inp["zcmb"])
pvs = pd.DataFrame.from_dict(pvs, orient="index", columns=inp["SNID"]).transpose()
print(pvs)

print("All supernova PVs for the 2MRS model:")
print(pvs["2MRS_redshift"].to_numpy())

print("All model PVs for SNe 2017cbv:")
print(pvs.loc["2017cbv"])