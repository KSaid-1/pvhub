from pvhub_new import *
import pandas as pd

inp = pd.read_csv("./inputs/example.csv")
SNID = inp["SNID"]
folder01 = inp["folder"]
ra00 = inp["RA_host"]
dec00 = inp["Dec_host"]
z00 = inp["zcmb"]
pv = calculate_pv(ra00, dec00, z00, extrapolation=False)
print(pv)
