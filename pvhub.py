import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

c = 299792.458
modelflag = 0

if modelflag == 0:
    model = "2MPP_SDSS.txt"
    model_ext = "2MPP_SDSS_out.txt"
elif modelflag == 1:
    model = "2MPP_SDSS_6dF.txt"
    model_ext = "2MPP_SDSS_6dF_out.txt"
elif modelflag == 2:
    model = "2MRS_redshift.txt"
    model_ext = "2MRS_redshift_out.txt"
elif modelflag == 3:
    model = "2MPP_redshift.txt"
    model_ext = "2MPP_redshift_out.txt"
else:
    print("Unknown Model")
    raise ValueError

inp = pd.read_csv(f"./data/{model}", delim_whitespace=True)
sgx_tm = inp["sgx"]
sgy_tm = inp["sgy"]
sgz_tm = inp["sgz"]
vpred_x = inp["vpred_x"]
vpred_y = inp["vpred_y"]
vpred_z = inp["vpred_z"]
vproj = inp["vproj_2MPP"]

dmin = -20000.0
dmax = 20000.0
nbins = 129
bsz = (dmax - dmin) / float(nbins - 1.0)

# Model beyond Vext
inp = pd.read_csv(f"./data/{model_ext}", delim_whitespace=True)
zcmb_m = inp["z"]
vx = inp["Vsgx"]
vy = inp["Vsgy"]
vz = inp["Vsgz"]


def calculate_pv(RA, DEC, z_cmb_in, extrapolation=True):
    cz = c * np.array(z_cmb_in)
    zcmb = z_cmb_in
    ccc = SkyCoord(
        RA * u.degree, DEC * u.degree, distance=cz * u.km / u.s, frame="icrs"
    )
    sgc = ccc.transform_to("supergalactic")
    sgc.representation_type = "cartesian"
    xbin = np.round(((sgc.sgx.value - dmin) / bsz), 0)
    ybin = np.round(((sgc.sgy.value - dmin) / bsz), 0)
    zbin = np.round(((sgc.sgz.value - dmin) / bsz), 0)
    binindex = (
        xbin.astype(int) * nbins * nbins + ybin.astype(int) * nbins + zbin.astype(int)
    )  # calculate bin index even if coords outside 2M++

    try:
        binindex[
            np.where((binindex < 0) | (binindex >= len(vproj)))
        ] = 0  # set indices outside 2M++ to 0
    except TypeError:  # For single input
        pass

    k = np.searchsorted(zcmb_m, zcmb)  # calculate bin index even if coords inside 2M++

    in2MPP = (
        (cz < 0.0665218237200617 * c)  # precise redshift of 2M++ boundary
        & ((dmin < sgc.sgx.value) & (sgc.sgx.value < dmax))
        & ((dmin < sgc.sgy.value) & (sgc.sgy.value < dmax))
        & ((dmin < sgc.sgz.value) & (sgc.sgz.value < dmax))
    )
    if extrapolation:
        pv = np.where(
            in2MPP,
            vproj.loc[binindex],
            (
                (sgc.sgx.value * vx[k])
                + (sgc.sgy.value * vy[k])
                + (sgc.sgz.value * vz[k])
            )
            / (
                np.sqrt(
                    np.sum(
                        np.square([sgc.sgx.value, sgc.sgy.value, sgc.sgz.value]), axis=0
                    )
                )
            ),
        )
        pv = np.round(pv, 0)
    else:
        pv = np.where(in2MPP, np.round(vproj.loc[binindex], 0), np.nan)
    return pv
