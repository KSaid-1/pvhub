import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd

c = 299792.458
global model_spec, map_inside, map_outside
model_spec = None


def choose_model(modelflag):
    """Choose and load one of the four models:
    0: 2M++_SDSS (Said et al. 2020, Peterson et al. 2021, Carr et al. 2021)
    1: 2M++_SDSS_6dF (Said et al. 2020)
    2: 2MRS (Lilow & Nusser 2021)
    3: 2M++ (Carrick et al. 2015)"""
    if modelflag == 0:
        modelname = "2M++_SDSS"
        model = "2MPP_SDSS.txt"
        model_ext = "2MPP_SDSS_out.txt"
    elif modelflag == 1:
        modelname = "2M++_SDSS_6dF"
        model = "2MPP_SDSS_6dF.txt"
        model_ext = "2MPP_SDSS_6dF_out.txt"
    elif modelflag == 2:
        modelname = "2MRS"
        model = "2MRS_redshift.txt"
        model_ext = "2MRS_redshift_out.txt"
    elif modelflag == 3:
        modelname = "2M++"
        model = "2MPP_redshift.txt"
        model_ext = "2MPP_redshift_out.txt"
    else:
        raise ValueError("Unknown Model")

    global model_spec, map_inside, map_outside

    print(f"Loading model {modelflag} ({modelname})")

    # Model inside reconstruction
    map_inside = pd.read_csv(f"data/{model}")

    # Model beyond reconstruction
    map_outside = pd.read_csv(f"data/{model_ext}", delim_whitespace=True)

    model_spec = modelflag
    return model, model_ext


def calculate_pv(RA, DEC, z_cmb_in, extrapolation=True):
    """Get peculiar velocities from a peculiar velocity map."""
    global model_spec, map_inside, map_inside

    if model_spec is None:
        print("No model specified with choose_model(); will load default model.")
        choose_model(0)
    else:
        print(f"Using model {model_spec}.")

    vproj = map_inside["vproj_2MPP"]

    zcmb_m = map_outside["z"]
    vx = map_outside["Vsgx"]
    vy = map_outside["Vsgy"]
    vz = map_outside["Vsgz"]

    dmin = -20000.0
    dmax = 20000.0
    nbins = 129
    bsz = (dmax - dmin) / float(nbins - 1.0)

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
        (cz < 19942)  # precise 2M++ boundary
        & ((dmin < sgc.sgx.value) & (sgc.sgx.value < dmax))
        & ((dmin < sgc.sgy.value) & (sgc.sgy.value < dmax))
        & ((dmin < sgc.sgz.value) & (sgc.sgz.value < dmax))
    )
    if extrapolation:
        r = np.sqrt(
            np.sum(np.square([sgc.sgx.value, sgc.sgy.value, sgc.sgz.value]), axis=0)
        )
        vdot = (
            (sgc.sgx.value * vx[k]) + (sgc.sgy.value * vy[k]) + (sgc.sgz.value * vz[k])
        )
        pv = np.where(
            in2MPP,
            vproj.loc[binindex],
            vdot / r,
        )
        # 2M++ is not completely spherical, so extrapolation must be extended inwards in some
        # parts of the sky. The PV of the centre of the reconstruction is undefined, so we only
        # use the extrapolation above a redshift securely inside 2M++ but outside the central cell
        pv = np.where((np.isnan(pv)) & (zcmb > 0.01), vdot / r, pv)
        pv = np.round(pv, 0)
    else:
        pv = np.where(in2MPP, np.round(vproj.loc[binindex], 0), np.nan)
    return pv
