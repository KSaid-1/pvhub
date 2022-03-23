import os
import inspect
import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const
from abc import ABC
from astropy.coordinates import SkyCoord


class Recon(ABC):
    """Abstract base class for the reconstruction data included in pvhub

    Meant to only be accessed through the various listed reconstruction subclasses. The attributes below are inherited
    by these subclasses

    Attributes
    ----------
    """

    def __init__(self):

        current_file = os.path.dirname(inspect.stack()[0][1])
        self.data_location = os.path.normpath(current_file + "/data/")

        # The properties of the grid on which the reconstructions are stored.
        self.dmin = -20000.0
        self.dmax = 20000.0
        self.nbins = 129
        self.bsz = (self.dmax - self.dmin) / float(self.nbins - 1.0)

        self.name = None
        self.vmodel, self.vextmodel = None, None

    def _load_data(self, modelfile, extfile):
        """Loads the relevant reconstruction model files and stores the attributes

        Args
        ----
            modelfile : string
                The filepath for the file containing the reconstructed field as a function of
                cartesian coordinates
            extfile : string
                The filepath for the file containing the external field prediction as a function of redshift

        Raises
        ------
        """
        # The reconstruction model
        self.vmodel = pd.read_csv(self.data_location + modelfile, delim_whitespace=True)

        # The model beyond Vext
        self.vextmodel = pd.read_csv(self.data_location + extfile, delim_whitespace=True)

        return

    def calculate_pv(self, ra, dec, zcmb, extrapolation=True):
        """Interpolates/extrapolates the reconstruction for a set of ra, dec and CMB-frame redshifts

        Args
        ----
            ra : float, np.ndarray
                A value or array of values for the Right Ascention to be queried.
            dec : float, np.ndarray
                A value or array of values for the Declination to be queried.
            z_cmb_in:  float, np.ndarray
                A value or array of values for the CMB-frame redshifts to be queried.
            extrapolation: bool
                Whether or not to extrapolate outside of the reconstructed field based on a smooth model of the
                velocity as a function of redshift. Default = True

        Returns
        -------
            pv : np.ndarray
                An array of peculiar velocities for each of the coordinates passed in

        Raises
        ------
        """

        ra, dec, zcmb = np.asarray(ra), np.asarray(dec), np.asarray(zcmb)
        ccc = SkyCoord(
            ra * u.degree,
            dec * u.degree,
            distance=const.c.value / 1000.0 * zcmb * u.km / u.s,
            frame="icrs",
        )

        sgc = ccc.transform_to("supergalactic")
        sgc.representation_type = "cartesian"

        xbin = np.round(((sgc.sgx.value - self.dmin) / self.bsz), 0)
        ybin = np.round(((sgc.sgy.value - self.dmin) / self.bsz), 0)
        zbin = np.round(((sgc.sgz.value - self.dmin) / self.bsz), 0)

        # calculate bin index even if coords outside 2M++
        binindex = xbin.astype(int) * self.nbins * self.nbins + ybin.astype(int) * self.nbins + zbin.astype(int)

        # set indices outside 2M++ to 0
        try:
            binindex[np.where((binindex < 0) | (binindex >= len(self.vmodel["vproj_2MPP"])))] = 0
        except TypeError:  # For single input
            pass

        k = np.searchsorted(self.vextmodel["z"], zcmb)  # calculate bin index even if coords inside 2M++

        in2MPP = (
            (zcmb < 0.0665218237200617)  # precise redshift of 2M++ boundary
            & ((self.dmin < sgc.sgx.value) & (sgc.sgx.value < self.dmax))
            & ((self.dmin < sgc.sgy.value) & (sgc.sgy.value < self.dmax))
            & ((self.dmin < sgc.sgz.value) & (sgc.sgz.value < self.dmax))
        )

        if extrapolation:
            r = np.sqrt(np.sum(np.square([sgc.sgx.value, sgc.sgy.value, sgc.sgz.value]), axis=0))
            vdot = (
                sgc.sgx.value * self.vextmodel["Vsgx"].loc[k]
                + sgc.sgy.value * self.vextmodel["Vsgy"].loc[k]
                + sgc.sgz.value * self.vextmodel["Vsgz"].loc[k]
            )
            vproj = np.where(
                self.vmodel["vproj_2MPP"].loc[binindex].isna(), 0.0, self.vmodel["vproj_2MPP"].loc[binindex].isna()
            )
            pv = np.round(np.where(in2MPP, self.vmodel["vproj_2MPP"].loc[binindex], vdot / r), 0)
        else:
            pv = np.where(in2MPP, np.round(self.vmodel["vproj_2MPP"].loc[binindex], 0), np.nan)

        return pv


class TwoMPP_SDSS(Recon):
    def __init__(self):
        super().__init__()

        self.name = "2M++_SDSS"
        model = "/2MPP_SDSS.txt"
        model_ext = "/2MPP_SDSS_out.txt"

        self._load_data(model, model_ext)


class TwoMPP_SDSS_6dF(Recon):
    def __init__(self):
        super().__init__()

        self.name = "2M++_SDSS_6dF"
        model = "/2MPP_SDSS_6dF.txt"
        model_ext = "/2MPP_SDSS_6dF_out.txt"

        self._load_data(model, model_ext)


class TwoMRS_redshift(Recon):
    def __init__(self):
        super().__init__()

        self.name = "2MRS"
        model = "/2MRS_redshift.txt"
        model_ext = "/2MRS_redshift_out.txt"

        self._load_data(model, model_ext)


class TwoMPP_redshift(Recon):
    def __init__(self):
        super().__init__()

        self.name = "2M++"
        model = "/2MPP_redshift.txt"
        model_ext = "/2MPP_redshift_out.txt"

        self._load_data(model, model_ext)


def get_concrete(baseclass):
    classes = baseclass.__subclasses__()
    for c in classes:
        classes += c.__subclasses__()
    final_classes = [c for c in classes if ABC not in c.__bases__]
    return final_classes


if __name__ == "__main__":

    from pvhub import *

    # Get a list of all models
    models = [c for c in get_concrete(Recon)]

    # Read in some test data
    inp = pd.read_csv("./inputs/example.csv")

    # Loop over all the models and return the predicted peculiar velocities.
    # Store them in a dataframe with the model name as rows, and each SNe as columns.
    pvs = {}
    for c in models:
        model = c()
        pvs[model.name] = model.calculate_pv(inp["RA_host"], inp["Dec_host"], inp["zcmb"])
    pvs = pd.DataFrame.from_dict(pvs, orient="index", columns=inp["SNID"])

    # Example of querying for a single object
    test_RA = 334.6
    test_Dec = 40.6
    test_zcmb = 0.0029
    pv = model.calculate_pv(test_RA, test_Dec, test_zcmb)
    print(f"PV of object at (RA, Dec, zcmb) = ({test_RA}, {test_Dec}, {test_zcmb}): {pv}")

    # Turning off extrapolation beyond 2M++, z>0.067
    test_zs = [0.05, 0.06, 0.07]
    test_RA = [334.6]
    test_Dec = [40.6]
    pv = model.calculate_pv(test_RA * 3, test_Dec * 3, test_zs, extrapolation=False)
    print("No extrapolation:")
    for z, p in zip(test_zs, pv):
        print(f"z={z}, vpec={p}")
