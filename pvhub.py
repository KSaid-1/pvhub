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

    def __init__(self, name=None, modelfile=None, extfile=None, verbose=False):

        current_file = os.path.dirname(inspect.stack()[0][1])
        self.data_location = os.path.normpath(current_file + "/data") + "/"

        self.name = name
        self.verbose = verbose

        # The properties of the grid on which the reconstructions are stored.
        self.dmin = -20000.0
        self.dmax = 20000.0
        self.nbins = 129
        self.bsz = (self.dmax - self.dmin) / float(self.nbins - 1.0)

        if modelfile is not None and extfile is not None:
            self._load_data(modelfile, extfile)
        else:
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
        try:

            # The reconstruction model
            self.vmodel = pd.read_csv(self.data_location + modelfile)

            # The model beyond Vext
            self.vextmodel = pd.read_csv(self.data_location + extfile, delim_whitespace=True)

            if self.verbose:
                print("Loaded model " + self.name)

        except:

            print("ERROR: Unable to load model " + self.name)
            print("Please check files " + modelfile + " and " + extfile + " are in " + self.data_location)

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
            # Note that cz is *not* a distance, but treating it as such is self-consistent with 2M++
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
            (const.c.value * zcmb < 19942)  # precise redshift of 2M++ boundary
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
            pv = np.where(
                in2MPP,
                self.vmodel["vproj_2MPP"].loc[binindex],
                vdot / r,
            )
            # 2M++ is not completely spherical, so extrapolation must be extended inwards in some
            # parts of the sky. The PV of the centre of the reconstruction is undefined, so we only
            # use the extrapolation above a redshift securely inside 2M++ but outside the central cell
            pv = np.where((np.isnan(pv)) & (zcmb > 0.01), vdot / r, pv)
            pv = np.round(pv, 0)
        else:
            pv = np.where(in2MPP, np.round(self.vmodel["vproj_2MPP"].loc[binindex], 0), np.nan)

        return pv


class TwoMPP_SDSS(Recon):
    def __init__(self, verbose=False):

        name = "2M++_SDSS"
        model = "2MPP_SDSS.txt"
        model_ext = "2MPP_SDSS_out.txt"
        super().__init__(name = name, modelfile = model, extfile= model_ext, verbose=verbose)

class TwoMPP_SDSS_6dF(Recon):
    def __init__(self, verbose=False):

        name = "2M++_SDSS_6dF"
        model = "2MPP_SDSS_6dF.txt"
        model_ext = "2MPP_SDSS_6dF_out.txt"
        super().__init__(name = name, modelfile = model, extfile= model_ext, verbose=verbose)


class TwoMRS_redshift(Recon):
    def __init__(self, verbose=False):

        name = "2MRS_redshift"
        model = "2MRS_redshift.txt"
        model_ext = "2MRS_redshift_out.txt"
        super().__init__(name = name, modelfile = model, extfile= model_ext, verbose=verbose)

class TwoMPP_redshift(Recon):
    def __init__(self, verbose=False):

        name = "2M++_redshift"
        model = "2MPP_redshift.txt"
        model_ext = "2MPP_redshift_out.txt"
        super().__init__(name = name, modelfile = model, extfile= model_ext, verbose=verbose)


def get_models():
    """Utility function to return all the reconstruction model subclasses

    Args
    ----
        None

    Returns
    -------
        final_classes: list
            A list of all the reconstruction model subclasses
    """

    classes = Recon.__subclasses__()
    for c in classes:
        classes += c.__subclasses__()
    final_classes = [c for c in classes if ABC not in c.__bases__]

    return final_classes


if __name__ == "__main__":

    from pvhub import *

    # Get a list of all models and check they can be read in and used
    # This is just a unit test. See examples.py for usage examples
    print("Checking models load and run:")
    for model in get_models():
        model(verbose=True).calculate_pv(334.6, 40.6, 0.0029)