# Some example code based around GW170817 and its host galaxy NGC4993 that uses pvhub.
# Includes:
# Routines to return 1) the luminosity distance to this galaxy given a value for H0 and peculiar velocity
#                    2) A sample of luminosity distances to this galaxy given a value for H0 and drawing PVs from the reconstruction
# 					 3) A likelihood for GW170817 given a luminosity distance
# 					 4) Constraints on H0 and the PV of NGC4993 given GW170817 and a reconstruction

from pvhub import *
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import gaussian_kde
from scipy.interpolate import CubicSpline
from chainconsumer import ChainConsumer, Chain

LightSpeed = const.c.value / 1000.0


# Heliocentric to CMB frame conversion routine. Adapted from Anthony Carr's code/work.
def perform_corr(z_helio, RA, Dec):
    """
    A function to perform heliocentric corrections on redshifts.
    Inputs:
        z_helio: float, input heliocentric
        RA: float, equatorial right ascension
        Dec: float, equatorial declination

    Outputs:
        z_CMB: float

    Notes:
      Co-ords of North Galactic Pole (ICRS): RA = 192.729 ± 0.035 deg, Dec = 27.084 ± 0.023 deg (https://doi.org/10.1093/mnras/stw2772)
      Co-ords of Galactic Centre (ICRS): RA = 17h45m40.0409s, Dec = −29d00m28.118s (see above reference)
                                         RA = 266.41683708 deg, Dec = -29.00781056 deg
      Ascending node of the galactic plane = arccos(sin(Dec_GC)*cos(Dec_NGP)-cos(Dec_GC)*sin(Dec_NGP)*cos(RA_NGP-RA_GC))
                                           = 122.92828126730255 = l_0
      Transform CMB dipole from (l,b) to (RA,Dec):
          Dec = arcsin(sin(Dec_NGP)*sin(b)+cos(Dec_NGP)*cos(b)*cos(l_0-l))
              = -6.9895105228347 deg
          RA = RA_NGP + arctan((cos(b)*sin(l_0-l)) / (cos(Dec_NGP)*sin(b)-sin(Dec_NGP)*cos(b)*cos(l_0-l)))
             = 167.81671014708002 deg
      Astropy co-ordinates:
      RA_NGP_J2000 = 192.8594812065348, Dec_NGP_J2000 = 27.12825118085622, which are converted from B1950
      RA_NGP_B1950 = 192.25, Dec_NGP_B1950 = 27.4
      l_0_B1950 = 123
      l_0_J2000 = 122.9319185680026
    """

    v_Sun_Planck = 369.82  # +/- 0.11 km/s
    RA_Sun_Planck = 167.81671014708002  # deg
    Dec_Sun_Planck = -6.9895105228347  # deg
    c = 299792.458  # km/s

    rad = np.pi / 180.0

    # Vincenty formula
    alpha = np.arctan2(
        np.hypot(
            np.cos(Dec_Sun_Planck * rad) * np.sin(np.fabs(RA - RA_Sun_Planck) * rad),
            np.cos(Dec * rad) * np.sin(Dec_Sun_Planck * rad)
            - np.sin(Dec * rad) * np.cos(Dec_Sun_Planck * rad) * np.cos(np.fabs(RA - RA_Sun_Planck) * rad),
        ),
        np.sin(Dec * rad) * np.sin(Dec_Sun_Planck * rad)
        + np.cos(Dec * rad) * np.cos(Dec_Sun_Planck * rad) * np.cos(np.fabs((RA - RA_Sun_Planck)) * rad),
    )

    v_Sun_proj = v_Sun_Planck * np.cos(alpha)
    z_Sun = -v_Sun_proj / c

    z_CMB = np.where(z_helio <= 0.0, z_helio, (1.0 + z_helio) / (1.0 + z_Sun) - 1)

    return z_CMB


def generate_luminosity_distance(ra, dec, cz_helio, cz_helio_err, recon, H0, om=0.31, dist_lookup_table=None, nsamples=1000):

    pverr = 250.0
    cz_helio_samples = np.random.normal(cz_helio, cz_helio_err, nsamples)
    z_cmb = perform_corr(cz_helio_samples / LightSpeed, ra, dec)
    pv = recon.calculate_pv(ra, dec, z_cmb)
    pvsamples = np.random.normal(pv, pverr)
    z_cosmo = (1.0 + z_cmb) / (1.0 + pvsamples / LightSpeed) - 1.0

    if dist_lookup_table is not None:
        d_L = (1.0 + cz_helio / LightSpeed) * 100.0 / H0 * dist_lookup_table(z_cosmo)
    else:
        cosmo = FlatLambdaCDM(H0=H0, Om0=om)
        d_L = (1.0 + cz_helio / LightSpeed) * cosmo.comoving_distance(z_cosmo).value

    return d_L


def get_GW170817(type="high_spin"):
    """Generate a KDE of cosi vs. d_L for GW170817"""

    current_file = os.path.dirname(inspect.stack()[0][1])
    data_location = os.path.normpath(current_file + "/../data") + "/"

    if type == "vlbi":
        data = np.flip(np.array(pd.read_csv(data_location + "gw_vlbi_d_cos.dat", sep="\s+")), axis=1)
    elif type == "low_spin":
        data = np.array(pd.read_csv(data_location + "low_spin_PhenomPNRT_posterior_samples.dat.gz", sep="\s+"))
    else:
        data = np.array(pd.read_csv(data_location + "high_spin_PhenomPNRT_posterior_samples.dat.gz", sep="\s+"))

    gwpdf = gaussian_kde(np.vstack([data[:, 0], data[:, 1]]))

    return gwpdf


def log_likelihood_GW170817(
    ra,
    dec,
    cz_helio,
    cz_helio_err,
    recon_pv,
    recon_pv_err,
    gwpdf,
    p,
    zcos_lookup_table=None,
):
    """Calculate the likelihood of a luminosity distance and true peculiar velocity for NGC4993 given a reconstruction and a value for H0"""

    cosi, dl, H0, cz = p

    # The GW likelihood
    like_gw = gwpdf.logpdf([cosi, dl])[0]

    # Calculate the log likelihood of the observed group redshift
    z_cmb = perform_corr(cz / LightSpeed, ra, dec)
    like_cz = -0.5 * (cz - cz_helio) ** 2 / cz_helio_err**2

    # The PV log likelihood
    z_cos = zcos_lookup_table(dl / (1.0 + cz / LightSpeed) * H0 / 100.0)
    z_pec = (1.0 + z_cmb) / (1.0 + z_cos) - 1.0
    v_pec = LightSpeed * ((1.0 + z_pec) ** 2 - 1.0) / ((1.0 + z_pec) ** 2 + 1.0)
    like_pv = -0.5 * np.sum((v_pec - recon_pv) ** 2 / recon_pv_err**2)

    return like_gw + like_cz + like_pv


def run_nautilus(ra, dec, cz_helio, cz_helio_err, recon_pv, recon_pv_err, gwpdf, zcos_lookup_table=None):
    """Fit for H0, inclination angle, distance and true group redshift of NGC4993 given the GW170817 likelihood and the reconstruction"""

    def prior_transform(u):
        """Transforms the unit cube to the parameter space according to the priors"""

        cosi = u[0] - 1.0
        dl = 30.0 * u[1] ** 2 + 20.0
        H0 = np.exp((np.log(120.0) - np.log(20.0)) * u[2] + np.log(20.0))
        cz = 2000.0 * u[3] + 2000.0

        return cosi, dl, H0, cz

    import nautilus

    sampler = nautilus.Sampler(
        prior_transform,
        log_likelihood_GW170817,
        n_dim=4,
        n_live=500,
        likelihood_args=(ra, dec, cz_helio, cz_helio_err, recon_pv, recon_pv_err, gwpdf),
        likelihood_kwargs={"zcos_lookup_table": zcos_lookup_table},
    )
    sampler.run(verbose=True, discard_exploration=True)

    print(f"Fit finished. Effective sample size: {sampler.n_eff}")

    chain, logw, likelihood = sampler.posterior()
    logz = sampler.log_z
    weights = np.exp(logw)
    max_weight = weights.max()
    trim = max_weight / 1e5
    mask = weights > trim
    return chain[mask], weights[mask], likelihood[mask], logz


if __name__ == "__main__":

    # Set the properties of NGC4993 based on different group catalogues. There are several different values for this in the literature,
    # So I've listed the values from Table 2 of Howlett & Davis 2019. The second one is that used in the original Abbott et. al., 2017 paper.
    group_choice = 1
    ra_group = np.array([199.096, 196.973, 197.038, 196.909, 197.029, 196.999, 197.013, 196.908])
    dec_group = np.array([-17.236, -23.757, -23.862, -23.675, -23.763, -23.789, -23.781, -23.865])
    cz_helio_group = np.array([2558.0, 3016.0, 3026.0, 2935.0, 2997.0, 2995.0, 2919.0, 3028.0])
    cz_cmb_group = np.array([perform_corr(cz, ra, dec) for cz, ra, dec in zip(cz_helio_group, ra_group, dec_group)])
    cz_group_err = np.array([72.0, 72.0, 53.0, 19.0, 51.0, 32.0, 13.0, 51.0])

    # Print out all the models one can use
    print("Available reconstructions:")
    print(get_models())
    print("")

    # Generating a lookup table of redshift vs. luminosity distance in Mpc/h
    print("Generating a lookup table for comoving distances")
    cosmo = FlatLambdaCDM(H0=100.0, Om0=0.31)
    z_lookup = np.logspace(-4.0, 0.0, 1000)
    dist_lookup = cosmo.comoving_distance(z_lookup).value
    dist_lookup_table = CubicSpline(z_lookup, dist_lookup)
    zcos_lookup_table = CubicSpline(dist_lookup, z_lookup)

    # Generate a reconstruction class instance and then use this alongside some ra, dec and z_cmb values to compute luminosity distances.
    recon = TwoMPP_SDSS_6dF()
    d_L_samples = generate_luminosity_distance(
        ra_group[group_choice],
        dec_group[group_choice],
        cz_helio_group[group_choice],
        cz_group_err[group_choice],
        recon,
        70.0,
        dist_lookup_table=dist_lookup_table,
    )
    print(np.mean(d_L_samples), np.std(d_L_samples))

    # Fit for H0, inclination angle and true PV of NGC4993 given the GW170817 likelihood and the reconstruction
    c = ChainConsumer()
    gwpdf = get_GW170817()
    ra, dec, cz_helio, cz_helio_err = (
        ra_group[group_choice],
        dec_group[group_choice],
        cz_helio_group[group_choice],
        cz_group_err[group_choice],
    )
    recons = [TwoMRS_redshift(), TwoMPP_redshift(), TwoMPP_SDSS(), TwoMPP_SDSS_6dF()]
    for i, recon in enumerate(recons):
        recon_pv, recon_pv_err = recon.calculate_pv(ra, dec, cz_helio / LightSpeed), 250.0
        print(f"Running for {recon.name}: PV = {recon_pv} +/- {recon_pv_err} km/s")
        chain, weights, likelihood, logz = run_nautilus(
            ra,
            dec,
            cz_helio,
            cz_helio_err,
            recon_pv,
            recon_pv_err,
            gwpdf,
            zcos_lookup_table=zcos_lookup_table,
        )
        print(f"Finished {recon.name}")

        c.add_chain(
            Chain(
                samples=pd.DataFrame(chain, columns=[r"$cos(\iota)$", r"$d_{L}$", r"$H_{0}$", "$cz_{helio}$"]),
                weights=weights,
                log_likelihoods=likelihood,
                name=recon.name,
            )
        )

    fig = c.plotter.plot()
    plt.show()
