import math
import numpy as np
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
c = 299792.458
#----------------------------------------------------------------#
modelflag = int(input("Which model should I use? 0. 2MPP-SDSS (Said+22), 1. 2MPP-SDSS-6dF (Said+22), 2. 2MRS (Lilow & Nusser 2021), 3. 2MPP (Carrick+2015) ---->"))
#modelflag = 1
#----------------------------------------------------------------#
if modelflag == 0:
    usemodel = '2MPP_SDSS.txt'
    usemodel2 = '2MPP_SDSS_out.txt'
if modelflag == 1:
    usemodel = '2MPP_SDSS_6dF.txt'
    usemodel2 = '2MPP_SDSS_6dF_out.txt'
if modelflag == 2:
    usemodel = '2MRS_redshift.txt'
    usemodel2 = '2MRS_redshift_out.txt'
if modelflag == 3:
    usemodel = '2MPP_redshift.txt'
    usemodel2 = '2MPP_redshift_out.txt'
if modelflag > 3:
    print("Unknown Model")
    raise ValueError
#-------------------------------------------------------------------------------------------#
sgx_tm = []; sgy_tm = []; sgz_tm = []; vpred_pr = []
vpred_x = []; vpred_y = []; vpred_z = []; vproj = []
for line in open('./data/'+usemodel,"r"):
    if line[0] != '#':
        sgx_tm.append(float(line.split()[0]))
        sgy_tm.append(float(line.split()[1]))
        sgz_tm.append(float(line.split()[2]))
        vproj.append(float(line.split()[3]))
        vpred_x.append(float(line.split()[4]))
        vpred_y.append(float(line.split()[5]))
        vpred_z.append(float(line.split()[6]))
        
dmin = -20000.
dmax = 20000.
nbins = 129
bsz = ((dmax-dmin)/float(nbins-1.))
#----------------------------------------------------------------------------------------------------------#
#--------------Model beyond Vext---------------------------------------------------------------------------#
zcmb_m=[];vx=[];vy=[];vz=[]
for line in open('./data/'+usemodel2,"r"):
    if line[0]!='#':
        zcmb_m.append(float(line.split()[0]))
        vx.append(float(line.split()[7]))
        vy.append(float(line.split()[8]))
        vz.append(float(line.split()[9]))
#--------------------------------------------------------------------#
def calculate_pv(mastername,mastersource,RA,DEC,z_cmb_in,extrapolation="Yes"):
    ra0 = RA
    dec0 = DEC
    cz  = c*z_cmb_in
    zcmb = z_cmb_in
    ccc = SkyCoord(ra0*u.degree, dec0*u.degree, distance=cz*u.km/u.s, frame='icrs')
    sgc = ccc.transform_to('supergalactic')
    sgc.representation_type = 'cartesian'
    if cz < 0.0665218237200617*c and dmin < sgc.sgx.value < dmax and dmin < sgc.sgy.value < dmax and dmin < sgc.sgz.value < dmax:
        xbin = int(round(((sgc.sgx.value - dmin)/bsz),0))
        ybin = int(round(((sgc.sgy.value - dmin)/bsz),0))
        zbin = int(round(((sgc.sgz.value - dmin)/bsz),0))
        binindex = int(xbin)*nbins*nbins + int(ybin)*nbins + int(zbin)
        pv = vproj[binindex]
        v_x = vpred_x[binindex]
        v_y = vpred_y[binindex]
        v_z = vpred_z[binindex]
#            if pv != -10000.:
        return 'pv=',int(round(pv,0)),'within reconstruction:','True'
    else:
        if extrapolation == 'No':
            return 'pv=',9999,'within reconstruction:','False'
        if extrapolation == 'Yes':
            k = np.searchsorted(zcmb_m, zcmb)
            vx01 = vx[k]
            vy01 = vy[k]
            vz01 = vz[k]
            vpred = ((sgc.sgx.value*vx01)+(sgc.sgy.value*vy01)+(sgc.sgz.value*vz01))/(np.sqrt(sgc.sgx.value**2.+sgc.sgy.value**2.+sgc.sgz.value**2.))
            return 'pv=',int(round(vpred,0)),'within reconstruction:','False'
        else:
            print("Unknown extrapolation type")
            raise ValueError
#-------------------------------------------------------------------------------------------------------------------#
