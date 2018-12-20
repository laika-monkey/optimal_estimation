#import ssp_db, psd_functions
import numpy as np
#from pyrtmwrap.util.processWrapper import processWrapper
import os
#from pyrtmwrap.aer import lblrtm_utils
#from scipy.io.netcdf import netcdf_file
from distutils.version import StrictVersion
#from scipy import version as scipy_version
from copy import deepcopy, copy
#import fdiff

_lbldis_req_key_list = ('v_start', 'v_end', 'v_delta',
                        'clddef', 'lblrtm_dir', 'solar_datafile', 
                        'ssp_db_files')
_lbldis_clddef_req_key_list = ('dbnum', 'r_eff', 'ref_wn', 'tau')

#_psd_function_list = {'gamma':psd_functions.gamma_psd, 
#                      'HT_gamma':psd_functions.HT_gamma_psd, 
#                      'mult_gamma':psd_functions.mult_gamma_psd, 
#                      'mult_HT_gamma':psd_functions.mult_HT_gamma_psd,
#                      'power_law':psd_functions.power_law_psd}

def get_example_kw(get_profile=False):
    """
    Return a sample set of keywords for use with the run() method.
    For illustrative purposes only.
    It generates a nadir pointed, upwelling radiance spectrum from 
    600 - 2000 wavenumbers at 0.5 resolution, using the standard 
    tropical atmosphere profile.

    A single cloud layer is included, at z=8 km, tau=2, using 
    mie ice properties; note the path to the ice property file 
    and solar spectrum files are hardcoded and will need to be 
    updated to reflect the actual data paths.

    Input keyword specifies whether profile data should be used 
    (it attempts to load the standard tropical atm profile from the 
    included testdata); if false, then the keyword set will include 
    just the minimum description, which includes a set of levels, 
    with the remaining information pulled from whatever data is 
    hardcoded into the LBLRTM for the standard atmosphere data.
    """

    kw = {}
    kw['v_start'] = 600.0
    kw['v_end'] = 2000.0
    kw['v_delta'] = 0.5
    kw['atm_flag'] = 1
    kw['nmol'] = 7

    kw['solar_datafile'] = (os.path.expanduser('~/bin/') + 
                            'unit_solar_spectrum.1cm.asc')
##  kw['ssp_db_files'] = ['/data/lbldis_ssp_db/' + 
##                          'ssp_db.mie_ice.gamma_sigma_0p100',]
    kw['ssp_db_files'] = ['/data/wsessions/lbldis_inputs/' + 
                          'ssp_db.shettle_dust.gamma_sigma_0p100',]
    
    if get_profile:
        # look for the testdata in this file's parent dir structure.
        # I'm not sure if this is "kosher" to do this via __file__, 
        # But I think this should work assuming setup.py was used.
        file_dir, file_name = os.path.split(__file__)
        main_root_dir, sub_dir = os.path.split(file_dir)
        testdata_dir = os.path.join(main_root_dir, 'testdata')
        dat = np.load(os.path.join(testdata_dir, 'std_atm_data.npz'))
        profile = dat['profile'][:,0].copy()
        dat.close()
        kw['profile'] = profile
        kw['levels'] = profile['alt']
    else:
        kw['levels'] = np.r_[np.linspace(0,15,31), np.linspace(16,32,17)]

    clddef = []
    clddef.append({'dbnum':0, 'r_eff':20.0, 'ref_wn':0.0, 
                   'z':8.0, 'tau':[2.0,]})
    kw['clddef'] = clddef

    return kw


def flux_spectral_integ(lbldis_data, wn_range=None):
    """
    Integrate the spectral flux profile from the LBLDIS run.
    Optionally, limit the calculation to a certain wavenumber range, 
    otherwise the calculation is done for the entire wavenumber coverage 
    of the run.

    input:
    lbldis_dat - the python dictionary returned from an LBLDIS run for 
        a spectral profile. This should contain spectral flux profiles with 
        dimension [nwave, nlevel, ncase].

    Returns a dictionary with four [nlevel, ncase] ndarrays:
    flux_up_profile_integ
    flux_down_profile_integ
    flux_net_profile_integ
    flux_beam_profile_integ
    heating_rate (in K/day)

    note that flux_net = flux_up - flux_down
    """

    in_flds = ('flux_up_profile','flux_down_profile','flux_beam_profile')
    if False in map(lbldis_data.has_key, in_flds):
        raise ValueError, 'Missing required flux data in lbldis_dat'

    nwave, nlevel, ncase = lbldis_data['flux_up_profile'].shape

    x = lbldis_data['wnum']
    if wn_range is not None:
        wn_slice = slice(x.searchsorted(wn_range[0]), 
                         x.searchsorted(wn_range[1])-1)
    else:
        wn_slice = slice(None)
    x = x[wn_slice]

    out_data = {}
    for f in in_flds:
        out_data[f+'_integ'] = np.zeros((nlevel, ncase))
        for k, c in np.ndindex((nlevel,ncase)):
            y = lbldis_data[f][wn_slice,k,c]
            out_data[f+'_integ'][k,c] = np.trapz(y,x)

    # convert mW -> W
    for f in in_flds:
        out_data[f+'_integ'] = out_data[f+'_integ'] * 1e-3

    # compute net flux
    out_data['flux_net_profile_integ'] = \
        out_data['flux_up_profile_integ'] - out_data['flux_down_profile_integ']

    # compute heating; use virtual temp. approx.
    # H = - (dF/dz) / (rho C_p)
    # H = - (dF/dz) * (R_d (1+0.6q) T / P) (1 / C_p)

    # Petty page 315. 
    C_p = 1005.0
    R_d = 287.047
    rho_inv = R_d * (1 + 0.6e-3 * lbldis_data['mixing_ratio']) * \
        lbldis_data['temperature'] / lbldis_data['pressure']

    # units: 1e-2 for hPa -> Pa; 86400. for K/s to K/d; 
    # 1e-3 for km to m in the dz
    A = rho_inv * (1e-2 * 86400 * 1e-3 / C_p)
    # central averages (layers)
    A = 0.5*(A[1:] + A[:-1])
    dz = np.diff(lbldis_data['height'])

    out_data['heating_rate'] = np.zeros((nlevel, ncase))
    for c in range(ncase):
        dFdz = np.diff(out_data['flux_net_profile_integ'][:,c]) / dz
        # prepend a zero so it is the same length as the other arrays, 
        # so we can be lazy while plotting.
        out_data['heating_rate'][:,c] = np.r_[0.0, - A * dFdz]

    return out_data

def fdrun(**kwargs):
    """
    Helper to run default finite diff derivatives.
    Uses new separate code in the lbldis_utils.fdiff module.

    This "standard" run takes the profile, cldprop, and bdryprop derivatives, 
    and then collects the results into one PyDict.
    It sets the deriv_list to: ['temp', 'H2O', 'cld_tau', 'cld_r_e', 
       'surf_temp', 'surf_emis']
    """
    K = {}
    deriv_list = ['temp','H2O','cld_tau','cld_r_e','surf_temp','surf_emis','z']

    out1 = fdiff.profile(deriv_list=deriv_list,**kwargs)
    for k in ['temp', 'H2O']:
        K[k] = out1['K'][k]

    out1 = fdiff.cldprop(deriv_list=deriv_list,**kwargs)
    for k in ['cld_r_e', 'cld_tau']:
        K[k] = out1['K'][k]

    out1 = fdiff.bndprop(deriv_list=deriv_list,**kwargs)
    for k in ['surf_temp', 'surf_emis', 'cld_z']:
        K[k] = out1['K'][k]

    return {'wn':out1['wn'], 'K':K}

def fdiff_psd_run(**inp_kwargs):
    """
    fdiff_psd_run(): finite difference derivative of LBLDIS forward model.
    This will take derivatives with respect to PSD parameters; the 
    implementation is sufficiently different from fdiff run that I decided 
    to make this a separate function. This is due to the very different 
    way that a PSD parameter is handled - a temporary ssp_db file needs to 
    be made, and the r_eff needs to be determined from that.

    Inputs are the same as run(), except that the ssp_db_files keyword is not 
    used. Note that derivatives will be taken for each psd parameter, 
    automatically.

    There are two additional input keywords:
    psd_descr: a PyDict containing a single psd description (see 
       integrate_ssp_db().)
    sc_db: a netcdf file of single-size scattering properties
        (specifically, a Baum & Yang database)

    Note that this does not work for multiple tau or clddef 
    at this point.
    
    """

    kwargs = deepcopy(inp_kwargs)
    for k in _lbldis_req_key_list:
        if not kwargs.has_key(k):
            # we don't care about the lblrtm_dir here, so don't raise 
            # exception for that case
            if (k != 'lblrtm_dir') & (k != ssp_db_files):
                raise ValueError, 'Missing required keyword ' + k

    n_cld = len(kwargs['clddef'])
    n_tau = len(kwargs['clddef'][0]['tau'])
    n_param = len(kwargs['clddef'][0]['psd_desc']['params'])

    if n_cld > 1:
        raise NotImpementedError, 'not implemented for more than ' + \
            'a single cloud layer'

    # not sure yet what makes sense for the increments - these are 
    # hard-coded, additive
    param_increment = 0.01
    norm = 1.0/param_increment

    base_run = psd_run(**kwargs)
    wn = base_run['wnum']
    K = np.zeros( (wn.shape[0], n_tau*n_param) )

    psd_param = kwargs['clddef'][0]['psd_desc']['params']
    for n,p in enumerate(psd_param):
        psd_param[p] += param_increment
        perturb_run = psd_run(**kwargs)        
        psd_param[p] -= param_increment
        K[:,n*n_param:(n+1)*n_param] = \
            norm*(perturb_run['radiance']-base_run['radiance'])

    return {'wn':wn, 'K':K}

def fdiff_ssp_run(**inp_kwargs):
    """
    fdiff_ssp_run(): finite difference derivative of LBLDIS forward model.
    This variant takes derivatives with respect to bulk ssp parameters
    (g or w), for a run using HG phase functions.
    
    Inputs are the same as run(), with one input parameter (which must be first):
    deriv_list: a list of strings denoting which derivatives should be taken.
        The strings are checked in the following order, so this order fixes 
        the output order of the jacobian array.
        must be a subset of these options: ['g', 'w']
    
    """

    # parse out the inputs to determine which derivatives we want; 
    # then, combine everything into what should be a single LBLDIS 
    # run. This should be the most efficient method, since then the 
    # LBLRTM load only occurs once.

    kwargs = deepcopy(inp_kwargs)
    for k in _lbldis_req_key_list:
        if not kwargs.has_key(k):
            # we don't care about the lblrtm_dir here, so don't raise 
            # exception for that case
            if k != 'lblrtm_dir':
                raise ValueError, 'Missing required keyword ' + k
    if not ('deriv_list' in kwargs.keys()):
        raise ValueError, 'Missing required keyword deriv_list'

    n_cld = len(kwargs['clddef'])
    n_ssp_db = len(kwargs['ssp_db_files'])
    n_tau = len(kwargs['clddef'][0]['tau'])

    # not sure yet what makes sense for the increments - these are 
    # hard-coded, additive
    g_increment = 0.01
    w_increment = 0.01
    norms = []

    # unfortunately, this repeat some of the top of .run(); need to set 
    # this up first, since we will create some "phantom" ssp_db files 
    # in the wrk directory.

    clean_temp_wrk_dir = kwargs.get('clean_temp_wrk_dir', True)
    wrk_dir = kwargs.get('wrk_dir', None)
    proc_wrap = kwargs.get('proc_wrap', None)

    # rest should lbldis params, or params for LBLDIS.
    if proc_wrap is None:
        proc_wrap = processWrapper(wrk_dir = wrk_dir, 
                                   clean_temp_wrk_dir = clean_temp_wrk_dir, 
                                   name = 'lbldis_run_tmp')
    kwargs['proc_wrap'] = proc_wrap

    tmp_dir = proc_wrap.getcwd()

    if 'g' in kwargs['deriv_list']:

        # make the incremented ssp data
        for n in range(n_ssp_db):

            base_ssp_file = kwargs['ssp_db_files'][n]
            incr_ssp_file = \
                os.path.join(tmp_dir, 
                             os.path.basename(base_ssp_file) + '.g_incr')

            s = ssp_db.ssp_db(base_ssp_file)
            s['ssp']['g'] = s['ssp']['g'] - g_increment
            ssp_db.write_ssp_db(incr_ssp_file, s, 'autogenerated by ' + 
                                'lbldis_utils.fdiff_ssp_run()')

            kwargs['ssp_db_files'].append(incr_ssp_file)

        for c in range(n_cld):

            # copy a clddef and append, this will contain the incremented 
            # ssp data for g' = g + dg
            kwargs['clddef'].append(deepcopy(kwargs['clddef'][c]))

            # append zero to the base state set tau list (as padding 
            # for the incremented set of cld properties)
            for n in range(n_tau):
                kwargs['clddef'][c]['tau'].append(0)

            # update ssp number, and prepend a padding zero for cld tau, 
            # for the new clddef.
            kwargs['clddef'][-1]['dbnum'] += n_ssp_db
            for n in range(n_tau):
                kwargs['clddef'][-1]['tau'].insert(0,0)

            # update norms
            for n in range(n_tau):
                norms.append(g_increment)


    if 'w' in kwargs['deriv_list']:

        raise NotImplementedError, "You broke this"

        # make the incremented ssp data
        norms.append(w_increment)
        for n in range(n_ssp_db):

            base_ssp_file = kwargs['ssp_db_files'][n]
            incr_ssp_file = \
                os.path.join(tmp_dir, 
                             os.path.basename(base_ssp_file) + '.w_incr')

            s = ssp_db.ssp_db(base_ssp_file)
            s['ssp']['w0'] = s['ssp']['w0'] - w_increment
            ssp_db.write_ssp_db(incr_ssp_file, s, 'autogenerated by ' + 
                                'lbldis_utils.fdiff_ssp_run()')

            kwargs['ssp_db_files'].append(incr_ssp_file)

        for c in range(n_cld):

            # copy a clddef and append, this will contain the incremented 
            # ssp data for w' = w + dw
            kwargs['clddef'].append(deepcopy(kwargs['clddef'][c]))

            # append zero to the base state set tau list (as padding 
            # for the incremented set of cld properties)
            kwargs['clddef'][c]['tau'].append(0)

            # if incremented cld defs were added for g, we need to append 
            # another zero to the cld tau list
            if 'g' in kwargs['deriv_list']:
                kwargs['clddef'][c+n_cld]['tau'].append(0)

                # update ssp number, and prepend a padding zero for cld tau, 
                # for the new clddef. 
                kwargs['clddef'][-1]['dbnum'] += 2*n_ssp_db
                kwargs['clddef'][-1]['tau'].insert(0,0)
                # this swaps the last 2 items
                kwargs['clddef'][-1]['tau'][-2:] = \
                    kwargs['clddef'][-1]['tau'][-1:-3:-1]
            else:
                kwargs['clddef'][-1]['dbnum'] += n_ssp_db
                kwargs['clddef'][-1]['tau'].insert(0,0)

    kwargs['use_HG'] = True
    dat = run(**kwargs)

    K = np.zeros( (dat['radiance'].shape[0], n_tau) )

    for c in range(n_tau):
        K[:,c] = (dat['radiance'][:,c+n_tau] - dat['radiance'][:,c]) / norms[c]

    return {'wn':dat['wnum'], 'K':K}    

def fdiff_run(**kwargs):
    """
    fdiff_run(): finite difference derivative of LBLDIS forward model.
    Inputs are the same as run(), with one additional keyword parameter:
    deriv_list: a list of strings denoting which derivatives should be taken.
        The strings are checked in the following order, so this order fixes 
        the output order of the jacobian array.
        must be a subset of these options:
        'cld_r_e' - r_eff parameters in the clddef
        'cld_tau' - tau parameters in the clddef
        'temp_profile' - the atmospheric temperature profile (per K)
        'h2o_profile' - the atmospheric water vapor profile (per ln change)
        'surf_temp' - surface temperature (per K)
        'surf_emis' - surface emissivity (per 1% change)

    Note that the stored_lblrtm_run keyword might be of interest here 
    if the finite difference is repeated called on the same LBLRTM 
    gas/temperature profile.

    """

    # parse out the inputs to determine which derivatives we want; 
    # then, combine everything into what should be a single LBLDIS 
    # run. This should be the most efficient method, since then the 
    # LBLRTM load only occurs once.

    lbldis_kwargs = deepcopy(kwargs)
    for k in _lbldis_req_key_list:
        if not lbldis_kwargs.has_key(k):
            # we don't care about the lblrtm_dir here, so don't raise 
            # exception for that case
            if k != 'lblrtm_dir':
                raise ValueError, 'Missing required keyword ' + k
    if not 'deriv_list' in kwargs.keys():
        raise ValueError, 'Missing required keyword deriv_list'


    # unfortunately, this repeats some of the top of .run(); need to set 
    # this up first, in order to use only one processWrapper for all 
    # the different runs

    clean_temp_wrk_dir = kwargs.get('clean_temp_wrk_dir', True)
    wrk_dir = kwargs.get('wrk_dir', None)
    proc_wrap = kwargs.get('proc_wrap', None)

    if proc_wrap is None:
        proc_wrap = processWrapper(wrk_dir = wrk_dir, 
                                   clean_temp_wrk_dir = clean_temp_wrk_dir, 
                                   name = 'lbldis_run_tmp', 
                                   debug = True)
    lbldis_kwargs['proc_wrap'] = proc_wrap

    ncld = len(lbldis_kwargs['clddef'])
    ntau = len(lbldis_kwargs['clddef'][0]['tau'])

    if (('cld_r_e' in lbldis_kwargs['deriv_list']) or 
        ('cld_tau' in lbldis_kwargs['deriv_list'])):

        # not sure yet what makes sense for the increments - these are 
        # hard-coded, multiplicative (e.g., 0.01 = 1% change)
        tau_increment = 0.01
        r_e_increment = 0.10
        base_clddef = lbldis_kwargs['clddef']

        if 'cld_r_e' in lbldis_kwargs['deriv_list']:
            r_e_cases, r_e_norm = \
                _add_r_e_cases(base_clddef, r_e_increment)
        else:
            r_e_cases, r_e_norm = [], []
            
        if 'cld_tau' in lbldis_kwargs['deriv_list']:
            tau_cases, tau_norm = \
                _add_tau_cases(base_clddef, tau_increment)
        else:
            tau_cases, tau_norm = [], []
        
        # the tau_cases are merged onto r_e_cases.
        lbldis_kwargs['clddef'] = _merge_cases(r_e_cases, tau_cases)
        norms = r_e_norm + tau_norm

        dat = run(**lbldis_kwargs)

        K = np.zeros( (dat['radiance'].shape[0], len(norms)) )

        if 'cld_r_e' in lbldis_kwargs['deriv_list']:
            for t in range(ntau):
                for c in range(ncld):
                    k_ind = t + c*ntau
                    r_ind = k_ind + ntau
                    K[:,k_ind] = (dat['radiance'][:,r_ind] -
                                  dat['radiance'][:,t]) / norms[k_ind]
        if 'cld_tau' in lbldis_kwargs['deriv_list']:
            for t in range(ntau):
                for c in range(ncld):
                    if 'cld_r_e' in kwargs['deriv_list']:
                        k_ind = t + ncld*ntau + c*ntau
                    else:
                        k_ind = t + c*ntau
                    r_ind = k_ind + ntau
                    K[:,k_ind] = (dat['radiance'][:,r_ind] -
                                  dat['radiance'][:,t]) / norms[k_ind]

        # after finishing, replace with the base_clddef (later derivs 
        # will need it)
        lbldis_kwargs['clddef'] = base_clddef

    else:
        # need the "base" run for the other derivs
        dat = run(**lbldis_kwargs)
        # make a place holder for the later concatenations for K
        K = np.zeros( (dat['radiance'].shape[0], 0) )

    # for profiles, modify data, using the base clddef
    profile = lbldis_kwargs['profile']
    nlvl = profile.shape[0]
    if 'temp_profile' in kwargs['deriv_list']:
        
        temp_incr = 0.1
        Ktemp = np.zeros( (K.shape[0], nlvl*ntau) )

        # otherwise, the jacobian at level 0 would also increase the
        # surface temperature (making a much larger signal there)
        lbldis_kwargs['surf_temp'] = profile['temp'][0]

        for n in range(nlvl):
            profile['temp'][n] = profile['temp'][n] + temp_incr
            dat_n = run(**lbldis_kwargs)
            for t in range(ntau):
                Ktemp[:,n + t*nlvl] = (dat_n['radiance'][:,t] -
                                       dat['radiance'][:,t])
            profile['temp'][n] = profile['temp'][n] - temp_incr

        Ktemp = Ktemp / temp_incr
        K = np.c_[ K, Ktemp ]

    if 'h2o_profile' in kwargs['deriv_list']:

        h2o_incr = 1.01
        profile = lbldis_kwargs['profile']
        Kh2o = np.zeros( (K.shape[0], nlvl*ntau) )

        for n in range(nlvl):
            profile['H2O'][n] = profile['H2O'][n] * h2o_incr
            dat_n = run(**lbldis_kwargs)
            for t in range(ntau):
                Kh2o[:,n + t*nlvl] = (dat_n['radiance'][:,t] -
                                      dat['radiance'][:,t])
            profile['H2O'][n] = profile['H2O'][n] / h2o_incr

        Kh2o = Kh2o / np.log(h2o_incr)
        K = np.c_[ K, Kh2o ]

    if 'surf_temp' in kwargs['deriv_list']:
        temp_incr = 0.1
        if 'surf_temp' in lbldis_kwargs:
            lbldis_kwargs['surf_temp'] = \
                lbldis_kwargs['surf_temp'] + temp_incr
        else:
            lbldis_kwargs['surf_temp'] = \
                lbldis_kwargs['profile']['temp'][0] + temp_incr

        dat_n = run(**lbldis_kwargs)
        lbldis_kwargs['surf_temp'] = lbldis_kwargs['surf_temp'] - temp_incr
        KTs = np.zeros((K.shape[0], ntau))
        for t in range(ntau):
            KTs[:,t] = dat_n['radiance'][:,t] - dat['radiance'][:,t]
        KTs = KTs / temp_incr
        K = np.c_[ K, KTs ]

    if 'surf_emis' in kwargs['deriv_list']:
        emis_incr = 0.001
        if 'surf_emis' in lbldis_kwargs:
            # might be a list? or ndarray
            # Note that the LBLRTM utils may not deal with the spectrally 
            # defined emissivity at this point...
            for n in range(len(lbldis_kwargs['surf_emis'])):
                lbldis_kwargs['surf_emis'][n,1] = \
                    lbldis_kwargs['surf_emis'][n,1] - emis_incr
        else:
            # beware - hardcoded, copy/paste from the default used in 
            # lbldis.run()
            surf_emis_arr = np.c_[[100.0,3000.0],[0.985,0.985]]
            surf_emis_arr[:,1] -= emis_incr
            lbldis_kwargs['surf_emis'] = surf_emis_arr

        dat_n = run(**lbldis_kwargs)
        for n in range(len(lbldis_kwargs['surf_emis'])):
            lbldis_kwargs['surf_emis'][n,1] = \
                lbldis_kwargs['surf_emis'][n,1] + emis_incr
        KTe = np.zeros((K.shape[0], ntau))
        for t in range(ntau):
            KTe[:,t] = dat_n['radiance'][:,t] - dat['radiance'][:,t]
        # want this relative to 1%; minus sign takes care of the fact that 
        # the emis_incr was subtracted (which was done to make sure emis 
        # does not increase over 1.0)
        KTe = KTe * (-0.01 / emis_incr)
        K = np.c_[ K, KTe ]

    return {'wn':dat['wnum'], 'K':K}

def _add_tau_cases(inp_clddefs, tau_incr):
    clddefs = deepcopy(inp_clddefs)
    ncld = len(clddefs)
    ntau = len(clddefs[0]['tau'])
    nadd = ncld*ntau
    tau_norm = []
    for c1 in range(ncld):
        for c2 in range(ncld):
            for t in range(ntau):
                newtau = clddefs[c1]['tau'][t]
                if c1==c2:
                    tau_norm.append(tau_incr * newtau)
                    newtau *= 1.0 + tau_incr
                clddefs[c1]['tau'].append(newtau)
    return clddefs, tau_norm

def _add_r_e_cases(inp_clddefs, r_e_incr):

    clddefs = deepcopy(inp_clddefs)
    ncld = len(clddefs)
    ntau = len(clddefs[0]['tau'])
    nadd = ncld*ntau
    r_e_norm = []
    dummy = [0.0 for n in range(ntau)]
    for c1 in range(ncld):
        for c2 in range(ncld):
            new_cld = deepcopy(clddefs[c2])
            for t in range(len(new_cld['tau'])):
                new_cld['tau'][t] = 0.0
            if c1==c2:
                r_e_norm.extend(
                    [r_e_incr*new_cld['r_eff'] for n in range(ntau)])
                new_cld['r_eff'] *= 1.0 + r_e_incr
            for c3 in range(ncld):
                if c1==c3:
                    new_cld['tau'].extend(clddefs[c2]['tau'])
                else:
                    new_cld['tau'].extend(dummy)
            clddefs.append(new_cld)

    for c1 in range(ncld):
        for c2 in range(ncld):
            clddefs[c1]['tau'].extend(dummy)

    return clddefs, r_e_norm

def _merge_cases(r_e_cases, tau_cases):
    # This merges the two lists of cases, and onto the r_e_cases variable
    # Three possibilities: 1) no r_e_cases, in which case the merged list is 
    # equal to tau_cases;
    # possibility 2) both case lists are non-empty - so do the 
    # complicated merging.
    # possibility 3): no tau_cases, so the r_e_cases already represent 
    # the merged list, so we don't have to do anything.
    if len(r_e_cases) == 0:
        merged_cases = tau_cases
    elif len(tau_cases) > 0:
        merged_cases = deepcopy(r_e_cases)
        # Basically, append the tau cases to the end of the first set of 
        # r_e cases (these have the unmodified r_e values). This required 
        # dummy zeros to pad out the added tau_cases for the modified r_e's.

        # Number of tau cases = ntau (original set) + ncld*ntau
        # only want to copy over the [ntau:] elements of each tau list 
        # in the clddef.
        ncld = len(tau_cases)
        ntau = len(tau_cases[0]['tau']) / (1 + ncld)
        for n, t in enumerate(tau_cases):
            merged_cases[n]['tau'].extend(t['tau'][ntau:])
        # padding zeros, with the same length as what was added to the lists 
        # just above.
        dummy = [0.0 for q in range(len(tau_cases[0]['tau'])-ntau)]
        for n in range(len(tau_cases), len(r_e_cases)):
            merged_cases[n]['tau'].extend(dummy)
    else:
        merged_cases = r_e_cases

    return merged_cases

def psd_run(**kwargs):
    """
    psd_run(): wrapper to run LBLDIS, using raw single-size 
       scattering properties and a size PSD as input, rather than 
       ssp_db integrated bulk ssp.

    inputs:
    same as run(), except for kw related to cloud definition and 
    ssp_db, as follows:

    instead of clddef containing:
    {'dbnum':0, 'z' :10.0, 'ref_wn':0.0, 'tau':[1.0,], 'r_eff':10.0}
    The r_eff field is replaced with 'psd_descr', containing 
    valid PSD information suitable for integrate_ssp_db.

    instead of a keyword for ssp_db_files, use a keyword nc_files 
    containing the single size scattering property netcdf file.

    These are translated into a clddef suitable for run().

    outputs:
    Dictionary with all fields written into the LBLDIS netcdf file.

    example:
    lbldis_data = lbldis_utils.run(output_file='example.cdf', 
                                   **lbldis_kwargs)

    """

    # pull out clean_temp_wrk_dir, wrk_dir, output_file, proc_wrap keywords, 
    # and delete if present (since we are sending 
    # these directly to the LBLRTM TAPE5 creator)
    clean_temp_wrk_dir = kwargs.get('clean_temp_wrk_dir', True)
    if kwargs.has_key('clean_temp_wrk_dir'):
        del kwargs['clean_temp_wrk_dir']

    wrk_dir = kwargs.get('wrk_dir', None)
    if kwargs.has_key('wrk_dir'):
        del kwargs['wrk_dir']

    output_file = kwargs.get('output_file', None)
    if kwargs.has_key('output_file'):
        del kwargs['output_file']

    proc_wrap = kwargs.get('proc_wrap', None)
    if kwargs.has_key('proc_wrap'):
        del kwargs['proc_wrap']

    # rest should lbldis params, or params for LBLDIS.
    if proc_wrap is None:
        proc_wrap = processWrapper(wrk_dir = wrk_dir, 
                                   clean_temp_wrk_dir = clean_temp_wrk_dir, 
                                   name = 'lbldis_run_tmp')

    tmp_dir = proc_wrap.getcwd()

    kwargs['filename'] = os.path.join(tmp_dir, 'lbldis_input')

    # now, translate psd_descr, etc, into ssp_db, and update the 
    # kwargs accordingly.
    new_clddef = []
    ssp_db_files = []
    for n,c in enumerate(kwargs['clddef']):
        nc_file = kwargs['nc_files'][c['dbnum']]
        tmp_ssp_db_file = os.path.join(tmp_dir, 'ssp_tmp_{0:02d}'.format(n))
        psd_desc = [deepcopy(c['psd_desc']), 
                    deepcopy(c['psd_desc']), 
                    deepcopy(c['psd_desc'])]
        # kludge here - this will break for anything other than 
        # gamma PDF / HT PDF (Im not sure of another better way to do this.)
        if psd_desc[0]['params'].has_key('De'):
            psd_desc[0]['params']['De'] = psd_desc[0]['params']['De'] - 0.01
            psd_desc[2]['params']['De'] = psd_desc[2]['params']['De'] + 0.01
        elif psd_desc[1]['params'].has_key('mu'):
            psd_desc[0]['params']['mu'] = psd_desc[0]['params']['mu'] - 0.001
            psd_desc[2]['params']['mu'] = psd_desc[2]['params']['mu'] + 0.001
        else:
            raise NotImplementedError, 'Unfortunately the different ' + \
                'psd method increments are hardcoded, and neither known ' + \
                'psds were used...'
        # issue here - most work is expected to be done inside the proc_wrap
        # dir, except that this integrate function is entirely within python,
        # so will be using the python current dir (which is not the proc_wrap
        # dir). So, change / change back here, so integrate runs in the same
        # place as proc_wrap executed programs (e.g. the LBLDIS); put this in
        # a finally so that there is no cd-side effect if the integrate crashes.
        saved_dir = os.getcwd()
        os.chdir(tmp_dir)
        try:
            data, bulk_data = \
                  integrate_ssp_db(nc_file, psd_desc, tmp_ssp_db_file, 
                                   wn_range = [kwargs['v_start'],
                                               kwargs['v_end']])
        finally:
            os.chdir(saved_dir)
        new_clddef.append({})
        for k in ['ref_wn', 'tau', 'z', 'zbnd']:
            # this is to account for the fact that z or zbnd may be defined
            # (perhaps not both)
            if k in c.keys():
                new_clddef[n][k] = c[k]
        new_clddef[n]['dbnum'] = n
        new_clddef[n]['r_eff'] = bulk_data['r_eff'][1]
        ssp_db_files.append(tmp_ssp_db_file)
        
    kwargs['clddef'] = new_clddef
    kwargs['ssp_db_files'] = ssp_db_files

    # call helper to adjust levels (this winds up being somewhat complicated, 
    # splitting into a function gives easier testing)
    clddef, levels = _add_cld_levels(kwargs['clddef'], kwargs['levels'])
    kwargs['levels'] = levels
    kwargs['clddef'] = clddef

    # call LBLRTM - will run inside the above processWrapper's tmp dir, 
    # in the 'lblrtm_output' subdirectory.
    # Note we need to run this just after the _add_cld_levels() runs, 
    # since that might add a level.
    if kwargs.has_key('stored_lblrtm_run'):
        # do one sanity check - does this dir exist, and contain a TAPE7?
        # otherwise you are on your own...
        stored_lblrtm_run = kwargs['stored_lblrtm_run']
        if not os.access(os.path.join(stored_lblrtm_run, 'TAPE7'), os.R_OK):
            raise ValueError, \
                'stored_lblrtm_run, ' + stored_lblrtm_run + \
                ' did not contain a readable TAPE7'
        kwargs['lblrtm_dir'] = stored_lblrtm_run
    else:
        kwargs['lblrtm_dir'] = os.path.join('.', 'lblrtm_output')
        _run_od_lblrtm(tmp_dir, proc_wrap = proc_wrap, **kwargs)

    create_lbldis_inputfile(**kwargs)

    lbldis_exec =  kwargs.get('lbldis_exec', None)
    if lbldis_exec is not None:
        r = proc_wrap.runcmd(lbldis_exec + ' lbldis_input 0 lbldis_output')
    else:
        r = proc_wrap.runcmd('lbldis lbldis_input 0 lbldis_output')

    output_data = _load_lbldis_nc_file(tmp_dir + os.sep + 'lbldis_output.cdf')

    if output_file is not None:
        os.rename(os.path.join(tmp_dir, 'lbldis_output.cdf'), output_file)

    return output_data

def run(**kwargs):
    """
    run(): wrapper to run LBLDIS and return the computed radiation data.

    inputs:    
    save_stdout: set to True to collect the stdout/err from LBLDIS.
    lbldis_exec: string path/name to lbldis executable. Optional, if
       left unspecified, 'lbldis' is assumed (e.g., it is defined on
       the path somewhere)
    clean_temp_wrk_dir: boolean specifying whether to clean up the temporary 
       working dir (default True)
    wrk_dir: string containing directory path if a certain work directory 
       is desired. default None causes a directory in the system tmp area 
       to be used.
    output_file: string path+filename desired for the output netcdf from 
       LBLDIS. This will copy the file out from the tmp wrk dir, into 
       the desired location. Note this is necessary if you want the actual 
       file but have clean_temp_wrk_dir set to True, since that stuff will be 
       deleted.
    proc_wrap: specify an already existing processWrapper object to use for 
       the LBLDIS run. Note, if this is input, then clean_temp_wrk_dir will 
       have no meaning, since the the wrk_dir will essentially be managed 
       by the given processWrapper object.
    stored_lblrtm_run: set to path of a precomputed LBLRTM optical depth 
        run that should be used by the LBLDIS.
    microwindow_list: a [N,2] shaped array, or N element iterable of pairs 
        of numbers, containing microwindow definitions (pairs of low, high, 
        in wavenumber). If this is specified, then v_start, v_end, v_delta, 
        need not be entered (these will be populated with the appropriate 
        dummy values.) This also implies an additional file, with the same 
        name as the LBLDIS input parameter file, but with "_mw" prepended, 
        will be created in the working directory.

    kw dict suitable for create_lbldis_inputfile
       (this should be input with **kwargs syntax, see example)
    
    Note that clddef inputs here can also have 'zbnd' instead of z 
    to define the cloud altitudes; this must be handled at the level 
    above create_lbldis_inputfile, since the zbnds change the level data sent 
    to LBLRTM as well.
    zbnd: top/bottom of cloud [km]; In this case, the leveling will be 
        modified to account for the cloud boundaries - so LBLRTM will run 
        a different set of levels than were specified with the input 
        keyword. If the zbnd brackets multiple levels, the cloud tau will 
        be linearly partitioned among the potential subcloud layers 
        (if the top/bottom enclose an already existing layer).
        Note, care must be taken in inputting multiple clouds - unexpected 
        results may occur if one cloud is specified with a single z versus 
        boundaries.
        either z or zbnd should be used (z is ignored if zbnd is used).

    outputs:
    Dictionary with all fields written into the LBLDIS netcdf file.

    if save_stdout is set to True, then the stdout from LBLRTM 
    is returned in the 'stdout' key. Note this is really a merge 
    of whatever is returned from the python subprocess module's stdout 
    and stderr.

    example:
    lbldis_data = lbldis_utils.run(output_file='example.cdf', 
                                   **lbldis_kwargs)

    """

    # pull out clean_temp_wrk_dir, wrk_dir, output_file, proc_wrap keywords, 
    # and delete if present (since we are sending 
    # these directly to the LBLRTM TAPE5 creator)
    clean_temp_wrk_dir = kwargs.get('clean_temp_wrk_dir', True)
    if kwargs.has_key('clean_temp_wrk_dir'):
        del kwargs['clean_temp_wrk_dir']

    wrk_dir = kwargs.get('wrk_dir', None)
    if kwargs.has_key('wrk_dir'):
        del kwargs['wrk_dir']

    output_file = kwargs.get('output_file', None)
    if kwargs.has_key('output_file'):
        del kwargs['output_file']

    proc_wrap = kwargs.get('proc_wrap', None)
    if kwargs.has_key('proc_wrap'):
        del kwargs['proc_wrap']

    save_stdout = kwargs.get('save_stdout', False)
    if kwargs.has_key('save_stdout'):
        del kwargs['save_stdout']

    # rest should lbldis params, or params for LBLDIS.
    if proc_wrap is None:
        proc_wrap = processWrapper(wrk_dir = wrk_dir, 
                                   clean_temp_wrk_dir = clean_temp_wrk_dir, 
                                   name = 'lbldis_run_tmp')

    tmp_dir = proc_wrap.getcwd()

    # Process microwindow data and make the temp mw file.
    if kwargs.has_key('microwindow_list'):
        microw_list = kwargs['microwindow_list']
        del kwargs['microwindow_list']
        microw_list = np.asarray(microw_list)
        # Valid dimensions: [N,2], or [2,N] (transpose), or [2,]
        if microw_list.ndim == 1:
            microw_list = microw_list[:, np.newaxis]
        if microw_list.ndim > 2:
            raise ValueError, 'microwindow_list cannot have more ' + \
                'than 2 dimensions'
        if microw_list.shape[1] != 2:
            if microw_list.shape[0] != 2:
                raise ValueError, 'microwindow_list must be [N,2] or [2,N]'
            else:
                microw_list = microw_list.T
        n_microw = microw_list.shape[0]
        # Note this does not change the values in the caller's kw dictionary
        kwargs['v_start'] = -1.0
        kwargs['v_end'] = 0.0
        kwargs['v_delta'] = 0.0
        kwargs['microwindow_file'] = 'lbldis_input_mw'
        # Note adding the microw_list to the kwargs is so that the 
        # call to _run_od_lblrtm() knows about the needed wn range
        kwargs['microwindow_list'] = microw_list
        _create_mw_list_file(kwargs['microwindow_list'], 
                             os.path.join(tmp_dir, kwargs['microwindow_file']))

    # Messy part:
    # redo layering if cloud thicknesses are given. This must be done before 
    # LBLRTM run, since it effects the leveling/layering done within LBLRTM.
    # First check to see if any cloud boundary definitions are present at all, 
    # to trap for error of the levels input being missing.
    cld_check = [cld.has_key('zbnd') for cld in kwargs['clddef']]
    if cld_check.count(True) > 0:
        if not kwargs.has_key('levels'):
            raise ValueError, 'LBLRTM levels keyword is required when ' + \
                'cloud boundaries are specified'

    # call helper to adjust levels (this winds up being somewhat complicated, 
    # splitting into a function gives easier testing)
    clddef, levels = _add_cld_levels(kwargs['clddef'], kwargs['levels'])
    kwargs['levels'] = levels
    kwargs['clddef'] = clddef

    # call LBLRTM - will run inside the above processWrapper's tmp dir, 
    # in the 'lblrtm_output' subdirectory.
    if kwargs.has_key('stored_lblrtm_run'):
        # do one sanity check - does this dir exist, and contain a TAPE7?
        # otherwise you are on your own...
        stored_lblrtm_run = kwargs['stored_lblrtm_run']
        if not os.access(os.path.join(stored_lblrtm_run, 'TAPE7'), os.R_OK):
            raise ValueError, \
                'stored_lblrtm_run, ' + stored_lblrtm_run + \
                ' did not contain a readable TAPE7'
        kwargs['lblrtm_dir'] = stored_lblrtm_run
    else:
        kwargs['lblrtm_dir'] = os.path.join('.', 'lblrtm_output')
        _run_od_lblrtm(tmp_dir, proc_wrap = proc_wrap, **kwargs)

    kwargs['filename'] = os.path.join(tmp_dir, 'lbldis_input')
    create_lbldis_inputfile(**kwargs)

    lbldis_exec =  kwargs.get('lbldis_exec', None)
    if lbldis_exec is not None:
        stdout, stderr = proc_wrap.runcmd(
            lbldis_exec + ' lbldis_input 0 lbldis_output')
    else:
        stdout, stderr = proc_wrap.runcmd(
            'lbldis lbldis_input 0 lbldis_output')

    stdout = stdout.split(os.linesep)
    stderr = stderr.split(os.linesep)

    # Define "normal" as a run that contains "DONE." at the end, inside 
    # the stdout (actually 3rd from end)
    if (stdout[-3].strip() != 'DONE.') or (stderr != ['']):
        print 'LBLDIS returned abnormally.'
        print 'Last 5 Console output lines from LBLDIS:'
        print ''
        for s in stdout[-5:]:
            print s.strip()
        print ''
        print 'Stderr:'
        print ''
        for s in stderr:
            print s.strip()
        print ''
        raise ValueError, 'LBLDIS Failed'

    # This is sort of a hack - older versions of the processWrapper didn't 
    # correctly hand back stdout/err so it was harded to determine if LBLDIS 
    # failed, so we just had to try to read the netCDF and see what happened.
    # I think that most cases should now be caught by the above test (since 
    # stdout should now be correctly returned), but leave this check in case 
    # there is still an unhandled exception.
    try:
        output_data = _load_lbldis_nc_file(
            tmp_dir + os.sep + 'lbldis_output.cdf')
    except IOError, theError:
        print 'LBLDIS netCDF file not created'
        print 'Last 5 Console output lines from LBLDIS:'
        print ''
        for s in stdout[-5:]:
            print s.strip()
        print ''
        raise IOError, theError

    if output_file is not None:
        os.rename(os.path.join(tmp_dir, 'lbldis_output.cdf'), output_file)

    if save_stdout:
        output_data['stdout'] = stdout
        output_data['stderr'] = stderr

    return output_data


def _create_mw_list_file(microw_list, microw_file):

    mfile = open(microw_file, 'w')
    mfile.write(str(microw_list.shape[0]) + os.linesep)
    for microw in microw_list:
        mfile.write('{0:f} {1:f}'.format(microw[0], microw[1]) + os.linesep)
    mfile.close()

    

def _load_lbldis_nc_file(ncfile):

    # simple wrapper around the ncfile method to load NetCDF variables.
    # note that copy is used, so there should be no issues with memmap()
    # handles left over.

    # Zara is running older scipy (0.7.1), and the netcdf wrapper has a
    # different syntax. So, here, if we have the older scipy (I don't actually
    # know what version holds the change), then automatically call the version
    # for the older netcdf wrapper.

    # again, this is just a guess - 0.7.1 and 0.9.x are the only versions I
    # have checked.
    if StrictVersion(scipy_version.version) < StrictVersion('0.8.0'):
        return _load_lbldis_nc_file_oldver(ncfile)

    # proceed with newer syntax.
    f = netcdf_file(ncfile, 'r')
    dat = {}
    for var_name in f.variables:
        # use copy to make sure there is no leftover memmap reference.
        dat[var_name] = f.variables[var_name].data.copy()
    f.close()
    return dat


def _load_lbldis_nc_file_oldver(ncfile):
    f = netcdf_file(ncfile, 'r')
    dat = {}
    for name, var in f.variables.items():
        if len(var.shape) == 0:
            dat[name] = var.getValue()
        else:
            # use copy to make sure there is no leftover memmap object
            # that will keep the netcdf file open. This is not the best
            # solution but I have not figured out how to correctly deal
            # with the memmaps in all SciPy versions.
            tmp_dat = var[:]
            dat[name] = tmp_dat.copy()
    f.close()
    return dat


def create_lbldis_inputfile(**kwargs):
    """
    create input text file for LBLDIS. inputs (via keywords), as follows:

    required input keywords:

    atm_flag: standard atmosphere flag for LBLRTM run
    v_start:  starting wavenumber
    v_end:    ending wavenumber
    v_delta:  wavenumber increment for spectral radiance calculation
    Note on v_start, v_end, v_delta: these can contain other values which 
        indicate microwindows, microwindow files, etc. See LBLDIS 
        instructions for more details. Note, if the run() method is 
        used, the microwindow input will he handled correctly.
    microwindow_file: Set to string containing filename of microwindow 
        definitions. This is a simple text file with the number of 
        microwindows (an integer) on the first line, and then low/high 
        wavenumber pairs, one per line, defining the microwindows.
    ssp_db_files: list of string path+filenames for ssp_db files
    solar_datafile: path+filename to solar spectrum file
    clddef: a list or tuple of dictionaries, containing the following 
        fields for each simulated cloud layer:

        dbnum: integer index into ssp_db_files list, specifying the ssp 
            database that should be used for scattering properties
        r_eff: effective radius of integrated particle mixture [um]
        tau: optical depth of cloud layer
        z: height of cloud [km]; In this case, the cloud will "fill" whatever 
            layer encloses this particular altitude.
        ref_wn: reference wavenumber for tau calculation; pick a non-positive 
            number to assume the short wavelength, geometric optics limit 
            (Q_ext = 2)

    optional input keywords (defaults noted where applicable):

    use_HG: boolean (default false), on whether to use HG phase functions 
        for all scattering calculations.
    filename: defaults to 'lbldis_input' in the current directory.
    solar_zen: solar zenith angle; value < 0 implies sun below horizon [deg]
    solar_azi: solar relative azimuth angle [deg]
    solar_dst: solar distance (in AU)
    obs_zen: observer zenith angle [deg]; 0 implies upwelling, 
        180 implies downwelling.
    nstreams: number of DISORT streams; default is 16
    lblrtm_dir: path to LBLRTM optical depth run. 
        default is './lblrtm_output'
    surf_temp: surface skin temperature; defaults to -1, which implies the 
        temperature at the bottom level of the LBLRTM run is reused.
    surf_emis: [n,2] shaped array, containing wavenumber (first column), 
        and emissivity (second column); this emissivity spectrum is 
        linearly interpolated to the simulated wavenumber range.
        default is a gray emis=0.985 surface; the array is defined as:
        [ [100, 0.985], [3000, 0.985] ]
    enable_flux_profile (default false): set to true to compute the per-level
        spectral flux profile.

    If any clouds were input with thicknesses (e.g., ztop/zbottom), then the 
    levels keyword is required; this keyword is normally just sent to LBLRTM

    """

    # defaults are implemented by hardcoding the defaults in calls 
    # to dict.get(); any params that have no sensible defaults will then 
    # cause ValueError exceptions if they are missing.

    ilines = []
    # to shorten syntax.
    l = ilines.append

    for k in _lbldis_req_key_list:
        if not kwargs.has_key(k):
            raise ValueError, \
                'input keyword list missing required key ' + k

    clddef = kwargs['clddef']

    for k in _lbldis_clddef_req_key_list:
        if not clddef[0].has_key(k):
            raise ValueError, \
                'input clddef dictionary missing required key ' + k
    # kludge to account for cloud altitude set via either 'z' or 'zbnd'
    if (not clddef[0].has_key('z')) and (not clddef[0].has_key('zbnd')):
        raise ValueError, 'input clddef dictionary must specify z or zbnd'

    use_HG = kwargs.get('use_HG', False)

    numtau = len(clddef[0]['tau'])
    numcld = len(clddef)
    numssp = len(kwargs['ssp_db_files'])    

    l('LBLDIS input file auto-generated ' + 
      'by lbldis_utils.create_lbldis_inputfile()')
    
    l( '{0:<16d}Number of streams'.format(kwargs.get('nstreams',16)) )
    l( '{0:4.1f} {1:4.1f} {2:4.1f} '.format(kwargs.get('solar_zen', -90.0), 
                                            kwargs.get('solar_azi', 0.0), 
                                            kwargs.get('solar_dst', 1.0)) + \
           'Solar zenith angle (deg), relative azimuth (deg), ' + \
           'solar distance (a.u.)' )

    l( '{0:<16.1f}'.format(kwargs.get('obs_zen', 0.0)) + \
           'Zenith angle (degrees): 0 -> upwelling, 180 -> downwelling' )
    if kwargs.has_key('microwindow_file'):
        l( '{0:7.2f} {1:7.2f} {2:6.3f} '.format(kwargs['v_start'], 
                                                kwargs['v_end'], 
                                                kwargs['v_delta']) + \
               kwargs['microwindow_file'])
    else:
        l( '{0:7.2f} {1:7.2f} {2:6.3f} '.format(kwargs['v_start'], 
                                                kwargs['v_end'], 
                                                kwargs['v_delta']) + \
               'v_start, v_end, and v_delta [cm-1]' )

    if kwargs.get('scattering_off', False):
        print 'Turning scattering off!'
        cloud_param_spec = -numtau
    else:
        cloud_param_spec = numtau

    l( '{0:<16d}'.format(cloud_param_spec) + 
       'Cloud parameter option flag: 0: ' + 
       'reff and numdens, >=1:  reff and tau' )
    l( '{0:<16d}Number of cloud layers'.format(numcld) )

    for n in range(numcld):

        dbnum = clddef[n]['dbnum']
        if use_HG:
            dbnum = dbnum + 1
            dbnum = -1 * dbnum
    
        p1 = '{0:1d} {1:8.4f} {2:8.5f} {3:5.0f}'.format(
            dbnum, clddef[n]['z'], clddef[n]['r_eff'], clddef[n]['ref_wn'])
        p2 = " ".join(['{0:10.5f}'.format(float(t))
                       for t in clddef[n]['tau']])
        l( p1 + p2 )

    l( kwargs['lblrtm_dir'] )
    l( kwargs['solar_datafile'] )

    l( '{0:<16d}Number of scattering property databases'.format(numssp) )
    for ssp_db_file in kwargs['ssp_db_files']:
        l(ssp_db_file)

    l( '{0:<12.2f}'.format(kwargs.get('surf_temp', -1.0)) + 
       'Surface temperature (specifying a negative value takes ' +
       'the value from profile)' )

    emis = kwargs.get('surf_emis', np.c_[[100.0,3000.0],[0.985,0.985]])
    numemis = emis.shape[0]

    l( '{0:<12d}Number of surface spectral '.format(numemis) + 
       'emissivity lines (wnum, emis)' )

    for n in range(numemis):
        l( '{0:8.3f} {1:8.6f}'.format(emis[n,0],emis[n,1]) )

    l('{0:d}'.format(kwargs.get('enable_flux_profile', False)))

    filename = kwargs.get('filename',  'lbldis_input')

    f = open(filename, 'w')
    for line in ilines:
        f.write(line + '\n')
    f.close()

def _add_cld_levels(clddef, levels):
    # helper for create_inputfile - changes leveling and clddef z values 
    # to account for cloud boundary definitions. Inputs should be the list 
    # of clddef, and the levels array (expected to be np.ndarray).

    # First, make a pass and collect all levels that need to be added; 
    # then recreate the levels array in one step.
    needed_levels = []
    for cld in clddef:
        if cld.has_key('zbnd'):
            z1, z2 = cld['zbnd']
            if np.all(levels != z1):
                if needed_levels.count(z1) == 0:
                    needed_levels.append(z1)
            if np.all(levels != z2):
                if needed_levels.count(z2) == 0:
                    needed_levels.append(z2)
    if len(needed_levels) > 0:
        # use PyList to make this syntax easy
        old_levels = levels.tolist()
        new_levels = old_levels + needed_levels
        new_levels.sort()
        new_levels = np.array(new_levels)
    else:
        # if there are no needed levels, we can just exit 
        return deepcopy(clddef), copy(levels)

    # Now, run through clddef list again, and either we have: 
    # a cloud that only fills one layer; in which case, set z to themidpoint;
    # a cloud that fills more than one layer; in that case, divide up the 
    # tau amount proportional to the Dz in each layer.
    new_clddef = []
    clddef_copy = deepcopy(clddef)
    for cld in clddef_copy:
        if cld.has_key('zbnd'):
            z1, z2 = cld['zbnd']
            i1 = new_levels.searchsorted(z1)
            i2 = new_levels.searchsorted(z2)
            # take abs, just in case z1,z2 are not in the order we expect 
            # (note the above code doesn't care which order)
            # first case - z bound does not contain any interior levels, 
            # so just set z to the midpoint.
            if abs(i1-i2) == 1:
                cld['z'] = 0.5*(z1 + z2)
            else:
                # otherwise, z bound contains >= 1 interior levels; so, 
                # linearly divide tau among the interior layers.
                # first, make sure we get top/bottom straight here.
                ibot = min(i1,i2)
                itop = max(i1,i2)
                nsubcld = itop-ibot
                zfractions = new_levels[ibot:itop+1].copy()
                zfractions -= zfractions[0]
                zfractions /= zfractions[-1]
                total_tau = copy(cld['tau'])
                # set original cloud to the bottom one; scale taus by the 
                # fraction of z bound within this layer.
                cld['z'] = 0.5*(new_levels[ibot]+new_levels[ibot+1])
                frac = zfractions[1]
                cld['tau'] = [t*frac for t in total_tau]
                # then loop through additional clouds, which are appended to 
                # new_clddef.
                for n in range(1,nsubcld):
                    new_cld = copy(cld)
                    new_cld['z'] = 0.5*(new_levels[ibot+n]+new_levels[ibot+n+1])
                    frac = zfractions[n+1]-zfractions[n]
                    new_cld['tau'] = [t*frac for t in total_tau]
                    new_clddef.append(new_cld)

    new_clddef = clddef_copy + new_clddef
    return new_clddef, new_levels

def write_ncdf(file, data):
    """
    write a netcdf file output from data dictionary.
    This is not a generalized routine, and should only be used with 
    the data dictionary returned from lbldis_utils.run().
    It is mainly intended to be a convenience function in case the user 
    would like the netcdf file after the run() method occurs; in addition, 
    this was a useful short test of the SciPy netcdf functions.

    Note if the LBLDIS file is really desired, see the output_file keyword 
    to lbldis_utils.run(). The netcdf file created here is stripped from 
    all of the nice metadata in the original file.
    """
    f = netcdf_file(file, 'w')
    f.createDimension('n_wnum', data['radiance'].shape[0])
    f.createDimension('n_cloud_layers', data['cld_tau'].shape[1])
    f.createDimension('n_instances', data['cld_tau'].shape[0])
    f.createDimension('n_height', data['height'].shape[0])
    f.createDimension('dummy', None)

    # hard coded mapping of the dimensions to key names - not ideal, but 
    # I'm not going to spend time to figure out a better method.
    var_dims = {'wnum': ('n_wnum',), 
                'radiance': ('n_wnum', 'n_instances'), 
                'flux_up': ('n_wnum', 'n_instances'), 
                'flux_down': ('n_wnum', 'n_instances'), 
                'flux_beam': ('n_wnum', 'n_instances'), 
                'cld_height': ('n_cloud_layers',), 
                'cld_tau': ('n_instances', 'n_cloud_layers'), 
                'cld_tau_wnum': ('n_cloud_layers',), 
                'cld_reff': ('n_cloud_layers',), 
                'cld_phase': ('n_cloud_layers',), 
                'sfc_temperature' : ('dummy',), 
                'sfc_emissivity': ('n_wnum',), 
                'pwv': ('dummy',), 
                'height': ('n_height',), 
                'pressure': ('n_height',), 
                'temperature': ('n_height',), 
                'mixing_ratio': ('n_height',)}

    for k in data.keys():
        if k == 'cld_phase':
            v = f.createVariable(k, 'short', var_dims[k])
        else:
            v = f.createVariable(k, 'single', var_dims[k])
        if var_dims[k][0] == 'dummy':
            v[0] = data[k]
        else:
            v[:] = data[k]

    f.close()

def _run_od_lblrtm(tmp_dir, **kwargs):
    # Helper to run LBLRTM to produce Optical Depth files.
    # The input can just be a copy of the kwargs into lbldis.run(); this 
    # function copies out the desired keys.
    # first parameter input must contain the desired directory parent 
    # where the LBLRTM run should exist - the subdir 'lblrtm_output' will 
    # be created under this name.
    lblrtm_subdir_name = 'lblrtm_output'
    lblrtm_subdir = os.path.join(tmp_dir, lblrtm_subdir_name)
    lblrtm_kwargs = {}
    lblrtm_kwargs['run_type'] = 'optdepth'
    if kwargs['v_start'] == -1:
        # if this run uses mw, then pick the wavenumber range as the range that 
        # closes all mw, plus the +/- 25 1/cm buffer
        mw_min = kwargs['microwindow_list'][:,0].min()
        mw_max = kwargs['microwindow_list'][:,1].max()
        lblrtm_kwargs['wmin'] = mw_min - 25.0
        lblrtm_kwargs['wmax'] = mw_max + 25.0
    else:
        # otherwise just copy the start/end with the 25 1/cm buffer
        lblrtm_kwargs['wmin'] = kwargs['v_start'] - 25.0
        lblrtm_kwargs['wmax'] = kwargs['v_end'] + 25.0
    lblrtm_kwargs['atm_flag'] = kwargs['atm_flag']

    for k in ['profile', 'levels', 'level_units', 'proc_wrap', 'nmol', 
              'continua_mult']:
        if kwargs.has_key(k):
            lblrtm_kwargs[k] = kwargs[k]

    # Make sure rayleigh is not included
    if lblrtm_kwargs.has_key('continua_mult'):
        if lblrtm_kwargs['continua_mult'][6] != 0.0:
            raise ValueError, 'continua_mult[6] represents rayleigh ' + \
                  'continuum which must be equal to zero since LBLDIS ' + \
                  'computes it independently'
    else:
        # set all continua to the default values (multiply by 1)
        # except for the rayleigh, which is disabled.
        lblrtm_kwargs['continua_mult'] = [1.0,1.0,1.0,1.0,1.0,1.0,0.0]
    if not os.access(lblrtm_subdir, os.R_OK):
        os.mkdir(lblrtm_subdir)

    cwd = kwargs['proc_wrap'].getcwd()
    kwargs['proc_wrap'].setcwd(lblrtm_subdir)
    lblrtm_utils.run(**lblrtm_kwargs)
    kwargs['proc_wrap'].setcwd(cwd)


def integrate_ssp_db(ncdb_file, psd_descriptions, ssp_db_file, 
                     wavenumbers=None, wn_range=None, 
                     save_npz=None):
    """
    integrate single-sized single scattering properties over a 
    particle size dist.

    saves a LBLDIS format ssp_db file.

    input single scattering properties are assumed to be in a particularly 
    defined netcdf file (variable names, dimensions, etc)

    Inputs:
    ncdb_file (Baum/Yang format netCDF database with single-size scattering 
        properties.
    psd_descriptions: list or tuple of PyDict, containing two fields:
        name: string name of PSD function ("HT_gamma", "gamma", etc. - 
            see psd_functions.py)
        params: PyDict containing the params for the associated PSD 
            function (see psd_functions.py)
    ssp_db_file: string path/filename for the created LBLDIS format ssp_db 
        file.

    optional input:
    wavenumbers: array of wavenumbers to extract from the single size 
        scattering database. Note that the closest wavenumber in the 
        single size database is selected, and then output in the ssp_db 
        file (meaning, the output wavnumber array won't match the requested 
        list of wavenumbers, since interpolation is not attempted)
        If more than one of the requested wavenumbers map to the same 
        wavenumber in this way, the duplicate will be dropped; in this 
        way, the output wavenumber list may be shorter than the requested
        wavenumber list.
        If not specified, a default list of wavenumbers is used.
    wn_range: a 2-element wavenumber range (low, high), that will trim 
        the wavenumber list. This is most useful when using the default 
        wavenumber list (e.g., not specifying the wavenumber keyword).
    save_npz: set to a string filename if you wish to write the full 
        data into a NumPy npz format file (in addition to the ssp).
        This contains additional information that does not fit into 
        the ssp_db format (since LBLDIS requires a particular text format 
        in that instance).

    """

    for psd in psd_descriptions:
        if psd['name'] not in _psd_function_list.keys():
            raise ValueError, 'psd_description ' + \
                str(psd['name']) + ' is not a known psd'

    # assign default wavenumber array, if needed.
    if wavenumbers is None:
        wavenumbers = np.r_[np.arange(100,600,10), np.arange(600,1700,20), 
                            np.arange(1700,3201,50)]

    if wn_range is not None:
        # must include 900.0 (this is required for LBLDIS, as it uses that wn
        # for reference for the geometric limit Q_e rescaling)
        if wn_range[0] > 900.0:
            wn_low = 900.0
        else:
            wn_low = wn_range[0]
        if wn_range[1] < 900.0:
            wn_high = 900.0
        else:
            wn_high = wn_range[1]
        wavenumbers = wavenumbers[np.searchsorted(wavenumbers,wn_low)-1 :
                                  np.searchsorted(wavenumbers,wn_high)+1]
    data = _load_mono_ssp_data(ncdb_file, wavenumbers)

    bulk_data = _compute_bulk_ssp(data, psd_descriptions)
    
    comments = \
        ['Auto generated by lbldis_utils.integrate_ssp_db()', 
         'using netcdf file ' + ncdb_file]
    
    # post processing to populate remaining fields (these are mostly 
    # derived from what already exists)
    ssp_db.write_ssp_db(ssp_db_file, bulk_data, comments)

    if save_npz is not None:
        np.savez(save_npz, bulk_data = bulk_data)

    return data, bulk_data

def _compute_wtd_var(v, wt, x):
    # computed a variance of a quantity (x) with some weighting (wt), 
    # both dependent on variable x, Here the variance is defined as 
    # 2nd integral moment, integ(x-xbar)/integ(wt), where xbar is the 
    # weighted first moment, integ(x*wt).

    m1 = np.trapz(v*wt, x) / np.trapz(wt, x)
    m2 = np.trapz( (v-m1)*(v-m1)*wt, x ) / np.trapz(wt,x)
    return m2


def _compute_s_ext_wt(N, D, s_ext):

    # compute extincton cross section weighting array; given input arrays for:
    # D, N(D), s_ext(D), the diameter, number density, 
    # and extinction cross section as a function of Diameter 

    # this calculation is intended to match the trapezoidal integration 
    # used in other sections (np.trapz); the result is that the 
    # un-normalized sum computed here should match np.trapz(N*s_ext, D)

    dD = np.diff(D)
    beta_e = N * s_ext
    avg_beta_e = 0.5*(beta_e[1:] + beta_e[0:-1])

    s_ext_wt = np.zeros(N.shape)
    s_ext_wt[0:-1] = dD * (avg_beta_e + beta_e[0:-1]) * 0.25
    s_ext_wt[1:] = s_ext_wt[1:] + dD * (avg_beta_e + beta_e[1:]) * 0.25

    #print s_ext_wt.sum(), np.trapz(N*s_ext, D)

    # normalize so the sum of the weights == 1
    s_ext_wt = s_ext_wt / s_ext_wt.sum()

    return s_ext_wt


def _load_mono_ssp_data(ncdb_file, wavenumbers):

    # see comments in _load_lbldis_nc_file for explanation.
    if StrictVersion(scipy_version.version) < StrictVersion('0.8.0'):
        return _load_mono_ssp_data_oldver(ncdb_file, wavenumbers)

    # Note - we are opening with mmap=True (which is actually the default, 
    # but being explicit here to describe the possible issues here)
    # So, later calls to specific variables need to use copy(), otherwise 
    # we are creating many many open file handles if this is called in a loop.

    # the design decision here is to open with mmap, since this specific 
    # function call will only want to load a tiny fraction of the NC file 
    # (it is just one index into the 3000 wavenumbers)

    f = netcdf_file(ncdb_file, 'r', mmap=True)

    db_wn = f.variables['wavenumber'].data.copy()

    actual_wavenumbers, wn_idx = _get_matching_wavenumbers(db_wn, wavenumbers)

    angles = f.variables['scattering_angles'].data.copy().astype('double')
    mu = np.cos(angles*np.pi/180.0)

    D = f.variables['particle_sizes'].data.copy()
    # Note, new versions (the Yang 2012 data conversion) don't have the
    # redundant (2-D) shape information. So check the ndim to determine
    # how to read this variable.
    if f.variables['area'].data.ndim > 1:
        # A,V, are present per wavenumber in the file, but I don't think 
        # they every vary (these are properties of the shape of the crystal), 
        # so just pick the first one.
        A = f.variables['area'].data[0, :].copy()
        V = f.variables['volume'].data[0, :].copy()
    else:
        A = f.variables['area'].data.copy()
        V = f.variables['volume'].data.copy()
    n_size = D.shape[0]
    n_wn = actual_wavenumbers.shape[0]
    n_angle = angles.shape[0]

    Q_ext = np.zeros((n_size, n_wn))
    w = np.zeros((n_size, n_wn))
    g = np.zeros((n_size, n_wn))
    P11 = np.zeros((n_angle, n_size, n_wn))

    for n, wn_i in enumerate(wn_idx):

        Q_ext[:,n] = f.variables['Q_ext'].data[wn_i, :]
        w[:,n] = f.variables['ssa'].data[wn_i, :]
        g[:,n] = f.variables['g'].data[wn_i, :]
        P11[:,:,n] = f.variables['P11'].data[:,wn_i,:]

    f.close()

    return {'wavenumber':actual_wavenumbers, 
            'D':D, 'A':A, 'V':V, 'Q_ext':Q_ext, 
            'g':g, 'w':w, 'P11':P11, 'mu':mu, 'angles':angles}

def _load_mono_ssp_data_oldver(ncdb_file, wavenumbers):

    # Note - there is no memmap option in the old version, but it
    # appears that memmap is used. for example:
    # >>> var = f.variables['something'][:]
    # >>> type(var.base.base)
    #     <type 'mmap.mmap'>
    # So, for each variable, call the copy method so that the mmap references
    # are not kept around.

    f = netcdf_file(ncdb_file, 'r')

    db_wn = f.variables['wavenumber'][:].copy()

    actual_wavenumbers, wn_idx = _get_matching_wavenumbers(db_wn, wavenumbers)

    angles = f.variables['scattering_angles'][:].copy().astype('double')
    mu = np.cos(angles*np.pi/180.0)

    D = f.variables['particle_sizes'][:].copy()
    # Note, new versions (the Yang 2012 data conversion) don't have the
    # redundant (2-D) shape information. So check the ndim to determine
    # how to read this variable.
    if f.variables['area'][:].ndim > 1:
        # A,V, are present per wavenumber in the file, but I don't think 
        # they every vary (these are properties of the shape of the crystal), 
        # so just pick the first one.
        A = f.variables['area'][0,:].copy()
        V = f.variables['volume'][0,:].copy()
    else:
        A = f.variables['area'][:].copy()
        V = f.variables['volume'][:].copy()
    n_size = D.shape[0]
    n_wn = actual_wavenumbers.shape[0]
    n_angle = angles.shape[0]

    Q_ext = np.zeros((n_size, n_wn))
    w = np.zeros((n_size, n_wn))
    g = np.zeros((n_size, n_wn))
    P11 = np.zeros((n_angle, n_size, n_wn))

    for n, wn_i in enumerate(wn_idx):

        # should not need to .copy() here since we are copying slices.
        Q_ext[:,n] = f.variables['Q_ext'][wn_i, :]
        w[:,n] = f.variables['ssa'][wn_i, :]
        g[:,n] = f.variables['g'][wn_i, :]
        P11[:,:,n] = f.variables['P11'][:,wn_i,:]

    f.close()

    return {'wavenumber':actual_wavenumbers, 
            'D':D, 'A':A, 'V':V, 'Q_ext':Q_ext, 
            'g':g, 'w':w, 'P11':P11, 'mu':mu, 'angles':angles}

def _get_matching_wavenumbers(db_wn, wavenumbers):

    # this is complicated by the fact that different netcdf DB's sort 
    # the wavelengths differently

    # first just do a trap to make sure the db encloses the requested 
    # wavenumbers.
    if np.all(wavenumbers[0] < db_wn):
        raise ValueError, 'Requested wavenumber ' + str(wavenumbers[0]) + \
            ' is below smallest wavenumber in db ' + str(db_wn.min())
    if np.all(wavenumbers[-1] > db_wn):
        raise ValueError, 'Requested wavenumber ' + str(wavenumbers[-1]) + \
            ' is above largest wavenumber in db ' + str(db_wn.max())

    # determine actual wavenumbers that can be used; we can use 
    # a Python list here to make the syntax easy

    # first find closest matching wavenumbers in the database array, 
    # taking care to skip duplicates.
    actual_wavenumbers = []
    wn_idx = []
    for n, wn in enumerate(wavenumbers):
        tmp_idx = np.argmin( abs(db_wn - wn) )
        closest_wn = db_wn[tmp_idx]
        if not closest_wn in actual_wavenumbers:
            actual_wavenumbers.append(closest_wn)
            wn_idx.append(tmp_idx)

    # 2nd, make sure that this list encloses the requested list, by potentially 
    # adding end points. If this doesnt work, throw an exception (LBLDIS 
    # will need the database to be "wider" than any requested wavenumbers)
    if wavenumbers[0] < actual_wavenumbers[0]:
        # this looks funny because db_wn could be sorted either small->large 
        # or large->small in wavenumber. We just need to add one additional 
        # wavenumber that is the smaller adjacent wavenumber to the first 
        # one we selected.
        tmp_idx = np.argmin( abs(db_wn - wavenumbers[0]) )
        if db_wn[tmp_idx-1] < db_wn[tmp_idx+1]:
            tmp_idx -= 1
        else:
            tmp_idx += 1
        wn_idx.insert(0,tmp_idx)
        start_wn = db_wn[tmp_idx]
        actual_wavenumbers.insert(0,start_wn)
    if wavenumbers[-1] > actual_wavenumbers[-1]:
        # similar to above, but we want to add one higher wavenumber.
        tmp_idx = np.argmin( abs(db_wn - wavenumbers[-1]) )
        if db_wn[tmp_idx-1] < db_wn[tmp_idx+1]:
            tmp_idx += 1
        else:
            tmp_idx -= 1
        wn_idx.append(tmp_idx)
        end_wn = db_wn[tmp_idx]
        actual_wavenumbers.append(end_wn)
    # convert this to float64 (e.g. normal python float) at the same 
    # time, to prevent downstream errors when we try to print the db.
    actual_wavenumbers = np.array(actual_wavenumbers, dtype=np.float)

    return actual_wavenumbers, wn_idx

def _compute_bulk_ssp(data, psd_descriptions, recompute_g=True):

    # Input: data from load_mono_ssp, which contains D, A, V, N, IWC
    # (n_size elements); Q_ext, g, w, wavenumber (n_size, n_wn elements); 
    # mu, angles (n_angle elements), and P11 (n_angle, n_size, n_wn)

    # output: a PyDict, that can be directly written as a ssp_db.
    # includes all required fields, plus w_wvar, g_wvar, Q_ext_wvar
    # (weighted variances, for testing/debugging), and IWC 
    # (normalized IWC as function of D)

    n_angle, n_size, n_wn = data['P11'].shape

    if type(psd_descriptions) == dict:
        n_r_eff = 1
    else:
        n_r_eff = len(psd_descriptions)

    bulk_data = {}
    bulk_data['wavenum'] = data['wavenumber']
    bulk_data['wavelen'] = 1e4/data['wavenumber']
    bulk_data['angles'] = data['angles']
    bulk_data['r_eff'] = np.zeros(n_r_eff)
    bulk_data['N'] = np.zeros((data['D'].shape[0], n_r_eff))
    bulk_data['IWC'] = np.zeros((data['D'].shape[0], n_r_eff))
    bulk_data['Nt'] = np.zeros(n_r_eff)
    bulk_data['IWCt'] = np.zeros(n_r_eff)

    if recompute_g:
        mu = np.cos(np.pi*data['angles']/180.0)

    ssp_dtype = ssp_db.get_ssp_dtype(n_angle)

    bulk_data['ssp'] = np.zeros( (n_wn, n_r_eff), dtype = ssp_dtype )

    # to simplify syntax
    D = data['D']
    T = np.trapz

    bulk_data['s_ext_wt'] = np.zeros((D.shape[0], n_r_eff))

    # hardcoded const for IWC, g/cm^3
    ice_density = 0.917

    for p, psd in enumerate(psd_descriptions):

        psd_func = _psd_function_list[psd['name']]
        N = psd_func(D, **psd['params'])
        N0 = T(N, D)
        N = N / N0
        total_A = T(data['A']*N, D)
        total_V = T(data['V']*N, D)
        bulk_data['r_eff'][p] = 0.75 * total_V / total_A

        # units conversions - [m^-3/um^-3], [g/m^-3]/[g/cm^-3]
        IWC = data['V'] * N * N0 * ice_density * 1e-18 * 1e6
        bulk_data['N'][:,p] = N
        bulk_data['IWC'][:,p] = IWC
        bulk_data['Nt'][p] = N0
        bulk_data['IWCt'][p] = T(IWC, D)

        for w, wn in enumerate(data['wavenumber']):

            s_ext = data['A'] * data['Q_ext'][:,w]
            s_sca = s_ext * data['w'][:,w]
            total_s_ext = T(s_ext*N, D)
            total_s_sca = T(s_sca*N, D)

            bulk_data['ssp'][w,p]['s_ext'] = T(s_ext*N, D)
            bulk_data['ssp'][w,p]['s_sca'] = T(s_sca*N, D)
            bulk_data['ssp'][w,p]['s_abs'] = \
                bulk_data['ssp'][w,p]['s_ext']-bulk_data['ssp'][w,p]['s_sca']

            bulk_data['ssp'][w,p]['Q_ext'] = total_s_ext / total_A
            bulk_data['ssp'][w,p]['Q_sca'] = total_s_sca / total_A
            bulk_data['ssp'][w,p]['Q_abs'] = \
                bulk_data['ssp'][w,p]['Q_ext']-bulk_data['ssp'][w,p]['Q_sca']

            bulk_data['ssp'][w,p]['w0'] = total_s_sca / total_s_ext

            bulk_data['ssp'][w,p]['vol'] = total_V
            bulk_data['ssp'][w,p]['proj_area'] = total_A

            for m in range(n_angle):
                bulk_data['ssp'][w,p]['phase_func'][m] = \
                    T(s_sca * N * data['P11'][m,:,w], D) / total_s_sca

            if recompute_g:
                P11 = bulk_data['ssp'][w,p]['phase_func']
                bulk_data['ssp'][w,p]['g'] = T(P11*mu, mu) / T(P11, mu)
            else:
                bulk_data['ssp'][w,p]['g'] = \
                    T(s_sca * N * data['g'][:,w], D) / total_s_sca

            s_ext_wt = _compute_s_ext_wt(N, D, s_ext)
            bulk_data['s_ext_wt'][:,p] = s_ext_wt
            # Not sure what this is doing anymore, seems dead code
#            bulk_data['ssp'][w,p]['w_wvar'] = \
#                _compute_wtd_var(data['w'][:,w], s_ext_wt, D)
#            bulk_data['ssp'][w,p]['g_wvar'] = \
#                _compute_wtd_var(data['g'][:,w], s_ext_wt, D)
#            bulk_data['ssp'][w,p]['Q_ext_wvar'] = \
#                _compute_wtd_var(data['Q_ext'][:,w], s_ext_wt, D)

    return bulk_data
