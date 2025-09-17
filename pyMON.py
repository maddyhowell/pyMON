import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, convolve, Gaussian1DKernel, convolve_fft
from scipy.signal import find_peaks
from tqdm import tqdm
from astropy.stats import mad_std
import os
from . import harvey_model as hm
from . import dnu_relations as dr

def ps_no_wnoise(frequency, power, time, Star_ID, verbose):

    # determine nyquist frequency from light curve
    dt = np.median(np.diff(time))  # median handles uneven sampling better
    nyquist_freq = 0.5 / (dt*24*60*60) *10**6 # units of muHz


    nu_lower = nyquist_freq - nyquist_freq/10

    if nu_lower > frequency[-1]:
        nu_lower = frequency[-1] - frequency[-1]/10

    idx = np.where(frequency > nu_lower)[0][0]
    w_noise = np.mean(power[idx:])

    if verbose:
        print(f'Parameters for {Star_ID}:')
        print('White noise metric= ', w_noise, ' ppm^2/muHz')

    power_no_wnoise = power - w_noise

    return power_no_wnoise, w_noise

def get_PS_mask(numax_est, lowerp, upperp):
    """ Get the boundaries of the power excess around the predicted numax """
    numax_sun = 3090.00
    width_sun = 1300.0
    #Get power excess mask
    width = width_sun*(numax_est/numax_sun)
    maxpower = [numax_est-(width), numax_est+(width)]

    if lowerp == None:
        lowerp = maxpower[0]

    if upperp == None:
        upperp = maxpower[1]

    ps_mask = [lowerp, upperp]

    return ps_mask

def get_linear_bg(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm, lowerp, upperp):
    """
    Estimate the BG with a straight line. Linear in log!

    """
    ps_mask = get_PS_mask(numax_est, lowerp, upperp)
    smps = ps_smooth(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm, lowerp, upperp)




    #---Get vaules for linear fit:
    # print('PSMASK: ', ps_mask[0], ps_mask[1])
    mask = np.ma.getmask(np.ma.masked_inside(frequency, ps_mask[0], ps_mask[1]))
    smps = smps[mask]
    pwr1 = smps[0]
    pwr2 = smps[len(smps)-1]  #SWC: Get power values at start and end of excess region.
    # print('LINEARBG: pwr1, pwr2=',pwr1,pwr2)

    nu1 = ps_mask[0]
    nu2 = ps_mask[1]
    nu1 = np.log10(nu1)
    nu2 = np.log10(nu2)
    # print('LINEARBG: nu1,nu2=',nu1,nu2)


    offsetfrac = 0.10

    pwr1 = pwr1 - offsetfrac*pwr1
    pwr2 = pwr2 - offsetfrac*pwr1  #Move line down a bit.
    pwr1 = np.log10(pwr1)
    pwr2 = np.log10(pwr2)

    grad = (pwr2-pwr1)/(nu2-nu1)
    c = pwr1 - grad*nu1

    # print(f'Parameters of Linear Bkg model: grad = {grad}, c = {c}')


    # print('LINEARBG: offsetfrac, pwr1,pwr2,grad,c=',offsetfrac,pwr1,pwr2,grad,c)

    linbg = np.zeros(len(frequency))

    for i in range(len(frequency)):
        linbg[i] = 10.0**(grad*np.log10(frequency[i]) + c)
        #print('Linear BG, freq, i=',linbg[i],self.frequency[i],i)

    return linbg, grad, c

def ps_no_slope(frequency, power, time, numax_est, dnu_coefficient, dnu_exponent, background_model, sm, lowerp, upperp, w_noise, verbose):
    """ This function should be executed after ps_no_wnoise
        Currently three background models are included: nuSYD (nu/numax)^-2 or a single harvey-like function or linear over FWHM of excess
    """


    if background_model == 'nuSYD':

        bg = (frequency/numax_est)**-2

        power_slope_removed = power/bg

    elif background_model == 'harvey':
        ps_mask = get_PS_mask(numax_est, lowerp, upperp)

        bin_freq, bin_pow, bin_err, guesses, a = hm.fit_background(frequency, power, time, numax_est, ps_mask, verbose)

        model = hm.harvey_one(frequency, guesses[0], guesses[1], guesses[2])

        model2, pars = hm.optimise_harvey(frequency,  bin_freq, bin_pow, bin_err, guesses, verbose)

        if verbose:
            hm.test_plot(frequency, power, bin_freq, bin_pow, model, improved_model = model2)

        pssm = ps_smooth(frequency, power, numax_est, Dnu_relation, sm, lowerp, upperp)

        bg = model2

        power_slope_removed = pssm - bg

    elif background_model == 'linear':
        bg, grad, c = get_linear_bg(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm, lowerp, upperp)

        pssm = ps_smooth(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm,lowerp, upperp)

        power_slope_removed = pssm - bg

    elif background_model == 'white':
        bg = np.zeros(len(frequency))
        bg[:] = w_noise
        # print(bg)
        power_slope_removed = ps_no_wnoise(frequency, power, time, Star_ID, verbose)

    # power_slope_removed = power/bg

    return power_slope_removed, bg

def boxcar_filter(frequency, power):

    resolution = frequency[1] - frequency[0]

    box_filter = 1
    boxkernel = Box1DKernel(int(np.ceil(box_filter/resolution)))
    smooth_pow = convolve(power, boxkernel)

    return smooth_pow

def ps_smooth(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm, lowerp, upperp):
    """ This function should be executed after ps_no_slope """

    dnu_est = dnu_coefficient*(numax_est)**dnu_exponent

    resolution = frequency[1] - frequency[0]

    ## following code taken from pysyd target function line 693
    if sm == None:
        numax_sun = 3090.00
        sm = 4.*(numax_est/numax_sun)**0.2
    if sm < 2.:   # update on 16th sept 2025 to address issues of measuring numax at low frequencies. the smoothing condition needs to be changed to a minimum of 2 Dnu for low numax stars
            sm = 2.


    sig = (sm*(dnu_est/resolution))/np.sqrt(8.0*np.log(2.0))
    # print(sig)   # sig is same as pySYD
    pssm = convolve_fft(np.copy(power), Gaussian1DKernel(int(sig)))


    return pssm

def find_parameters(frequency, power, numax_est, lowerp, upperp, verbose):
    """ This function should be executed after ps_smooth
        I have altered the range to search for the peak of the smoothed ps to +/-0.5*numax_est

    """

    if lowerp != None:
        ii = np.where(frequency > lowerp)[0][0]
    else:
        ii = np.where(frequency > numax_est - 0.5*numax_est)[0][0]

    if upperp != None:
        jj = np.where(frequency < upperp)[0][-1]
    else:
        jj = np.where(frequency < numax_est + 0.5*numax_est)[0][-1]

    idx_numax = np.argmax(power[ii:jj])

    freq_window = frequency[ii:jj]
    power_window = power[ii:jj]

    numax = freq_window[idx_numax]

    # peak power measured as the amplitude of the smoothed power spectra
    peak_power = power_window[idx_numax]

    #  for bad PS, the following line might not work.
    if len(np.where(np.logical_and(power < peak_power/2, frequency > numax))[0]) == 0:
        nn, mm = 0,0
        # print('Cannot find the width!')
    elif frequency[np.where(np.logical_and(power < peak_power/2, frequency > numax))[0][0]] > numax_est + 0.5*numax_est:
    ### if the second point doesn exist within  then just extend the left side
        if verbose:
            print('Right hand side of width cannot be found')

        nn = np.where(power > peak_power/2)[0][0]

        freq_rightside = numax + (numax-frequency[nn])

        mm = np.where(frequency > freq_rightside)[0][0]
    else:
        nn, mm = np.where(power > peak_power/2)[0][0], np.where(np.logical_and(power < peak_power/2, frequency > numax))[0][0]

    width_idx = [mm, nn]

    return numax, peak_power, width_idx

def plot_power_excess(ax, frequency, power, pssm, numax_est, numax_measured, ps_mask, width_idx):

    # factor of 1.8 comes from the aspect ratio of the figure
    left, bottom, width, height = [0.75, 0.56, 0.22, 0.22*1.8]
    axin = ax.inset_axes([left, bottom, width, height])

    mask = np.ma.getmask(np.ma.masked_inside(frequency, ps_mask[0], ps_mask[1]))
    mask_frequency = frequency[mask]
    mask_power = power[mask]
    mask_pssm = pssm[mask]

    # ax.plot(mask_frequency, mask_power, c = 'black', alpha = 1, zorder = 0)

    axin.plot(mask_frequency, mask_pssm, c = 'mediumpurple', alpha = 0.8, zorder = 1)


    axin.axvline(numax_measured, color = 'mediumpurple', ls = 'dashed', lw = 3, zorder = 5)

    axin.axvline(numax_est, color = 'black', ls = 'solid', lw = 3, zorder = 5)

    ii = np.where(frequency > numax_measured)[0][0]
    axin.scatter(frequency[ii], pssm[ii], c = 'mediumpurple', s = 100, zorder = 5)


    xticks = np.linspace(int(ps_mask[0]), int(ps_mask[1]), 5)

    axin.set_xticks(xticks, )
    axin.set_xlim(int(ps_mask[0]), int(ps_mask[1]))
    peaks, _ = find_peaks(mask_power)

    axin.set_title('Smoothed Power Excess')

    axin.scatter(frequency[width_idx[0]],pssm[width_idx[0]], c = 'blue', s = 100)
    axin.scatter(frequency[width_idx[1]],pssm[width_idx[1]], c = 'blue', s = 100)
    axin.plot([frequency[width_idx[0]], frequency[width_idx[1]]], [pssm[width_idx[0]], pssm[width_idx[1]]],
            c = 'blue', ls = 'dashed', zorder = 5)

def plot_PS(scale, frequency, power, pssm, bg, numax_measured, numax_est, ps_mask, Star_ID, peak_measured, width_idx, width_measured, save, background_model):

    fig, ax = plt.subplots(figsize = (18,10));
    ax.plot(frequency, power, c = 'black', alpha = 0.5, zorder = 0);

    if background_model == 'linear' or background_model == 'harvey':
        # To plot for a linear background, the smoothed curved needs to add background to match pySYD

        ax.plot(frequency, pssm, c = 'mediumpurple', alpha = 1, label = 'Smoothed ps', zorder = 1, lw = 3);
    else:
        ax.plot(frequency, pssm, c = 'mediumpurple', alpha = 1, label = 'Smoothed ps - bg', zorder = 1, lw = 3);


    mask = np.ma.getmask(np.ma.masked_inside(frequency, ps_mask[0], ps_mask[1]))
    mask_freq = frequency[mask]
    mask_power = power[mask]

    ax.plot(mask_freq, mask_power, c = 'black', alpha = 1, zorder = 0);


    ax.plot(frequency, bg, c = 'r', ls = 'solid', alpha = 0.8, label = 'Background model')

    # ps_mask = get_PS_mask(numax_est, lowerp, upperp)
    ax.axvline(ps_mask[0], c = 'orange', ls = 'dashed')
    ax.axvline(ps_mask[1], c = 'orange', ls = 'dashed', label = 'Power excess envelope')


    if scale == 'linear':
        peaks, _ = find_peaks(power[mask])
        if numax_est > 40:
            ax.set_xlim(1,280)
        else:
            ax.set_xlim(1,50)

        ax.set_ylim(0,np.max(power[mask][peaks]) + 10**int(np.log10(np.max(power[mask][peaks]))))
    else:
        peaks, _ = find_peaks(power)
        ax.set_xlim(1,280)
        oom_min = np.floor(np.log10(np.min(power[peaks])))
        oom_max = np.floor(np.log10(np.max(power[peaks])))
        # ax.set_ylim(np.min(power[peaks]) + 10**(oom_min-1),np.max(power[peaks]) + 10**(oom_max+1))


    peaks, _ = find_peaks(power)

    ax.axvline(numax_measured, color = 'mediumpurple', label = r'Measured $\nu_{\rm max}$ = ' + str(round(numax_measured,2)) + r' $\mu$Hz',
                ls = 'dashed', lw = 3, zorder = 5)

    ax.axvline(np.nan, color = 'black', ls = 'solid', lw = 3, zorder = 5, label = r'Predicted $\nu_{\rm max}$ = ' + str(round(numax_est,2)) + r' $\mu$Hz')

    ax.scatter(np.nan,np.nan, c = 'blue', s = 100)
    ax.scatter(np.nan,np.nan, c = 'blue', s = 100)
    ax.plot(np.nan, np.nan,
            c = 'blue', ls = 'dashed', label = f'Width = {round(width_measured,2)}' + r' $\mu$Hz', zorder = 5)


    ii = np.where(frequency > numax_measured)[0][0]
    ax.scatter(np.nan, np.nan, c = 'mediumpurple', label = r'$P_{\rm peak} =$ '+ str(round(peak_measured,0)) + r'$\rm ppm^2$/$\mu$Hz',
               s = 100, zorder = 5)


    ax.set_ylabel(r'Power [$\rm ppm^2$/$\mu$Hz]', fontsize = 20);
    ax.set_xlabel(r'Frequency [$\mu$Hz]', fontsize = 20);

    ax.set_yscale(scale)
    ax.set_xscale(scale)
    ax.tick_params(axis='both', which='major', labelsize=18)


    if scale == 'linear':
        plot_power_excess(ax, frequency, power, pssm-bg, numax_est, numax_measured, ps_mask, width_idx)

    ax.legend(loc = 'upper left', fontsize = 15)


    if save:
        if not os.path.isdir(f'./results/{Star_ID}/{background_model}'):
            os.mkdir(f'./results/{Star_ID}/{background_model}')

        plt.savefig(f'./results/{Star_ID}/{background_model}/{Star_ID}_ps_{background_model}_bkg_{scale}.png')

    if scale == 'log':
        plt.close()

def MeasureNumax(Star_ID, frequency, power, time, inputs, verbose, plot, save):
    """ Code to run nusyd and functions above
    numax_est is an estimate for numax
    returns measurement for numax """

    sm = inputs['sm']
    lowerp = inputs['lowerp']
    upperp = inputs['upperp']
    background_model = inputs['background_model']
    numax_est = inputs['numax_est']

    Dnu_relation = inputs['Dnu_relation']
    if type(Dnu_relation) == str:
        ## go to file and find Dnu_relation coefficient and exponent
        dnu_coefficient, dnu_exponent = dr.Find_Dnu_relations(Dnu_relation)
    else:
        dnu_coefficient, dnu_exponent = Dnu_relation

    _____, w_noise = ps_no_wnoise(frequency, power, time, Star_ID, verbose)

    # print(f'Using the {background_model} background model')
    power_slope_removed, bg = ps_no_slope(frequency, power, time, numax_est, dnu_coefficient, dnu_exponent, background_model, sm, lowerp, upperp, w_noise, verbose)


    if background_model == 'nuSYD' or background_model == 'white':
        ## nuSYD smoothes after background correction. Linear is before
        pssm = ps_smooth(frequency, power_slope_removed, numax_est, dnu_coefficient, dnu_exponent, sm, lowerp, upperp)
    else:
        pssm = power_slope_removed



    numax_measured, peak_measured, width_idx = find_parameters(frequency, pssm, numax_est, lowerp, upperp, verbose)
    if width_idx == [0,0]:
        width_measured = np.nan
    else:
        width_measured = frequency[width_idx[0]] - frequency[width_idx[1]]

    if background_model == "nuSYD":
        if verbose:
            print('-------------------------------------------------------------')
            print('Correction applied from Sreenivas et al. 2024 for nuSYD model!')
            print('-------------------------------------------------------------')
            print('\n')
        sigma_actual = 0.248*((width_measured)**1.08)
        numax_measured =  (numax_measured**2 - 2*(sigma_actual)**2)/numax_measured

    if plot:
        ps_mask = get_PS_mask(numax_est, lowerp, upperp)

        if background_model == 'linear':
            pssm = ps_smooth(frequency, power, numax_est, dnu_coefficient, dnu_exponent, sm,lowerp, upperp)

        plot_PS('linear', frequency, power, pssm, bg, numax_measured, numax_est, ps_mask, Star_ID, peak_measured, width_idx, width_measured, save, background_model)
        plot_PS('log', frequency, power, pssm, bg, numax_measured, numax_est, ps_mask, Star_ID, peak_measured, width_idx, width_measured, save, background_model)

    if verbose:
        print('Measured numax = ', numax_measured)
        print('Measured Amplitude = ', peak_measured)
        print('Measured width = ', width_measured)


    return numax_measured, peak_measured, width_measured

def MeasureNumaxUncertainty(Star_ID, frequency, power, time, numax_measured, peak_measured, width_measured, inputs, verbose, plot, save):
    """ Calculates the uncertainty on numax using the same method as pySYD
        params includes [numax, peak, width]"""

    sm = inputs['sm']
    lowerp = inputs['lowerp']
    upperp = inputs['upperp']
    background_model = inputs['background_model']
    mc_iters = inputs['mc_iters']

    # determine nyquist frequency from light curve
    dt = np.median(np.diff(time))  # median handles uneven sampling better
    nyquist_freq = 0.5 / (dt*24*60*60) *10**6 # units of muHz

    bg_mask = [1, nyquist_freq]
    i = 0
    numax_sampling, peak_sampling, width_sampling = [], [], []
    numax_temp = numax_measured

    while i < mc_iters:
        if i == 0:
            # Switch to critically-sampled PS if sampling
            ## mhow: I dont know how that critically samples the array (just creates a mask between 1 to nyquist_freq freq)
            mask = np.ma.getmask(np.ma.masked_inside(frequency, bg_mask[0], bg_mask[1]))

            frequency_cp, power_cp = np.copy(frequency[mask]), np.copy(power[mask])
            resolution = frequency_cp[1] - frequency_cp[0]
            random_pow = (np.random.chisquare(2, len(frequency_cp))*power_cp)/2.0

            # update numax values in inputs
            inputs2 = inputs
            inputs2['numax_est'] = numax_temp

            numax_temp, peak_temp, width_temp = MeasureNumax(Star_ID, frequency_cp, random_pow, time, inputs2, verbose = False, plot = False, save = False)

            numax_sampling.append(numax_temp)
            peak_sampling.append(peak_temp)
            width_sampling.append(width_temp)


            print('------------------------------------------------------\nRunning sampling routine:')
            pbar = tqdm(total = mc_iters)
            pbar.update(1)
        else:
            random_pow = (np.random.chisquare(2, len(frequency_cp))*power_cp)/2.0
            pbar.update(1)

            inputs2 = inputs
            inputs2['numax_est'] = numax_temp

            numax_temp, peak_temp, width_temp = MeasureNumax(Star_ID, frequency_cp, random_pow, time, inputs2, verbose = False, plot = False, save = False)

            numax_sampling.append(numax_temp)
            peak_sampling.append(peak_temp)
            width_sampling.append(width_temp)

        i += 1


    numax_uncertainty = mad_std(numax_sampling)
    if numax_uncertainty == 0:
        print('MAD std didn\'t work. Standard error used instead')
        numax_uncertainty =np.std(numax_sampling)/np.sqrt(mc_iters)
    peak_uncertainty = mad_std(peak_sampling, ignore_nan = True)
    width_uncertainty = mad_std(width_sampling, ignore_nan = True)

    if verbose:
        print(f'Parameters for {Star_ID}:')
        print(f'Nyquist frequency (from light curve) =', nyquist_freq, 'muHz')
        print('Measured numax = ', numax_measured, '+/-', numax_uncertainty, 'muHz (', numax_uncertainty/numax_measured*100, '%)')
        print('Measured Amplitude = ', peak_measured, '+/-', peak_uncertainty, 'ppm^2/muHz')
        print('Measured width = ', width_measured, '+/-', width_uncertainty, 'muHz')



    if plot:
        fig, (ax, ax1, ax2) = plt.subplots(1, 3, figsize = (25,8));

        ax.hist(numax_sampling, color = 'mediumpurple', bins=20, histtype='step', lw=3, facecolor='0.75');
        ax.hist(numax_sampling, color = 'mediumpurple', bins=20, alpha = 0.6);

        ax.axvline(numax_measured, c = 'k', ls = 'dashed', lw = 2)

        ax1.hist(peak_sampling, color = 'mediumpurple', bins=20, histtype='step', lw=3, facecolor='0.75');
        ax1.hist(peak_sampling, color = 'mediumpurple', bins=20, alpha = 0.6);

        ax1.axvline(peak_measured, c = 'k', ls = 'dashed', lw = 2)

        ax2.hist(width_sampling, color = 'mediumpurple', bins=20, histtype='step', lw=3, facecolor='0.75');
        ax2.hist(width_sampling, color = 'mediumpurple', bins=20, alpha = 0.6);

        ax2.axvline(width_measured, c = 'k', ls = 'dashed', lw = 2, label = 'Measured parameter')



        ax.set_xlabel(r'Sampled $\nu_{\rm max}$ [$\mu$Hz]', fontsize = 20);
        ax1.set_xlabel(r'Sampled Amplitude [ppm$^2/\mu$Hz]', fontsize = 20);
        ax2.set_xlabel(r'Sampled width [$\mu$Hz]', fontsize = 20);

        ax2.legend(loc = 'best')
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax1.tick_params(axis='both', which='major', labelsize=18)
        ax2.tick_params(axis='both', which='major', labelsize=18)


        ax.set_yticks([]), ax1.set_yticks([]), ax2.set_yticks([]);
        ax.set_yticklabels([]), ax1.set_yticklabels([]), ax2.set_yticklabels([]);


        if save:
            if not os.path.isdir(f'./results/{Star_ID}/{background_model}'):
                os.mkdir(f'./results/{Star_ID}/{background_model}')

            plt.savefig(f'./results/{Star_ID}/{background_model}/{Star_ID}_samples_{background_model}_bkg.png')




    return numax_uncertainty, peak_uncertainty, width_uncertainty

def pyMON(frequency, power, time, Star_ID, inputs, verbose = True, plot = True, save = True):
    """ Main pyMON function. Calculates the numax of a star.

    Inputs
    -------
    frequency : array-lie
        Frequency values in units of muHz
    power : array-like
        Power values in units of ppm^2/muHz
    time : array-like
        Time values in units of days
    Star_ID : int or str
        Identificaiton for target star. Used for directory names
    inputs : dictionary
        Dictionary of input parameters. Must include {'sm': value, 'lowerp': value, 'upperp': value, 'background_model': value, 'numax_est': value, 'mc_iters': value, 'Dnu_relation': value}
    verbose : bool
        If True, print results. Default is True
    plot : bool
        If True, plot power spectrum figure. Default is True
    save : bool
        If True, save data to a csv file and three figures (power spectrum in a linear and log scale, distributions of the mc sampling) to a directory './results/{Star_ID}/{background_model}/'. Default is True.

    Returns
    -------
    pyMON_df: pandas dataframe with the inputs and outputs of the pyMON function

    """

    file_path = f'./results/{Star_ID}/'

    if not os.path.isdir('./results'):
        os.mkdir('./results')

    if not os.path.isdir(file_path):
        os.mkdir(file_path)


    numax_measured, peak_measured, width_measured = MeasureNumax(Star_ID, frequency, power, time, inputs, verbose, plot, save)

    mc_iters = inputs['mc_iters']
    if mc_iters != 0:
        numax_uncertainty, peak_uncertainty, width_uncertainty = MeasureNumaxUncertainty(Star_ID, frequency, power, time, numax_measured, peak_measured, width_measured, inputs, verbose, plot, save)
    else:
        numax_uncertainty, peak_uncertainty, width_uncertainty = None, None, None


    pyMON_dict = {'Star_ID': [Star_ID], 'numax_predicted': [inputs['numax_est']], 'Dnu_relation': [inputs['Dnu_relation']],
              'lowerp': [inputs['lowerp']], 'upperp': [inputs['upperp']], 'sm': [inputs['sm']], 'background_model': [inputs['background_model']],
              'numax_measured': [numax_measured], 'numax_uncertainty': [numax_uncertainty],
              'peak_measured': [peak_measured], 'peak_uncertainty': [peak_uncertainty],
              'width_measured': [width_measured], 'width_uncertainty': [width_uncertainty], 'mc_iters': [inputs['mc_iters']]
             }

    pyMON_df = pd.DataFrame(pyMON_dict)

    if save:
        if not os.path.exists(file_path+f'{Star_ID}_pyMON_results.csv'):
            # If the file doesn't exist, write with header and no index
           pyMON_df.to_csv(file_path+f'{Star_ID}_pyMON_results.csv', mode='a', index=False, header=True)
        else:
            # If the file exists, append without header and no index
            pyMON_df.to_csv(file_path+f'{Star_ID}_pyMON_results.csv', mode='a', index=False, header = False)

    return pyMON_df
