import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, convolve, Gaussian1DKernel, convolve_fft
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def bin_data(x, y, width = 1, log=False, mode='mean'):
    """
    Bins a series of data.

    Parameters
    ----------
    x : np.ndarray
        the x values of the data
    y : np.ndarray
        the y values of the data
    width : float
        bin width in muHz
    log : bool
        creates bins by using the log of the min/max values (i.e. not equally spaced in log if `True`)

    Returns
    -------
    bin_x : np.ndarray
        binned frequencies
    bin_y : np.ndarray
        binned power
    bin_yerr : numpy.ndarray
        standard deviation of the binned y data

    """
    if log:
        mi = np.log10(min(x))
        ma = np.log10(max(x))
        no = np.int_(np.ceil((ma-mi)/width))
        bins = np.logspace(mi, mi+(no+1)*width, no)
    else:
        bins = np.arange(min(x), max(x)+width, width)

    digitized = np.digitize(x, bins)
    if mode == 'mean':
        bin_x = np.array([x[digitized == i].mean() for i in range(1, len(bins)) if len(x[digitized == i]) > 0])
        bin_y = np.array([y[digitized == i].mean() for i in range(1, len(bins)) if len(x[digitized == i]) > 0])
    elif mode == 'median':
        bin_x = np.array([np.median(x[digitized == i]) for i in range(1, len(bins)) if len(x[digitized == i]) > 0])
        bin_y = np.array([np.median(y[digitized == i]) for i in range(1, len(bins)) if len(x[digitized == i]) > 0])
    else:
        pass
    bin_yerr = np.array([y[digitized == i].std()/np.sqrt(len(y[digitized == i])) for i in range(1, len(bins)) if len(x[digitized == i]) > 0])

    return bin_x, bin_y, bin_yerr

def get_initial(frequency, power, time, max_laws = 1):
    """
    Gets initial guesses for granulation components (i.e. timescales and amplitudes) using
    solar scaling relations. This resets the power spectrum and has its own independent
    filter (i.e. [lower,upper] mask) to use for this subroutine.

    Parameters
    ----------
    star : target.Target
        pySYD target object
    star.oversample : bool
        if `True`, it will use an oversampled power spectrum for the first iteration or 'step'
    minimum_freq : float
        minimum frequency to use for the power spectrum if `None` is provided (via info/star_info.csv). Default = `10.0` muHz. Please note: this is typically sufficient for most stars but may affect evolved stars!
    maximum_freq : float
        maximum frequency to use for the power spectrum if `None` is provided (via info/star_info.csv). Default = `5000.0` muHz.

    Returns
    -------
    star : target.Target
        updated pySYD target object

    """

    resolution = frequency[1] - frequency[0]

    lower_bg = 1
    upper_bg = 274.7
    bg_mask = [lower_bg,upper_bg]

    # Mask power spectrum for fitbg module
    mask = np.ma.getmask(np.ma.masked_inside(frequency, bg_mask[0], bg_mask[1]))
    frequency, power = np.copy(frequency[mask]), np.copy(power[mask])
    random_pow = np.copy(power)


    # Estimate granulation time scales
    tau_sun_single = [3.8e6,2.5e5,1.5e5,1.0e5,230.,70.]        # default from pySYD
    scale = 1.0                                                # default from pySYD
    taus = np.array(tau_sun_single)*scale

    baseline = (max(time)-min(time))*24.*60.*60.

    taus = taus[taus <= baseline]
    b = taus*10**-6.
    mnu = (1.0/taus)*10**5.

    b = b[mnu >= min(frequency)]
    mnu = mnu[mnu >= min(frequency)]


    if len(mnu)==0:
        b = b[mnu >= 10.]
        mnu = mnu[mnu >= 10.]
    elif len(mnu) > max_laws:
        b = b[mnu >= min(frequency)][-max_laws:]
        mnu = mnu[mnu >= min(frequency)][-max_laws:]
    else:

        pass


    return b, mnu

def estimate_initial_red(frequency, power, time, numax_est, ps_mask, w_noise, nlaws = 1):
    """
    Estimates amplitude of red noise components by using a smoothed version of the power
    spectrum with the power excess region masked out. This will take the mean of a specified
    number of points (via -nrms, default=20) for each Harvey-like component.

    Parameters
    ----------
    a : List[float]
        initial guesses for the amplitudes of all Harvey components

    Returns
    -------
    None

    """
    a = []

    box_filter = 1   #default from pySYD
    resolution = frequency[1] - frequency[0]
    n_rms = 20       #default from pySYD

    # Exclude region with power excess and smooth to estimate red noise components
    boxkernel = Box1DKernel(int(np.ceil(box_filter/resolution)))
    mask = (frequency >= ps_mask[0]) & (frequency <= ps_mask[1])

    smooth_pow = convolve(power, boxkernel)

    # Temporary array for inputs into model optimization
    guesses = np.zeros((nlaws*2+1))

    b, mnu = get_initial(frequency, power, time, max_laws = nlaws)

    # Estimate amplitude for each harvey component
    for n, nu in enumerate(mnu):
        diff = list(np.absolute(frequency - nu))
        idx = diff.index(min(diff))

        if idx < n_rms:
            guesses[2*n+1] = np.sqrt((np.mean(smooth_pow[~mask][:n_rms]))/(4.*b[n]))

        elif (len(smooth_pow[~mask])-idx) < n_rms:
            guesses[2*n+1] = np.sqrt((np.mean(smooth_pow[~mask][-n_rms:]))/(4.*b[n]))

        else:
            guesses[2*n+1] = np.sqrt((np.mean(smooth_pow[~mask][idx-int(n_rms/2):idx+int(n_rms/2)]))/(4.*b[n]))

        guesses[2*n] = b[n]
        a.append(guesses[2*n+1])

    guesses[-1] = w_noise
    a = np.array(a)

    numax_sun = 3090.0
    guesses[0] = guesses[0]*numax_sun/numax_est

    return guesses, a


def fit_background(frequency, power, time, numax_est, ps_mask, verbose):
    # Bin power spectrum to model stellar background/correlated red noise components
    bin_freq, bin_pow, bin_err = bin_data(frequency, power)

    # allocated the min and max freq of the power excess envelope


    # Mask out region with power excess
    mask = np.ma.getmask(np.ma.masked_outside(bin_freq, ps_mask[0], ps_mask[1]))
    bin_freq, bin_pow, bin_err = bin_freq[mask], bin_pow[mask], bin_err[mask]

    if verbose:
        print('------------------------------------------------------')
        print('Determining background model:')
        print('PS binned to %d data points'%len(bin_freq))


    # Estimate white noise level
    idx = np.where(frequency > 274.7)[0][0]
    w_noise = np.mean(power[idx:])

    # Get initial guesses for the optimization of the background model
    guesses, a = estimate_initial_red(frequency, power, time, numax_est, ps_mask, w_noise)

    if verbose:
        print('Guesses [tau, sigma, w_noise]: ', guesses)
    # print('Amplitude guesses: ', a)

    return bin_freq, bin_pow, bin_err, guesses, a

def harvey_one(frequency, tau_1, sigma_1, white_noise):
    """
    One Harvey model

    Parameters
    ----------
    frequency : numpy.ndarray
        the frequency array
    tau_1 : float
        timescale of the first harvey component [Ms]
    sigma_1 : float
        amplitude of the first harvey component
    white_noise : float
        the white noise component

    Returns
    -------
    model : numpy.ndarray
        the one-Harvey model

    """

    model = np.zeros_like(frequency)

    model += (4.*(sigma_1**2.)*tau_1)/(1.0+(2.*np.pi*tau_1*frequency)**2.0+(2.*np.pi*tau_1*frequency)**4.0)

    model += white_noise

    return model

def optimise_harvey(frequency,  bin_freq, bin_pow, bin_err, guesses, verbose):

    pars, _ = curve_fit(harvey_one, bin_freq, bin_pow, p0=guesses, sigma=bin_err)

    if verbose:
        print('Optimised parameters [tau, sigma, w_noise]: ', pars)

    model = harvey_one(frequency, pars[0], pars[1], pars[2])

    return model, pars


def test_plot(frequency, power, bin_freq, bin_pow, model, improved_model = []):

    fig, ax = plt.subplots(figsize = (18,10));
    ax.plot(frequency, power, c = 'black', alpha = 0.5, zorder = 0);

    ax.scatter(bin_freq, bin_pow, c = 'r', label = 'Binned PS')

    ax.plot(frequency, model, c = 'b', label = 'Harvey')

    if len(improved_model) != 0:
        ax.plot(frequency, improved_model, c = 'g', label = 'Improved Harvey')



    ax.set_xlim(1,280)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.legend(loc = 'best')
