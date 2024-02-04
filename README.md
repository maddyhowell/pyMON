# GetNumax (name work in progress) or pyMELB (play on pySYD?)

This pipeline estimates numax and the uncertainty from a power spectra. $\nu_{\rm max}$ is measured by from the frequency of the peak amplitude from the smoothed power spectra within an envelope around a predicted numax estimate. Uncertainties are determined by peturbing the power spectrum and measuring $\nu_{\rm max}$ again, and then finding the scatter in the distributions of the measured $\nu_{\rm max}$. This code will also estimate the peak amplitude of the smoothed power spectra and the width (and uncertainties for both using the same method). Smoothing is done by colvolving the power with a Gaussian kernel (based on the nuSYD method). 

This code includes four ways to deal with the background due to granulation (white noise is always subtracted):

1. white: This will only remove the white noise from the power spectrum before estimating $\nu_{\rm max}$
2. nuSYD: divides by $(\nu/\nu_{\rm max})^{-2}$ to the power before smoothing. Removes power at lower frequencies (i.e. just removed Granulation noise). 
3. linear: estimates a linear model between the power excess envelope (same method from pySYD implement by Simon Campbell). Subtracts the linear background from power after smoothing.
4. harvey: fits a harvey-like function to the granulation noise (based on pySYD method). Subtracts the linear background from power after smoothing.

Can access the code by running the jupyter notebook. A prediction for nuamx is needed initially. For M9, M80 & M19, the relation $\nu_{\rm max} = 1.229\times 10^{-19}G_{\rm mag}^{16.89}$. The evolutionary stage is needed as well. The code is set up for metal-poor low mass red giants (RGB, HB, AGB).
