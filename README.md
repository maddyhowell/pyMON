# pyMON

This pipeline estimates numax and the uncertainty from a power spectra. $\nu_{\rm max}$ is measured by from the frequency of the peak amplitude from the smoothed power spectra within an envelope around a predicted numax estimate. Uncertainties are determined by peturbing the power spectrum and measuring $\nu_{\rm max}$ again, and then finding the scatter in the distributions of the measured $\nu_{\rm max}$. This code will also estimate the peak amplitude of the smoothed power spectra and the width (and uncertainties for both using the same method). Smoothing is done by colvolving the power with a Gaussian kernel (based on the `pySYD` method). 

This code includes four ways to model the background noise:

1. white: This will only remove the white noise from the power spectrum before estimating $\nu_{\rm max}$
2. nuSYD: divides by $(\nu/\nu_{\rm max})^{-2}$ to the power before smoothing. Removes power at lower frequencies (i.e. just removed Granulation noise). See [Sreenivas et al. 2024](https://arxiv.org/abs/2401.17557) for more details. **Note: if the power spectrum is noisy, the width measurement can fail. This will affect the numax measurement using this model**
3. linear: estimates a linear model between the power excess envelope (implemented by Simon Campbell). Divides the linear background from power after smoothing. See [Howell et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.515.3184H/abstract) for more details
4. harvey: fits a harvey-like function to the granulation noise. Divides the background from power after smoothing. **Note: harvey-like function is not fully implemented**

`pyMON` was developed to use as a back-up option to `pySYD` for when the signal-to-noise of the star's power spectrum is low and delta nu cannot be easily measured. We recommend using the `linear` model for these cases. 

## Usage

Run code from the same directory as where it is located.

`pyMON` requires the following inputs:
- light curve (in units of flux and days)
- power spectrum (in units of ppm^2/muHz and muHz)
- an initial estimate for $\nu_{\rm max}$
- an identification for the star (Star_ID)

Other required inputs to be passed in as a python dictionary include:
- sm: smoothing parameter (same as pySYD). Default will be 2, to ensure sufficient smoothing for low numax stars.
- lowerp: the lower frequency of the power excess. If set to None, the lowerp frequency will be determined from a full width half max estimate
- upperp: the lower frequency of the power excess. If set to None, the upperp frequency will be determined from a full width half max estimate
- background_model: chosen from one of the options above. Note: there is a known bug for the harvey background model.
- mc_iters: number of iterations for the mc sampling uncertainty routine. Set to '0' to not estimate uncertainties.
- Dnu_relation: parameter to determine the $\Delta\nu$-$\nu_{\rm max}$ scaling relation to use, in the form $\Delta\nu = \rm coefficient\times \nu_{\rm max}^ {\rm exponent}$. There are two options for this parameter: i) provide a list in the form ```[coefficient, exponent]```, or ii) provide a string keyword that corresponds to a two element list in the dnu_relations.py file. This file can be edited to include your own relations. For a general relation derived from 16,000 Kepler stars, use the key 'Yu18' (see [Yu et al., 2018](https://arxiv.org/abs/1802.04455) for details)

See code for other default inputs.

Run `pyMON` in the same directory as the code, via the following:
```python
from pyMON import pyMON
inputs = {'sm': 2, 'lowerp': None, 'upperp': None, 'background_model': 'linear', 
          'numax_est': 40, 'mc_iters': 500, 'Dnu_relation': 'Yu18'}
pyMON_df = pyMON(star_psd.frequency, star_psd.power, star_lc.time, Star_ID, inputs)
```

This code returns a pandas dataframe with the inputs and measured numax, width and amplitude. A directory will be created called './results/{Star_ID}/{background_model}' where three plots have been saved and a csv file (same information as in pyMON_df).

## Example Jupyter Notebook
An example for red giant star KIC 2707716 using the KEPSEISMIC light curve and power spectrum is provided in the jupyter notebook pyMON_tutorial.ipynb. 

## Citing
If you use `pyMON`, please include the following citation [Howell et al., 2025](https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.1389H/abstract)
```tex
@ARTICLE{Howell-2025,
       author = {{Howell}, Madeline and {Campbell}, Simon W. and {Kalup}, Csilla and {Stello}, Dennis and {De Silva}, Gayandhi M.},
        title = "{Asteroseismic masses of red giants in the Galactic Globular Clusters M9 and M19}",
      journal = {\mnras},
     keywords = {Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Astrophysics of Galaxies},
         year = 2025,
        month = jan,
       volume = {536},
       number = {2},
        pages = {1389-1407},
          doi = {10.1093/mnras/stae2686},
archivePrefix = {arXiv},
       eprint = {2412.01089},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025MNRAS.536.1389H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
AND include that " `pyMON` is adapted from `pySYD`" with the following cititation [Chontos et al., 2022](https://joss.theoj.org/papers/10.21105/joss.03331)
```tex
@ARTICLE{2022JOSS....7.3331C,
       author = {{Chontos}, Ashley and {Huber}, Daniel and {Sayeed}, Maryum and {Yamsiri}, Pavadol},
        title = "{pySYD: Automated measurements of global asteroseismic parameters}",
      journal = {The Journal of Open Source Software},
     keywords = {Python, fundamental stellar properties, solar-like oscillations, stellar oscillations, stellar astrophysics, asteroseismology, astronomy, global asteroseismology, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2022,
        month = nov,
       volume = {7},
       number = {79},
          eid = {3331},
        pages = {3331},
          doi = {10.21105/joss.03331},
archivePrefix = {arXiv},
       eprint = {2108.00582},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022JOSS....7.3331C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```


If you use the `nuSYD` model, please include the following citation [Sreenivas et al., 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.3477S/abstract)
```tex
@ARTICLE{2024MNRAS.530.3477S,
       author = {{Sreenivas}, K.~R. and {Bedding}, Timothy R. and {Li}, Yaguang and {Huber}, Daniel and {Crawford}, Courtney L. and {Stello}, Dennis and {Yu}, Jie},
        title = "{A simple method to measure {\ensuremath{\nu}}$_{max}$ for asteroseismology: application to 16 000 oscillating Kepler red giants}",
      journal = {\mnras},
     keywords = {asteroseismology, stars: late-type, stars: oscillations, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2024,
        month = may,
       volume = {530},
       number = {3},
        pages = {3477-3487},
          doi = {10.1093/mnras/stae991},
archivePrefix = {arXiv},
       eprint = {2401.17557},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024MNRAS.530.3477S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Known Issues that are Under Development

Width measurement of the power excess: if you have noisy data, `pyMON` might be fail to identify the position of the FWHM on the left side, and hence will not measure an accurate width. This does not affect the measurement of numax or its uncertainty (unless you are using the `nuSYD` model)

Harvey-like function model: this currently doesn't work on low numax stars. Further tests are needed.
