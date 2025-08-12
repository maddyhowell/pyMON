# pyMON

This pipeline estimates numax and the uncertainty from a power spectra. $\nu_{\rm max}$ is measured by from the frequency of the peak amplitude from the smoothed power spectra within an envelope around a predicted numax estimate. Uncertainties are determined by peturbing the power spectrum and measuring $\nu_{\rm max}$ again, and then finding the scatter in the distributions of the measured $\nu_{\rm max}$. This code will also estimate the peak amplitude of the smoothed power spectra and the width (and uncertainties for both using the same method). Smoothing is done by colvolving the power with a Gaussian kernel (based on the nuSYD method). 

This code includes four ways to deal with the background due to granulation:

1. white: This will only remove the white noise from the power spectrum before estimating $\nu_{\rm max}$
2. nuSYD: divides by $(\nu/\nu_{\rm max})^{-2}$ to the power before smoothing. Removes power at lower frequencies (i.e. just removed Granulation noise). See Sreenivas et al. 2024 (https://arxiv.org/abs/2401.17557) for more details
3. linear: estimates a linear model between the power excess envelope (implemented by Simon Campbell). Divides the linear background from power after smoothing. See Howell et al. 2022 (https://academic.oup.com/mnras/article/515/3/3184/6649831?login=false) for more details
4. harvey: fits a harvey-like function to the granulation noise. Divides the background from power after smoothing.

## Usage
`pyMON` requires the following:
- light curve (in units of flux and days)
- power spectrum (in units of ppm^2/muHz and muHz)
- an initial estimate for $\nu_{\rm max}$
- an identification for the star (Star_ID)

Other required inputs to be passed in as a python dictionary include:
- sm: smoothing parameter (same as pySYD). Set to unity for no extra smoothing
- lowerp: the lower frequency of the power excess. If set to None, the lowerp frequency will be determined from a full width half max estimate
- upperp: the lower frequency of the power excess. If set to None, the upperp frequency will be determined from a full width half max estimate
- background_model: chosen from one of the options above. Note: there is a known bug for the harvey background model.
- mc_iters: number of iterations for the mc uncertainty calculation. Set to '0' to not estimate uncertainties.
- Dnu_relation: parameter to determine the $\Delta\nu$-$\nu_{\rm max}$ scaling relation to use, in the form $\Delta\nu = \rm coefficient\times \nu_{\rm max}^ {\rm exponent}$. There are two options for this parameter: i) provide a list in the form ```[coefficient, exponent]```, or ii) provide a string keyword that corresponds to a two element list in the dnu_relations.py file. This file can be edited to include your own relations. For a general relation, use the keyword 'pySYD'.

See code for other default inputs.

To run pyMON
```python
from pyMON.pyMON import pyMON
inputs = {'sm': 1, 'lowerp': None, 'upperp': None, 'background_model': 'linear', 
          'numax_est': 40, 'mc_iters': 500, 'Dnu_relation': 'pySYD'}
pyMON_df = pyMON(star_psd.frequency, star_psd.power, star_lc.time, Star_ID, inputs)
```

This code returns a pandas dataframe with the inputs and measured numax, width and amplitude. A directory would of be created called './results/{Star_ID}/{background_model}' were three plots have been saved and a csv file (same information as in pyMON_df).

## Example Jupyter Notebook
An example is also provided in the jupyter notebook pyMON_example.ipynb. 

## Citing
If you use `pyMON`, please include the following citation (https://academic.oup.com/mnras/article/536/2/1389/7916660)
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

