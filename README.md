# dem_error_study

Code and results of [**Hugonnet et al. (2022), *Uncertainty analysis of digital elevation models by spatial inference from stable terrain***](https://doi.org/10.1109/jstars.2022.3188922). :globe_with_meridians: :mount_fuji: 

Below a short guide to: perform the uncertainty analysis of your own DEMs, retrieve the case study datasets, reproduce the processing steps with the case studies, reproduce the figures and tables of the paper.

![alt text](https://github.com/rhugonnet/dem_error_study/blob/main/figures/fig_2.png?raw=true)

## Uncertainty analysis of your own data with xDEM

The results of this study are based on routines implemented in [xDEM](https://github.com/GlacioHack/xdem) with documentation at https://xdem.readthedocs.io/.
In particular, xDEM includes tutorials on DEM uncertainty analysis to:
- Estimate the heteroscedasticity of elevation errors ([Heteroscedasticity gallery example here](https://xdem.readthedocs.io/en/latest/auto_examples/plot_nonstationary_error.html#sphx-glr-auto-examples-plot-nonstationary-error-py)),
- Standardize elevation differences to use stable as an error proxy ([Standardization gallery example here](https://xdem.readthedocs.io/en/latest/auto_examples/plot_standardization.html#sphx-glr-auto-examples-plot-standardization-py)),
- Estimate spatial correlation of errors and propagate to uncertainties of volume change ([Spatial correlation gallery example here](https://xdem.readthedocs.io/en/latest/auto_examples/plot_vgm_error.html#sphx-glr-auto-examples-plot-vgm-error-py)).

***Note at the date of 08.07.22:** xDEM is still in development (version 0.0.6), and its documentation in construction. 
**Some changes are ongoing or planned for the summer 2022, including:***

- *Construction of an error pipeline that combines all steps,*
- *Automate the retrieval of ICESat-2 data (using [icepyx](https://github.com/icesat2py/icepyx)) to perform local uncertainty analysis*,
- *Streamlining of existing gallery example,*
- *Minor fixes and improvements of routines.*


## Retrieve the case study datasets

The dataset consists of:
1. **Nearly-simultaneous Pléiades–SPOT-6 elevation differences at the Mont-Blanc massif and the Pléiades DEM used as a reference for
alignment and deriving terrain attributes** (.tif, *~200 MB*) at 5 m posting available at [TBC](TBC).
2. **Nearly-simultaneous ASTER–SPOT-5 elevation differences at the Northern Patagonian Icefield** at 30 m posting (.tif, *~50 MB*)
 available at [TBC](TBC).


## Reproduce the processing steps with the case studies

### Setup environment

Most scripts rely on the code assembled in the package [xDEM](https://github.com/GlacioHack/xdem) which in turns relies on [SciKit-GStat](https://github.com/mmaelicke/scikit-gstat).
Some routines also rely on [GSTools](https://github.com/GeoStat-Framework/GSTools). You can rapidly 
install a working environment containing those packages and their dependencies with the 
*environment.yml* file, located at the root of the repository, using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```sh
conda env create -f environment.yml
```

*Note: Due to continuous development changes, xDEM is here set to v0.0.6 to exactly reproduce the processing steps as in the paper.* 

### How to use

Scripts for reproducing the processing steps are located in *case_study_montblanc/* or *case_study_npi/*. Those are generally quite short as they use one-liner routines of xDEM.
Some computations (e.g., simulation of correlated field) are performed only in the *figures/* scripts. In the future, those might be integrated in xDEM.

While all scripts are commented, the details on the functions used are available through **[xDEM's documentation](https://xdem.readthedocs.io/)**,
 **[SciKit-GStat's documentation](https://mmaelicke.github.io/scikit-gstat/)** and **[GSTools's documentation](https://geostat-framework.readthedocs.io/projects/gstools/en/stable/)**.


## Reproduce the figures and tables of the paper

Scripts for reproducing the figures and tables are located in *figures/*. These scripts also depend on the environment 
file `environment.yml` as they rely on [Cartopy](https://scitools.org.uk/cartopy/docs/latest/) and 
[Seaborn](https://seaborn.pydata.org/) in addition to Matplotlib.  In some occasions, the figure scripts duplicate the 
processing steps done in *case_study_montblanc/* or *case_study_npi/* for plotting purposes (e.g., violin plots require 
the full distribution of samples, not only the binned estimates of dispersion).

For plotting figures of your own data, xDEM provides simpler plotting tools of binned data and variograms 
(see [example gallery](https://xdem.readthedocs.io/en/latest/auto_examples/index.html)).


Enjoy! :volcano: 
