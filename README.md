# pvhub
Function to predict peculiar velocities given RA (right ascension), Dec (declination), and *CMB-frame* redshift. All maps are in redshift-space.
All maps are limited to z < 0.067 with a flag in the function for extrapolation option. Conversion from real-space to redshift space as well as 
the extrapolation option are explained by [Carr et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv211201471C).  
## Maps
The number of each map is the corresponding flag in `pvhub.py`
### default
0. 2M++_SDSS ([Said et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.1275S); [Peterson et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv211003487P); [Carr et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv211201471C)) 
### Other available maps
1. 2M++_SDSS_6dF ([Said et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.1275S))
2. 2MRS ([Lilow & Nusser 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1557L))
3. 2M++ ([Carrick et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..317C))
## Cloning
The PV maps are large files (156 MB each), so to properly clone this repository you must use Git Large File Storage (or download them individually from this webpage). If you have not used Git LFS before, see https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage or follow the instructions below.
### Mac 
    brew install git-lfs
    git lfs install
    git clone https://github.com/KSaid-1/pvhub.git
### Linux
    sudo apt install git-lfs
    git lfs install
    git clone https://github.com/KSaid-1/pvhub.git
### Windows
You can use the [GitHub Desktop](https://desktop.github.com/) GUI application or [Git for Windows](https://git-scm.com/download/win), both of which come with Git LFS. If not using a GUI, then continue as normal in Git Bash:

    git clone https://github.com/KSaid-1/pvhub.git

## Running
Simply `import pvhub` into your python code, execute the `choose_model()` function that accepts an integer from 0 to 3 for the maps as listed above, then run `calculate_pv()` with your RA, Dec and redshift:

    import pvhub
    pvhub.choose_model(flag)
    pv = pvhub.calculate_pv(RA, Dec, zcmb, extrapolation=True)

A few examples of how to use this code are shown in `example.py`. 

You only need to run `choose_model()` once and `calculate_pv()` will continue to use what was selected. 
If you do not choose a model, the default model 0 will be selected.
You can input single objects or lists of RA, Dec, and the redshift to calculate for any number of objects at once.

Coordinates are expected to be in decimal degrees, and output peculiar velocities are in km/s.
We provide an example set of input objects in `inputs/example.csv`.
