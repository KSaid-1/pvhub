# pvhub
Function to predict peculiar velocities given RA (right ascension), Dec (declination), and z (redshift). All maps are in redshift-space.
All maps are limited to z < 0.067 with a flag in the function for extrapolation option. Conversion from real-space to redshift space as well as 
the extrapolation option are explained by [Carr et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021arXiv211201471C).  
## Maps
The number of each map is the corresponding flag in `pvhub.py`
### default
0. 2MPP-SDSS ([Said et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.1275S); [Peterson et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv211003487P); [Carr et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv211201471C)) 
### Other available maps
1. 2MPP-SDSS_6dF ([Said et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.1275S))
2. 2MRS ([Lilow & Nusser 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1557L))
3. 2MPP ([Carrick et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.450..317C))
## Cloning
The PV maps are large files, so to properly clone this repository you must download them individually from this webpage or use Git Large File Storage. If you have not used Git LFS before, see https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage or follow the instructions below.
### Mac 
    brew install git-lfs
    git lfs install
    git lfs clone https://github.com/KSaid-1/pvhub.git
### Linux
    sudo apt install git-lfs
    git lfs install
    git lfs clone https://github.com/KSaid-1/pvhub.git
### Windows
If you have Git Bash, install Git LFS using the [installer](https://git-lfs.github.com), then continue as normal in Git Bash:

    git lfs install
    git lfs clone https://github.com/KSaid-1/pvhub.git

