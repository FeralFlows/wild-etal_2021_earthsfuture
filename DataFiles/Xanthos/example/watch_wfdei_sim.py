"""
Runs Xanthos across4 GCMs over history with WATCH + WFDEI reanalysis product.

5/28/2019

"""
from xanthos import Xanthos
import glob

# The short names of all the gcms to run. (Note that watch+wfdei is not a gcm, it's a bias correction target data set that we use to simulate history).
GCMS = ['watch+wfdei']

# Name of the .ini file containing the configuration to run
ini = "watch_impacts.ini"
base_dir = r'E:/NEXO-UA/Xanthos/example'
comb_yr = '1970_2010'
StartYear = 1970
EndYear = 2010
for gcm in GCMS:
    # Directory containing .npy files for each GCM. File names must contain the
    # GCM and either 'pr' or 'tas'. For example:
    #   pr_bced_1960_1999_ipsl-cm5a-lr_historical_mmpermth_1950_2005.npy
    CLIMATE_DATA_DIR = base_dir + "/input/climate/" + gcm + '/'
    climate_files = glob.glob(CLIMATE_DATA_DIR + "*.npy")
    # Filter to current GCM
    gcm_pr_future = [cf for cf in climate_files if gcm in cf and 'pr' in cf and comb_yr in cf]
    pm_tas = [cf for cf in climate_files if gcm in cf and 'tas_' in cf and comb_yr in cf]
    pm_tmin = [cf for cf in climate_files if gcm in cf and 'tasmin' in cf and comb_yr in cf]
    TempMaxFile = [cf for cf in climate_files if gcm in cf and 'tas_' in cf and comb_yr in cf]
    pm_rhs = [cf for cf in climate_files if gcm in cf and 'rhs' in cf and comb_yr in cf]
    pm_rlds = [cf for cf in climate_files if gcm in cf and 'rlds' in cf and comb_yr in cf]
    pm_rsds = [cf for cf in climate_files if gcm in cf and 'rsds' in cf and comb_yr in cf]
    pm_wind = [cf for cf in climate_files if gcm in cf and 'wind' in cf and comb_yr in cf]

    # Directory for outputting the Xanthos results
    output_dir = base_dir + "/output" + "/clim_impacts" + "_" + gcm + '_' + comb_yr
    args = {
        "OutputFolder": output_dir,
        "OutputNameStr": gcm + '_' + comb_yr,
        "ProjectName": gcm + "_" + comb_yr,
        "PrecipitationFile" : gcm_pr_future[0],
        "pm_tas": pm_tas[0],
        "pm_tmin": pm_tmin[0],
        "TempMinFile": pm_tmin[0],
        "TempMaxFile": TempMaxFile[0],
        "pm_rhs": pm_rhs[0],
        "pm_rlds": pm_rlds[0],
        "pm_rsds": pm_rsds[0],
        "pm_wind": pm_wind[0],
        "StartYear": StartYear,
        "EndYear": EndYear
    }
    xth = Xanthos(ini)
    xth.execute(args)