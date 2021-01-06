# ArgentinaNexus
<!--your zenodo badge here-->

Quantifying the value of emulation for estimating drought statistics

<br />

<!-------------------------->
<!-------------------------->
## Contents
<!-------------------------->
<!-------------------------->
- [Abstract](#abstract)
- [Code Reference](#code-reference)
- [Journal Reference](#journal-reference)
- [Contributing Models](#contributing-models)
- [Data Reference](#data-reference)
- [Reproduce My Experiement](#reproduce-my-experiement)

<br />

<!-------------------------->
<!-------------------------->
## Abstract
<!-------------------------->
<!-------------------------->

Your abstract here.

[Back to Contents](#contents)

<br />

<!-------------------------->
<!-------------------------->
## Code Reference
<!-------------------------->
<!-------------------------->

<!--References for each minted software release for all code involved.  If you have modified a codebase that is outside of a formal release, and the modifications are not planned on being merged back into a version, fork the parent repository and add a `.<shortname>` to the version number of the parent and conduct your own name.  For example, `v1.2.5.hydro`.-->
[1] **Metis:** Khan, Z., Wild, T., Vernon, C., Miller, A., Hejazi, M., Clarke, L., Miralles-Wilhelm, F., Castillo, R.M., Moreda, F., Bereslawski, J.L., Suriano, M. and Casado, J., (2020). Metis v1.1.0. Github. [![GitHub tag](https://img.shields.io/github/v/release/JGCRI/metis)](https://github.com/JGCRI/metis/releases/tag/v1.0.0)


[2] **GCAM:** Khan, Zarrar. (2019, November 5). gcam-v5.1.3LAC_khan_et_al_2020_Uruguay (Version 5.3.1LAC). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3897519.svg)](https://doi.org/10.5281/zenodo.3897519)

[3] **Xanthos:** Braun Caleb, Vernon Chris, Link Robert, Evanoff Jason, & Khan Zarrar. (2020, December 30). xanthos-v2.3.1 for Wild_et_al_2020_ArgentinaNexus (Version v2.3.1-wild2020-ArgentinaNexus). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4404834.svg)](https://doi.org/10.5281/zenodo.4404834)

[4] **Demeter:** Vernon Chris, & Braun Caleb. (2020, December 30). demeter-v1.1.0 for Wild_et_al_2020_ArgentinaNexus (Version v1.1.0-wild2020-ArgentinaNexus). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4404738.svg)](https://doi.org/10.5281/zenodo.4404738)

[5] **Tethys:** Vernon Chris, Link Robert, & Braun Caleb. (2020, December 31). tethys-v1.2.0 for Wild_et_al_2020_ArgentinaNexus (Version v1.2.0-wild2020-ArgentinaNexus). Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4405008.svg)](https://doi.org/10.5281/zenodo.4405008)

<!--
#### Example:

Human, I.M. (2020, January 1). human/myrepo: v1.2.5.hydro (Version v1.2.5.hydro). Zenodo. https://doi.org/some-doi-number
-->

[Back to Contents](#contents)

<br />

<!-------------------------->
<!-------------------------->
## Journal Reference
<!-------------------------->
<!-------------------------->
Update your journal reference here after acceptance.

[Back to Contents](#contents)

<br />

<!-------------------------->
<!-------------------------->
## Contributing Models
<!-------------------------->
<!-------------------------->

Please note that the models used in this research are the versions labeled in the table below.

| Model | Version | Repository Link | DOI |
|:-:|:-:|---|---|
| Xanthos | <v2.3.1> | <https://github.com/mengqi-z/xanthos/tree/v2.3.1-wild2020-ArgentinaNexus> | <https://doi.org/10.5281/zenodo.4404834> |
| AgMIP | <version> | <link to code repository> | <DOI link> |
| GCAM | <v5.1.3LAC> | <https://doi.org/10.5281/zenodo.3897519> | <https://doi.org/10.5281/zenodo.3897519> |
| Tethys | <v1.2.0> | <https://github.com/mengqi-z/tethys/tree/v1.2.0-wild2020-ArgentinaNexus> | <http://doi.org/10.5281/zenodo.4405008> |
| Demeter | <v1.1.0> | <https://github.com/mengqi-z/demeter/tree/v1.1.0-wild2020-ArgentinaNexus> | <https://doi.org/10.5281/zenodo.4404738> |

[Back to Contents](#contents)

<br />

<!-------------------------->
<!-------------------------->
## Data Reference
<!-------------------------->
<!-------------------------->

### 1. Input Data

#### Forcing Data
In the study of Argentina Energy-Water-Land systems, we selected one climate impact scenario from Global Climate Model (GCM) MIROC-ESM-CHEM forced by Representative Concentration Pathway (RCP) 6.0. The source of climate data for Xanthos is obtained from ISIMIP Fast Track Dataset (citation). Input data for GCAM, Demeter, and Tethys are outputs from their feeding models described in Figure 1. Generally, those outputs need to be post-processed to required formats in order to feed into other models as inputs. We provide R scripts in section [Reproduce My Experiment](#reproduce-my-experiment) for reproducing the post-processed input data.

For broader use of these data, we also provide post-processed input dataset directly in the table below. These data includes all 20 combinations of GCM/RCP scenarios. There are 5 GCMs (i.e., GFDL-ESM2M, HadGEM2-ES, IPSL-CM5A-LR, MIROC-ESM-CHEM, and NorESM1-M) and 4 RCPs (i.e., rcp2.6, rcp4.5, rcp6.0, and rcp8.5).

| Data Category | Model | DOI | Description |
|---------------|-------|-----|-------------|
| Climate | Xanthos | <link to DOI dataset> | [NPY files] 20 GCM/RCP climate projections. |
| Hydrobiologic Data | GCAM | <link to DOI dataset> | [XML files] created from Xanthos and AgYield outputs, including runoff, hydropower, and crop yields. |
| Land Allocation | Demeter | <link to DOI dataset> | [CSV files] created from land use land cover projection from GCAM output. |

#### Files Replaced for Argentina Study
All models come with default dataset and supporting files to run the associated example. It is necessary to test each model in python or java environment by running the default example. To run all the models for Argentina study, we replaced default files with modified files we provide in the table below. For configuration and model run files, you will need to modify the directories based on the location of your models. More detailed summary of data and files can be found in [here](https://).

| File Category | Model | Directory |
|---------------|-------|-----------|
| Runoff Module | Xanthos | ArgentinaNexus/Xanthos/example/input/runoff/ |
| Xanthos Configuration and Model Run | Xanthos | ArgentinaNexus/Xanthos/example/ |
| GCAM Configuration and Batch Files | GCAM | ArgentinaNexus/GCAM/ |
| Allocation | Demeter | ArgentinaNexus/Demeter/example/inputs/allocation |
| Constraints | Demeter | ArgentinaNexus/Demeter/example/inputs/constraints |
| Observation | Demeter | ArgentinaNexus/Demeter/example/inputs/observed |
| Reference | Demeter | ArgentinaNexus/Demeter/example/inputs/reference |
| Demeter Configuration | Demeter | ArgentinaNexus/Demeter/example/ |
| Region Grids | Tethys | ArgentinaNexus/Tethys/example/Input/harmonized_inputs/ |
| Region Names | Tethys | ArgentinaNexus/Tethys/example/Input/rng33/ |
| Livestock Fraction | Tethys | ArgentinaNexus/Tethys/example/Input/rng33/ |
| Electricity | Tethys | ArgentinaNexus/Tethys/example/Input/rng33/TD_Elec_paras/ |

<!--
#### Example:

Human, I.M. (2020). My dataset name [Data set]. DataHub. https://doi.org/some-doi-number
-->

### 2. Output Data
For broader use, we provide output dataset from model runs with all 20 combinations of GCM/RCP scenarios (See Table below).

| Model | DOI |
|-------|-----|
| Xanthos | <link to DOI dataset> |
| GCAM | <http://doi.org/10.5281/zenodo.4420154> | 
| Demeter | <http://doi.org/10.5281/zenodo.4420156> | 
| Tethys | <http://doi.org/10.5281/zenodo.4321776> | 

[Back to Contents](#contents)

<br />

<!-------------------------->
<!-------------------------->
## Reproduce My Experiment
<!-------------------------->
<!-------------------------->

### 1. Run Preparation

#### (A) Argentina Nexus Repository

Clone ArgentinaNexus reproducible repository into your desired location.

```
git clone https://github.com/FeralFlows/ArgentinaNexus.git
```

#### (B) GCAM LAC

***Pre-Requirements***
  
  * Downlaod GCAM v5.1.3-LAC https://doi.org/10.5281/zenodo.3897519
  * Install Java 64 http://openjdk.java.net/
  * Install Windows XML Maker http://symbolclick.com/xmlmarker_1_1_setup.exe
  
***File Replacement***

  * Replace gcam-core_LAC_v02_5Nov2019/exe/configuration_LAC.xml with ArgentinaNexus/GCAM/configuration_LAC.xml in the cloned ArgentinaNexus repository
  * Replace gcam-core_LAC_v02_5Nov2019/exe/batch_LAC.xml with ArgentinaNexus/GCAM/batch_LAC.xml
  
**Notes:** GCAM v5.1.3-LAC is a modified version from GCAM-Core-v5.1.3 for the study in Latin America and the Caribbean (LAC) Region. A 64-bit Java is required to run GCAM. We recommend the open source version of Java ([OpenJDK](http://openjdk.java.net/)). More details on GCAM installation, setting up, and trouble shooting, please refer to [GCAM Documentation](https://github.com/JGCRI/gcam-core).



#### (C) Xanthos, Demeter, and Tethys

***Pre-requirements***

  * Download Xanthos https://doi.org/10.5281/zenodo.4404834
  * Download Demeter https://doi.org/10.5281/zenodo.4404738
  * Download Tethys http://doi.org/10.5281/zenodo.4405008
  * Install PyCharm Professional version (optional) https://www.jetbrains.com/pycharm/download/#section=windows
  * Install Python 2.7 and Python 3.8 https://www.python.org/downloads/
  
***Model Installation and Setup***

The user is able to install each model as a Python package from terminal or command line. Navigate to one of the downloaded nodel folders (Xanthos, Demeter, or Tethys) and run ```python setup.py install```.

**Notes:** More details of setting up Python using PyCharm for Xanthos, Demeter, and Tethys can be found in this [PyCharm setup tutorial](https://github.com/FeralFlows/ArgentinaNexus/tree/master/docs/Python_setup_for_xanthos_demeter_tethys.md).

***File Modification***

Check each file listed in the table below and modify every directory to the directory that holds your data. For example, in configuration file 'pm_abcd_mrtm_future_impacts.ini' for xanthos model, change the directory of 'RootDir' to 'your-xanthos-location\example'.

| Model | Programming Language | Files to be Modified |
|---|:-:|---|
| Xanthos | Python 3.3+ | under xanthos/example: (1) pm_abcd_mrtm_future_impacts.ini; (2) future_sim_IDBFinalRuns.py; (3) watch_impacts.ini; (4) watch_wfdei_sim.py |
| Demeter | Python 2.7 | under demeter/example: (1) config_LAC.ini; (2) example_LAC.py |
| Tethys | Python 3 | under tethys/example: (1) config.ini; (2) example.py |




### 2. Output Processing and Model Integration

Figure 1 details the workflow for reproducing all model outputs once the input data have been prepared. 

<p align="center"> <img src="extras/paper_figs/FIG1.png", width = '700'></p>

**Figure 1.** The multi-model, multi-scale, multi-sector analysis framework.
<br /><br />

| Script | Description | Directory |
|--------|-------------|-----------|
| basin_runoff_analysis_plotting.R | Convert Xanthos runoff outputs to XML files for GCAM | ArgentinaNexus/Figures/XanthosProcessing/ |
| hydro_analysis_plotting.R | Convert Xanthos hydropower outputs to XML files for GCAM | ArgentinaNexus/Figures/XanthosProcessing/ |
| gcam_to_demeter_land_allocation_rgcam.R | Select projected landuse from GCAM output database and convert to required input format for Demeter. Put created files under /Demeter/example/inputs/projected/ | ArgentinaNexus/Figures/DemeterProcessing/ |
| aggregate_5arcmin_to_0p5degree.py | Aggregate Demeter output from 5 arcmin to 0.5 degree for further spatial landuse map plotting with metis | ArgentinaNexus/Figures/DemeterProcessing |

### 3. Reproduce Figures
We also provide scripts (ArgentinaNexus/Figures/) for reproducing figures in our paper.

| Script | Corresponding Figures | Description |
|--------|-----------------------|-----------|
| metis.masterX_Argentina.R | Figure 1 and Figure 6 | <> |
| basin_runoff_analysis_plotting.R | Figure 5 | <> |
| hydro_analysis_plotting.R | Figure 5 | <> |
| metis_plot_argentina.R | Figure 7 and Figure 8 | <> |

[Back to Top](#argentinanexus)