# ArgentinaNexus
<!--your zenodo badge here-->

Quantifying the value of emulation for estimating drought statistics

## Contents
- [Abstract](#abstract)
- [Code Reference](#code-reference)
- [Journal Reference](#journal-reference)
- [Contributing Models](#contributing-models)
- [Data Reference](#data-reference)
- [Reproduce My Experiement](#reproduce-my-experiement)

<!----------------------------------------------------------------------------->
## Abstract
Your abstract here.

[Back to Contents](#contents)

<!----------------------------------------------------------------------------->
## Code Reference
References for each minted software release for all code involved.  If you have modified a codebase that is outside of a formal release, and the modifications are not planned on being merged back into a version, fork the parent repository and add a `.<shortname>` to the version number of the parent and conduct your own name.  For example, `v1.2.5.hydro`.

#### Example:

Human, I.M. (2020, January 1). human/myrepo: v1.2.5.hydro (Version v1.2.5.hydro). Zenodo. https://doi.org/some-doi-number

[Back to Contents](#contents)

<!----------------------------------------------------------------------------->
## Journal Reference
Update your journal reference here after acceptance.

[Back to Contents](#contents)

<!----------------------------------------------------------------------------->
## Contributing Models

<p align="center"> <img src="extras/paper_figs/FIG2.png"></p>

| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| Xanthos | <v2.3.1> | <https://github.com/JGCRI/xanthos> | <link to DOI dataset> |
| AgMIP | <version> | <link to code repository> | <link to DOI dataset> |
| GCAM | <v5.1> | <https://zenodo.org/record/3897519#.X20P-mhKiUk> | <link to DOI dataset> |
| Tethys | <v1.2.0> | <https://github.com/JGCRI/tethys> | <link to DOI dataset> |
| Demeter | <v1.1.0> | <https://github.com/JGCRI/demeter> | <link to DOI dataset> |

[Back to Contents](#contents)

<!----------------------------------------------------------------------------->
## Data Reference

### 1. Input Data
Reference for each minted data source for your input data.

| Data Category | Model | DOI |
|---------------|-------|-----|
| Climate | Xanthos | <link to DOI dataset> |

#### Example:

Human, I.M. (2020). My dataset name [Data set]. DataHub. https://doi.org/some-doi-number

### 2. Output Data
Reference for each minted data source for your output data.

[Back to Contents](#contents)


<!----------------------------------------------------------------------------->
## Reproduce My Experiement
Fill in detailed info here or link to other documentation that is a thorough walkthrough of how to use what is in this repository to reproduce your experiment.

### 1. Replace Files for Argentina Study
All models come with default dataset and supporting files to run the associated example. It is necessary to test each model in python or java environment by running the default example. To run all the models for Argentina study, please replace default files with modified files we provide in the table below. For configuration and model run files, you will need to modify the directories based on the location of your models. More detailed summary of data and files can be found in [here](https://).

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

### 2. Output Processing and Model Integration

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

[Back to Contents](#contents)