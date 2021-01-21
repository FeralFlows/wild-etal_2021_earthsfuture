# Summary of Modified Model Files for Argentina Nexu Study

<br />

## Contents
- [Xanthos](#xanthos)
- [GCAM LAC](#gcam-lac)
- [Demeter](#demeter)
- [Tethys](#tethys)

<br />

## Xanthos

'~/' in Table 1 represents your-xanthos-location/.

**Table 1:** File modification for Xanthos model.

| Data and File Category | File Name | Directory | Changes | Replace Original Files | Add-on | Notes |
|-|-|-|-|-|-|-|
| Climate data (5GCMs x 4RCPs) | gfdl-esm2m<br>hadgem2-es<br>ipsl-cm5a-lr<br>miroc-esm-chem<br>noresm1-m<br>watch+wfdei | ~/example/input/climate/ | climate input | Y |  | pr, rhs, rlds, rsfs, tas,tasmin,   wind |
| Runoff Module files | pars_watch_1971_1990_decadal_lc.npy | ~/example/input/runoff/abcd | model calibration | Y |  |  |
| Hydropower Files | simulated_cap_by_country<br>gridData.csv | ~/example/input/hydropower_actual | include Uruguay | Y |  |  |
| Reference Files | Rgn32Names.csv<br>region32_grids.csv | ~/example/input/input/reference | include Uruguay | Y |  |  |
| Configuration File | pm_abcd_mrtm_future_impacts.ini<br>watch_impacts.ini | ~/example/ | configuration | Y |  |  |
| Model Run file | future_sim_IDBFinalRuns.py<br>watch_wfdei_sim.py | ~/example/ | for climate change run<br>for historical run | Y |  | run with 5 GCMs data<br>     run with the watch+wfdei data |
| Post Processing R Script | xanthos_postprocessing_fns.R<br>basin_runoff_analysis_plotting_MZ.R<br>hydro_analysis_plotting_MZ.R | Your choice | post-processing functions<br>post-processing runoff output<br>post-processing hydropower output |  | Y |  |
| Files Used in R Scripts | gcam_basin_id.csv<br>gcam_country_id.csv<br>GCAM_region_names.csv<br>gcam_xanthos_basin_mapping.csv<br>GCAMBasin_country.csv<br>iso_GCAM_regID.csv<br>Rgn32Names.csv<br>runoff_max_wfdei_1970_2010.csv<br>basin_ID.csv<br>L201.RenewRsrcCurves_calib_watergap.csv<br>L201.GrdRenewRsrcMax_runoff.csv<br>L103.water_mapping_R_GLU_B_W_Ws_share.csv<br>L103.water_mapping_R_B_W_Ws_share.csv<br>reference_gcam_hydro_all_regions.csv<br>basin_to_country_mapping.csv | Your choice | Reference files for post processing |  | Y |  |

<br />

## GCAM LAC
'~/' in Table 2 represents your-gcam-lac-location/input/.

**Table 2:** File modification for GCAM LAC.

| Config and Batch File | Scenario | Region | Sector | XML | Description |
|-|-|-|-|-|-|
| Configuration_LAC.xml | Reference | Global |  | LACglobal_ref_gas_cc_sw_interp.xml | ~/idb/reference/ |
| Configuration_LAC.xml | Reference | Global |  | LACglobal_ref_electricity_water.xml | ~/idb/reference/ |
| Configuration_LAC.xml | Reference | Global |  | LACglobal_ref_transportation_UCD_CORE.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Socioeconomics | Colombia_refPopGDP_gSSP2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Socioeconomics | Colombia_refPopGDP_demand_input.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Socioeconomics | Colombia_refPopGDP_cement_incelas_gssp2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Socioeconomics | Colombia_refPopGDP_industry_incelas_gssp2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Energy | Colombia_ref_hydro.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Water | Colombia_ref_livestock_water_demand_coeff.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Water | Colombia_ref_L210.TechCoef_mod.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Water | Colombia_ref_L232.TechCoef_mod.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Colombia | Water | Colombia_ref_L2072.AgCoef_IrrWaterWdraw_ag_mgmt_mod2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_socioeconomics_gSSP2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_industry_incelas_gssp2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_cement_incelas_gssp2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_bld_agg_gSSP2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_trn_agg_gSSP2.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Socioeconomics | Uruguay_refPopGDP_demand_input.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Energy | Uruguay_refFinalNrg.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Electricity | Uruguay_refElecShareWeights_L223.SubsectorShrwt.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Electricity | Uruguay_refElecHydro.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Water | Uruguay_refWater.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Water | Uruguay_refWat_L2072.AgCoef_IrrWaterCons_ag_mgmt.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Uruguay | Water | Uruguay_refWat_L2072.AgCoef_IrrWaterWdraw_ag_mgmt.xml | ~/idb/reference/ |
| batch_LAC.xml | Reference | Argentina | Socioeconomics | Argentina_refPop_gSSP2.xml | ~/idb/reference/ |
| batch_LAC.xml | Impacts | Global | Energy | hydro_impacts_HadGEM2-ES_rcp8p5.xml | ~/idb/impacts/Hydro/ |
| batch_LAC.xml | Impacts | Global | Water | runoff_impacts_HadGEM2-ES_rcp8p5.xml | ~/idb/impacts/Water/ |
| batch_LAC.xml | Impacts | Global | Agriculture | ag_prodchange_HadGEM2-ES_rcp8p5.xml | ~/idb/impacts/Ag/ |
| batch_LAC.xml | Policy | Global | Policy | CarbonTax_CO2_NonCO2_CO_ARG.xml | ~/idb/policy/ |
| batch_LAC.xml | Policy | Global | Policy | CarbonTax_CO2_NonCO2_LINK_CO_ARG.xml | ~/idb/policy/ |

<br />

## Demeter

'~/' in Table 3 represents your-demeter-location/example/.

**Table 3:** File modification for Demeter model.

| File or Folder Category | File Name | Directory | Changes | Replace Original Files | Add-on |
|-|-|-|-|-|-|
| Allocation directory | gcam_regbasin_modis_v6_type5_mirca_5arcmin_constraint_alloc.csv<br>gcam_regbasin_modis_v6_type5_mirca_5arcmin_kernel_weighting.csv<br>gcam_regbasin_modis_v6_type5_mirca_5arcmin_observed_alloc.csv<br>gcam_regbasin_modis_v6_type5_mirca_5arcmin_order_alloc.csv<br>gcam_regbasin_modis_v6_type5_mirca_5arcmin_projected_alloc.csv<br>gcam_regbasin_modis_v6_type5_mirca_5arcmin_transition_alloc.csv | ~/inputs/allocation | 5arcmin, 35 crops | Y |  |
| Constraint data directory | 000_nutrientavail_hswd_5arcmin.csv<br>001_soilquality_hswd_5arcmin.csv | ~/inputs/constraints | 5arcmin | Same |  |
| Observed spatial data directory | gcam_reg32_basin235_modis_v6_2010_mirca_2000_5arcmin_sqdeg_wgs84_11Jul2019 | ~/inputs/observed | 5arcmin | Y |  |
| Projected GCAM land allocation directory | DemeterDownscaled_33Regions_MIROC-ESM-CHEM_rcp6p0_Reference.csv<br>DemeterDownscaled_33Regions_MIROC-ESM-CHEM_rcp6p0_Impacts.csv<br>DemeterDownscaled_33Regions_MIROC-ESM-CHEM_rcp6p0_Policy.csv | ~/inputs/projected | 33 regions, GCAM GLU (235   Basins), crops | Y |  |
| Reference directory | gcam_regions_33.csv | ~/inputs/reference | 33 gcam regions (Uruguay   included) | Y |  |
| Configuration | config_LAC.ini<br>example_LAC.py | ~/ | updated configuration file<br>     updated run file | Y |  |
| GCAM output to Demeter | query_demeter_33regions_3scenarios.xml<br>gcam_to_demeter_land_allocation_rgcam.R | Figures/DemeterProcessing |  |  | Y |
| Post Processing | aggregate_5arcmin_to_0p5degree.py<br>gcam_regbasin_5arcmin_to_0p5deg.txt<br>coords_0p5_deg_wgs84.csv | Figures/DemeterProcessing |  |  |  |

<br />

## Tethys

'~/' in Table 4 represents your-tethys-location/example/.

**Table 4:** File modification for Tethys model.

| File or Folder Category | File Name | Directory | Changes | Replace Original Files |
|-|-|-|-|-|
| gcam5p1_ref_db folder | atv.basex<br>inf.basex<br>tbl.basex<br>tbli.basex<br>txt.basex<br>txtl.basex<br>txtr.basex | ~/Input/GCAM | GCAM output | Y |
| harmonized_inputs folder | region33_grid.csv | ~/Input/harmonized_inputs | Added Uruguay as individual region | Y |
| rgn33 folder | bfracFAO2005.csv<br>gfracFAO2005.csv<br>RgnNames.csv<br>ElecBuilding_1971_2010.csv<br>ElecBuldingCool_1971_2010.csv<br>ElecBuildingHeat_1971_2010.csv<br>ElecBuildingOthers_1971_2010.csv<br>ElecIndustry_1971_2010.csv | ~/Input/rgn33<br> ~/Input/rgn33<br> ~/Input/rgn33<br> ~/Input/rgn33/TD_Elec_paras<br> ~/Input/rgn33/TD_Elec_paras<br> ~/Input/rgn33/TD_Elec_paras<br> ~/Input/rgn33/TD_Elec_paras<br> ~/Input/rgn33/TD_Elec_paras | Added Uruguay as individual region | Y |
| configuration | config.ini | ~/ | Directories and parameters | Y |

