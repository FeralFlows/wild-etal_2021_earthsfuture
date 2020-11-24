import os
import pandas as pd


def agg_5arcmin_to_0p5degree(demeter_data_dir, f_5arcmin, f_0p5_coords, yr, pfts, n_grids):
    """Aggregate Demeter fractional output from 5 arcmin to 0.5 degree cells"""
    
    f_dem = os.path.join(demeter_data_dir, 'landcover_{}_timestep.csv')

    # target_fid of each 5 arcmin grid cell to nearest OBJECTID of 0p5 degree cell
    near_df = pd.read_csv(f_5arcmin, usecols=['target_fid', 'NEAR_FID', 'NEAR_DIST'])
    near_df.rename(columns={'NEAR_FID': 'fid_0p5_deg', 'NEAR_DIST': 'dist_0p5_deg'}, inplace=True)

    # coordinates to 0p5 degree cells
    coord_df = pd.read_csv(f_0p5_coords, usecols=['pkey', 'OBJECTID', 'latitude', 'longitude'])
    # change latitude and longitude name to fit with the metis grid2poly
    coord_df.rename(columns={'pkey': 'pkey_0p5_deg', 'OBJECTID': 'fid_0p5_deg'}, inplace=True)

    near_df = pd.merge(near_df, coord_df, on='fid_0p5_deg', how='left')
    near_df.reset_index(inplace=True)

    # demeter output data in units -> fraction
    dem_df = pd.read_csv(f_dem.format(yr))
    dem_df['target_fid'] = dem_df.index + 1

    mdf = pd.merge(dem_df, near_df[['target_fid', 'pkey_0p5_deg', 'dist_0p5_deg']].copy(), on='target_fid', how='left')
    mdf = mdf.loc[mdf['dist_0p5_deg'] == 0]
    mdf['tally'] = 1
    mdf = mdf.groupby('pkey_0p5_deg').sum()
    mdf.reset_index(inplace=True)
    mdf.drop(['latitude', 'longitude', 'basin_id', 'region_id', 'target_fid'], axis=1, inplace=True)

    # adjust fraction from 5 arcmin to 0.5 degree grid
    mdf[pfts] /= n_grids

    xdf = pd.merge(mdf, coord_df, on='pkey_0p5_deg', how='left')
    xdf.drop(['fid_0p5_deg', 'dist_0p5_deg', 'tally'], axis=1, inplace=True)

    xdf.to_csv(os.path.join(demeter_data_dir, 'demeter_{}_0p5deg_{}.csv'.format(scenario, yr)), index=False)


if __name__ == '__main__':
    
    root_dir = 'E:/Nexo-UA/Demeter/example/outputs'
    agg_dir = 'E:/Nexo-UA/Results/downscaling/land/demeter/PostProcess'
    
    f_5arcmin = os.path.join(agg_dir, 'gcam_regbasin_5arcmin_to_0p5deg.txt')
    f_0p5_coords = os.path.join(agg_dir, 'coords_0p5_deg_wgs84.csv')

    # the number of 5 arcmin grid cells in each 0.5 degree cell
    n_grids = 36

    runs = ['reference_2020-10-28_20h29m05s',
            'impacts_2020-10-29_08h31m08s',
            'policy_2020-10-29_09h21m42s']
    
    pfts = ['water', 'forest',
            'shrub', 'grass', 'urban', 'snow', 'sparse', 'corn_irr',
            'fibercrop_irr', 'foddergrass_irr', 'fodderherb_irr', 'misccrop_irr',
            'oilcrop_irr', 'othergrain_irr', 'palmfruit_irr', 'rice_irr',
            'root_tuber_irr', 'sugarcrop_irr', 'wheat_irr', 'corn_rfd',
            'fibercrop_rfd', 'foddergrass_rfd', 'fodderherb_rfd', 'misccrop_rfd',
            'oilcrop_rfd', 'othergrain_rfd', 'palmfruit_rfd', 'rice_rfd',
            'root_tuber_rfd', 'sugarcrop_rfd', 'wheat_rfd', 'otherarableland',
            'biomass_grass_irr', 'biomass_grass_rfd', 'biomass_tree_irr', 'biomass_tree_rfd']
    
    gcam_years = range(2005, 2055, 5)

    
    for run in runs:
    
        demeter_data_dir = os.path.join(root_dir, run, 'spatial_landcover_tabular')

        scenario = run.split('_')[0]

        for yr in gcam_years:
            
            agg_5arcmin_to_0p5degree(demeter_data_dir, f_5arcmin, f_0p5_coords, yr, pfts, n_grids)
