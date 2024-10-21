
edit CP_general_parameters.m

%% preprocessing
% in-situ observations
edit CP_prep_saildrone.m
edit CP_prep_PIRATA_buoys.m
edit CP_prep_NDBC_buoys.m

% IBTrACS storm tracks
edit CP_prep_IBTrACS_raw.m
edit CP_prep_IBTrACS_final.m

% SHIPS storm tracks
edit CP_prep_SHIPS_raw.py % use python to run this code
edit CP_prep_SHIPS_final.m

% heat flux
edit CP_prep_HF.m
   
%% cold pool definition
edit CP_definition_grad_vs_diff_case_study.m
edit CP_slope_computation.m
edit CP_definition_Ta_wind_thresholds.m

%% possibility of using buoys for CP study
edit CP_detection_saildrone.m
edit CP_detection_grad_vs_diff10.m
edit CP_detection_grad_vs_diff2.m
edit CP_plot_CP_characteristics_saildrone.m


%% compute CPs from buoys
edit CP_PIRATA_coldpools_detection.m
edit CP_NDBC_coldpools_detection.m

%% CPs vs. TCs
edit CP_plot_datacoverage.m % plot TC tracks and buoys locations

% cold pool maps in TC
edit CP_plot_CPs_in_TC_coordinates_example.m
edit CP_prep_CPs_and_active_gauges_IBTRACS.m
edit CP_prep_CPs_and_active_gauges_in_VWS.m


edit CP_plot_CP_map_in_IBTrACS_motion.m % quadrant maps and histogram
edit CP_plot_CP_map_in_SHIPS_VWS.m
edit CP_plot_CP_map_in_SHIPS_motion.m
edit CP_plot_CP_map_in_VWS_and_motion.m
edit CP_plot_CP_map_land_effects.m

% probability maps
edit CP_plot_CP_probabililty_maps_in_TCs.m

% CP periods & CP counts in different quadrants
edit CP_count_and_period_wrt_TC_quadrant.m

% heatflux
edit CP_plot_HF_map_in_TCs.m
edit CP_plot_HF_map_in_TCs_SHIPS.m

%% functions
edit CP_HF_coare35vn.m % bulk formula V3.5
edit CP_function_histogram_bootstrap.m % CP time
edit CP_function_histogram_bootstrap_v2.m % scaled CP time


    
    
    
    
    
    
    
    
    
    
    
    
    
        
