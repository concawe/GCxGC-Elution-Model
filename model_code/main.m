% The file main.m represents the principal code of the GCxGC elution model. It calls several functions
% which implement the algorithms of the model and then plot the results.

% ------------------------------------------------------------------------------------------------------------

% BEGINNING OF SECTION: LOAD INPUT DATA AND MODEL PARAMETERS

% Load measured data of the retention times (floating point) and carbon number values (integers) of the n-alkanes.
% The retention time data of the n-alkanes are used to conduct the interpolations represented by equations S1 and S2 
% in Arey et al (2022). 
calibration_analytes.alkanes = load(strcat('../',input_path,alkanes_data_file));

% Load measured GCxGC retention time data and UFZ-LSER parameters (floating point) of all calibration analytes 
% (including n-alkanes). These data are used to fit the tunable parameters (alpha1, alpha2, alpha3) in model 
% equations 3 and 4.
calibration_analytes.numdata = importdata(strcat('../',input_path,calibration_data_file));

% Load the slope and intercept parameters of eq 5 in Arey et al (2022).
calibration_analytes.logL2_alkane_Nstar_parms = load('../model_parameters/logL2_alkane_Nstar_parms.dat');

% Load the assigned class values (string) for calibration analytes
fid1 = fopen(strcat('../',input_path,calibration_classes_file));
calibration_analytes.text = textscan(fid1,'%s');
fclose(fid1);
calibration_analytes.class = string(calibration_analytes.text{1});

% Load the UFZ-LSER-estimated values of E, S, A, B, V, and L (floating point) for the Constituent Library.
library.ESABVL = importdata(strcat('../model_parameters/',library_ESABVL_file),'\t');

% Load the assigned values of class (text) and carbon number (integers) for the Constituent Library.
fid2 = fopen(strcat('../model_parameters/',library_class_CN_file));
library.text = textscan(fid2,'%u%s');
fclose(fid2);

library.CN = library.text{1};
library.class = string(library.text{2});
library.HCBs = unique([library.class library.CN],'rows');

N_calibration = length(calibration_analytes.numdata(:,1));
N_library = length(library.ESABVL.data(:,1));

% END OF SECTION: LOAD INPUT DATA AND MODEL PARAMETERS

% ------------------------------------------------------------------------------------------------------------

% BEGINNING OF SECTION: PRINCIPAL ALGORITHMS OF THE GCXGC ELUTION MODEL 

% Predict values of logL1 and Delta_logL12_library for constituent library of 15495 structures, and determine
% the beta value of this library by Schmidt orthogonalization if needed
[logL1_DeltalogL21_library, beta] = determine_logL1_Delta_logL21(library.ESABVL.data, beta_input, output_path, prompt_output);

% Fit the tunable model parameters alpha1 and alpha2 by eq 3 in Arey et al (2022) with the set of calibration analytes.
alpha1_alpha2 = fit_alpha1_alpha2(calibration_analytes, generic_plot_flag, prompt_output);

% Fit the tunable model parameter alpha3 by eq 4 in Arey et al (2022) with the set of calibration analytes.
alpha3 = fit_alpha3(calibration_analytes, alpha1_alpha2, beta, generic_plot_flag, output_path, prompt_output);

% The array logL1_DeltalogL21_fitted represents model-fitted values of logL1 and Delta_logL21 that were determined from 
% regression of eqs 3 and 4 to measured retention time data (t1, t2) during the calibration of alpha1_alpha2 and alpha3.
logL1_DeltalogL21_fitted = load(strcat('../',output_path,'logL1_DeltalogL21_fitted.dat'));

% The array logL1_DeltalogL21_fitted_test represents the model-predicted values of logL1 and Delta_logL21 that are determined 
% from a forward calculation of eqs 3 and 4, using the measured retention time data (t1, t2) as input together with the 
% previously calibrated values of alpha1_alpha2 and alpha3. This is a simple test of the forward model.
logL1_DeltalogL21_fitted_test = test_fitted_logL_values(calibration_analytes, alpha1_alpha2, alpha3, beta, output_path);

% We should find that logL1_DeltalogL21_fitted is equivalent to logL1_DeltalogL21_fitted_test.
% This test of equivalency employs rounded data to avoid spurious failures arising from floating point precision.
if round(logL1_DeltalogL21_fitted_test,3) ~= round(logL1_DeltalogL21_fitted,3)
 disp('WARNING. Failed to produced consistent logL1 and logL2 values.');
end

% The file logL1_DeltalogL21_calibration.dat was created by the function fit_alpha3, and these data represent the logL1 and 
% DeltalogL21 values that are determined by the combination of UFZ LSER-estimated ESABVL parameters (Ulrich et al 2017) with 
% published gas-stationary phase Abraham coefficients for the appropriate phases.
logL1_DeltalogL21_calibration = load(strcat('../',output_path,'logL1_DeltalogL21_calibration.dat'));

% We determine model-fitted values of retention time 1 and retention time 2 for the set of calibration analytes, t1_calib and 
% t2_calib, by applying LSER-determined values of logL1 and DeltalogL21 to the elution model equations.
[t1_calib t2_calib] = convert_logL1_DeltalogL21_to_t1_t2(logL1_DeltalogL21_calibration, alpha1_alpha2, alpha3, beta, calibration_analytes);

% Save the model-fitted t1_calib and t2_calib values of the calibration analytes to file
t1t2_calib = [t1_calib t2_calib];
t1t2_calib_savefilename = strcat('../',output_path,'t1_t2_fitted.dat');
save(t1t2_calib_savefilename,'t1t2_calib','-ascii');

% We can determine model-fitted values of retention time 1 and retention time 2 for the set of calibration analytes
% by applying model-fitted values of logL1 and DeltalogL21 to the elution model equations. This quantity is not
% reported, but a comparison of t1_fit_test and t2_fit_test with the conventional model-fitted values (t1_calib and t2_calib) 
% can provide an estimate of the uncertainty in modeled retention times arising from the model fit itself.
[t1_fit_test t2_fit_test] = convert_logL1_DeltalogL21_to_t1_t2(logL1_DeltalogL21_fitted, alpha1_alpha2, alpha3, beta, calibration_analytes);

% Predict retention time 1 and retention time 2 for the constituent library, by applying the tuned calibration parameters 
% (alpha1_alpha2, alpha3) together with the measured retention time data of the n-alkanes.
[t1_lib t2_lib] = convert_logL1_DeltalogL21_to_t1_t2(logL1_DeltalogL21_library, alpha1_alpha2, alpha3, beta, calibration_analytes);

% Save the simulated t1_lib and t2_lib values of the constituent library to file
t1t2_lib = [t1_lib t2_lib];
t1t2_lib_savefilename = strcat('../',output_path,'t1_t2_library.dat');
save(t1t2_lib_savefilename,'t1t2_lib','-ascii');

% END OF SECTION: PRINCIPAL ALGORITHMS OF THE GCXGC ELUTION MODEL 

% ------------------------------------------------------------------------------------------------------------

% BEGINNING OF SECTION: PLOT RESULTS OF THE GCXGC ELUTION MODEL 

if generic_plot_flag

% Plot measured retention time data (triangles) versus model-fitted retention time data (circles) for the set of 
% calibration analytes, color-coded by class assignment.
 plot_fitted_retention_times_calibration(calibration_analytes, t1_calib, t2_calib, modulation_period, prompt_output);

% Plot the predicted retention time data of the constituent library (black dots) overlaid with positions of the 
% calibration analytes (red circles). 
 plot_predicted_retention_times_library(calibration_analytes, t1_calib, t2_calib, t1_lib, t2_lib, modulation_period, prompt_output);

% Plot the predicted retention time data of the constituent library, color-coded according to class assignment. 
 plot_predicted_retention_times_library_carbon_number_coded(calibration_analytes, library, t1_lib, t2_lib, modulation_period, prompt_output);

% Plot the predicted retention time data of the constituent library, color-coded according to class assignment. 
 plot_predicted_retention_times_library_class_coded(calibration_analytes, library, t1_lib, t2_lib, modulation_period, prompt_output);

% Plot the predicted retention time data of the constituent library, color-coded according to class assignment with maximum-area polygons
 plot_predicted_retention_times_library_class_coded_polygons(calibration_analytes, library, t1_lib, t2_lib, modulation_period, prompt_output);

end

if single_CN_plot_flag

% Plot the predicted retention time data of members of the constituent library sharing a single carbon number, color-coded according to class assignment
 plot_predicted_retention_times_library_class_single_CN(calibration_analytes, library, t1_lib, t2_lib, modulation_period, single_CN_value, prompt_output);

end

