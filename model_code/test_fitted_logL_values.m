function logL1_DeltalogL21_fitted_test = test_fitted_logL_values(calibration_analytes,alpha1_alpha2,alpha3,beta,output_path);

% Determine the number of structures in the set of calibration analytes
N_calibration = length(calibration_analytes.numdata.data(:,1));

% These are the Abraham solute parameters of the set of calibration analytes
ESABVLc = [calibration_analytes.numdata.data(:,3:8) ones(N_calibration,1)];

%{
% Set the LSER coefficients for column stationary phase 1, in the order [e s a b v l c]
% The default values are SE-30 coefficients taken from the 1999 review by Abraham
LSERcoeffs_1 = [0.024 0.190 0.125 0.0 0.0 0.498 -0.194];

% Set the LSER coefficients for column stationary phase 2, in the order [e s a b v l c]
% The default values are OV-17 coefficients taken from the 1999 review by Abraham 
LSERcoeffs_2 = [0.071 0.653 0.263 0.0 0.0 0.518 -0.372];

logL1_cal = sum(ESABVLc.*(ones(N_calibration,1)*LSERcoeffs_1),2);
logL2_cal = sum(ESABVLc.*(ones(N_calibration,1)*LSERcoeffs_2),2);
%}

% Determine the value of N_i_star given by eq S1 in Arey et al (2022)
N_star_test = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,1),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

% This is the fitted set of logL1 values of the set of calibration analytes, according to eq 3 in Arey et al (2022)
% The tunable parameters, alpha1 and alpha2, were determined in a previous function
logL1_test = alpha1_alpha2(1)*N_star_test + alpha1_alpha2(2);

% Determine the value of t_2,i_star given by eq S2 in Arey et al (2022)
t2_star_test = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,3),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

logL2ref_To = calibration_analytes.logL2_alkane_Nstar_parms(1)*N_star_test + calibration_analytes.logL2_alkane_Nstar_parms(2);

logL2_test = logL2ref_To + log10((calibration_analytes.numdata.data(:,2) - alpha3)./(t2_star_test - alpha3));

%{
u_1_cal = logL1_cal;
u_2_cal = logL2_cal-beta*logL1_cal;
u12_cal = [u_1_cal u_2_cal];
savefilename = strcat('../',output_path,'u1_u2_values_red_calib_set.dat');
save(savefilename,'u12_cal','-ascii');
%}

logL1_DeltalogL21_fitted_test = [logL1_test logL2_test-beta*logL1_test];
savefilename = strcat('../',output_path,'logL1_DeltalogL21_fitted_test.dat');
save(savefilename,'logL1_DeltalogL21_fitted_test','-ascii');
