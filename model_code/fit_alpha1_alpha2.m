function [alpha1_alpha2] = fit_alpha1_alpha2(calibration_analytes,plotflag,prompt_output)

% Determine the number of structures in the set of calibration analytes
N_calibration = length(calibration_analytes.numdata.data(:,1));

% Determine the value of N_i_star given by eq S1 in Arey et al (2022)
N_star = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,1),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

class_number = create_class_index(calibration_analytes.class);

% These are the Abraham solute parameters of the set of calibration analytes
ESABVLc = [calibration_analytes.numdata.data(:,3:8) ones(N_calibration,1)];

% Set the LSER coefficients for column stationary phase 1, in the order [e s a b v l c]
% The default values are SE-30 coefficients taken from the 1999 review by Abraham
LSERcoeffs_1 = [0.024 0.190 0.125 0.0 0.0 0.498 -0.194];

% Determine the gas-stationary phase partition coefficients at 121 C according to the LSER model
% for GCxGC column 1
logL1 = sum(ESABVLc.*(ones(N_calibration,1)*LSERcoeffs_1),2);

% Set the number of synthetic replicates to be applied in the bootstrap procedure below.
N_boot = 10000;

% Solve the values of alpha1 and alpha2 by linear regression of eq 3 in Arey et al (2022).
% Estimate the uncertainties of alpha1 and alpha2, reported as the 95% interval of bootstrapped parameter values.
[alpha1_alpha2, alpha1_alpha2_unc, logL1_pred, r2_logL1, r2_logL1_unc] = svd_regress_boot([N_star ones(N_calibration,1)],logL1,N_boot);

if ( 0.5*mean(alpha1_alpha2_unc(:,2)) > abs(alpha1_alpha2(2)) )
 disp('During the regression of eq 6, the intercept is found to be statistically equivalent to zero.');
 disp('Refitting with alpha2 set to zero.');
 [alpha1_alpha2, alpha1_alpha2_unc, logL1_pred, r2_logL1, r2_logL1_unc] = svd_regress_boot(N_star,logL1,N_boot);
 alpha1_alpha2 = [alpha1_alpha2 0];
 alpha1_alpha2_unc = [alpha1_alpha2_unc zeros(2,1)];
end

alpha1_alpha2_unc = mean(alpha1_alpha2_unc);
rmse_L1 = sqrt(sum((logL1_pred-logL1).^2)./length(logL1));

if strcmp(prompt_output,'verbose')
 disp('Fitted alpha_1 and alpha_2 values by eq 3 are:');
 disp(alpha1_alpha2);

 disp('Bootstrap uncertainty estimates of alpha_1 and alpha_2 are:');
 disp(alpha1_alpha2_unc);

 disp('The resulting r^2 and RMSE values of model-fitted logL1 values by eq 3 are:'); 
 disp([r2_logL1 rmse_L1])
end

% Plot the model-fitted logL1 values for the set of calibration analytes
if plotflag > 0
 colorlist = load('../model_parameters/colorlist_class.dat');
 figure;
 ax = axes;
 % Plot classes according to the following color order, except nP class which is plotted separately below
 for class_ind = unique(class_number')
  colorindex = colorlist(class_ind);
  hold on;
  ax.ColorOrderIndex = colorindex;
  plot(logL1(find(class_number==class_ind)),logL1_pred(find(class_number==class_ind)),'o');
 end 
 % Exceptionally, plot nP class with black circles
 plot(logL1(find(class_number==11)),logL1_pred(find(class_number==11)),'ko');
 plot(min(logL1):0.1:max(logL1),min(logL1):0.1:max(logL1),'k--');
 hold off;
 xlabel('LSER-Estimated log L_1 Values');
 ylabel('Elution Model Fitted log L_1 Values');
 box on;
end
