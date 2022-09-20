function [alpha3] = fit_alpha3(calibration_analytes,alpha1_alpha2,beta,plotflag,output_path,prompt_output)

% Determine the value of N_i_star given by eq S1 in Arey et al (2022)
N_star = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,1),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

% This is the fitted set of logL1 values of the set of calibration analytes, according to eq 3 in Arey et al (2022)
% The tunable parameters, alpha1 and alpha2, were determined in a previous function
logL1_pred = polyval(alpha1_alpha2,N_star);

% Determine the number of structures in the set of calibration analytes
N_calibration = length(calibration_analytes.numdata.data(:,1));

class_number = create_class_index(calibration_analytes.class);

% These are the Abraham solute parameters of the set of calibration analytes
ESABVLc = [calibration_analytes.numdata.data(:,3:8) ones(N_calibration,1)];

% Set the LSER coefficients for column stationary phase 1, in the order [e s a b v l c]
% The default values are SE-30 coefficients taken from the 1999 review by Abraham
LSERcoeffs_1 = [0.024 0.190 0.125 0.0 0.0 0.498 -0.194];

% Set the LSER coefficients for column stationary phase 2, in the order [e s a b v l c]
% The default values are OV-17 coefficients taken from the 1999 review by Abraham 
LSERcoeffs_2 = [0.071 0.653 0.263 0.0 0.0 0.518 -0.372];

% Determine the gas-stationary phase partition coefficients at 121 C according to the LSER model
% for both column 1 and column 2
logL1 = sum(ESABVLc.*(ones(N_calibration,1)*LSERcoeffs_1),2);
logL2 = sum(ESABVLc.*(ones(N_calibration,1)*LSERcoeffs_2),2);

% Determine Delta_logL21 for the calibration analytes according to eq 6 in Arey et al (2022)
Delta_logL21 = logL2-beta*logL1;

% The trial guess for alpha3 is set to 1/4 of the smallest value among the available measured 
% 2nd dimension retention time values of the set of calibration analytes
alpha3_guess = min([calibration_analytes.numdata.data(:,2); calibration_analytes.alkanes(:,2)])/4;

if strcmp(prompt_output,'verbose')
 disp('Now fitting alpha_3 with a nonlinear optimization of eq 4.');
end

% Solve the value of alpha3 by implementing a nonlinear fit of eq 4 in Arey et al (2022).
% To improve the robustness of the parameter fit, the algorithm minimizes residuals of the 
% fitted value of Delta_logL21 (eq 6) which incorporates logL2 (eq 4).
alpha3 = lsqcurvefit(@(x,calibration_analytes) convert_t2_to_logL2(x,calibration_analytes)-beta*logL1_pred, alpha3_guess, calibration_analytes, Delta_logL21);

logL2_pred = convert_t2_to_logL2(alpha3,calibration_analytes);

% Report the model-fitted values of Delta_logL21 for the calibration analytes according to eq 6 in Arey et al (2022)
Delta_logL21_pred = logL2_pred - beta*logL1_pred;

% Estimate the uncertainty of alpha3, reported as the 95% interval of bootstrapped parameter values.
N_boot = 10000;

if strcmp(prompt_output,'verbose')
 disp('Conducting a bootstrap uncertainty analysis of alpha_3. This may take a minute.');
end

rng('shuffle');
for ind_b = 1:N_boot
 ind_r = randi(N_calibration,1,N_calibration);
 calibration_analytes_b = calibration_analytes;
 calibration_analytes_b.numdata.data(:,1:2) = calibration_analytes.numdata.data(ind_r,1:2);
 N_star_b = N_star(ind_r);
 logL1_b = logL1(ind_r);
 p1_b = polyfit(N_star_b,logL1_b,1);
 logL1_pred_b = polyval(p1_b,N_star_b);
 logL2_b = logL2(ind_r);
 Delta_logL21_b = logL2_b-beta*logL1_b;
 options = optimset('display','off');
 alpha3_b(ind_b) = lsqcurvefit(@(x,calibration_analytes_b) convert_t2_to_logL2(x,calibration_analytes_b)-beta*logL1_pred_b, alpha3_guess, calibration_analytes_b, Delta_logL21_b, [], [], options);
end

alpha3_median_boot = quantile(alpha3_b,0.5);
alpha3_unc(1) = alpha3-quantile(alpha3_b,0.025);
alpha3_unc(2) = quantile(alpha3_b,0.975)-alpha3;
alpha3_unc = mean(alpha3_unc);

if strcmp(prompt_output,'verbose')
 disp('The fitted alpha_3 value is:');
 disp(alpha3);

 disp('The bootstrap uncertainty estimate of alpha_3 is:');
 disp(alpha3_unc);
end

% Plot the model-fitted Delta_logL21 values for the set of calibration analytes
if plotflag > 0
 colorlist = load('../model_parameters/colorlist_class.dat');
 figure;
 ax = axes;
 % Plot classes according to the following color order, except nP class which is plotted separately below
 colorlist = [6 5 7 1 4 6 4 3 5 2 1];
  for class_ind = unique(class_number')
  colorindex = colorlist(class_ind);
  hold on;
  ax.ColorOrderIndex = colorindex;
  plot(Delta_logL21(find(class_number==class_ind)), Delta_logL21_pred(find(class_number==class_ind)),'o');
  end
 plot(Delta_logL21(find(class_number==11)), Delta_logL21_pred(find(class_number==11)),'ko');
 plot(min(Delta_logL21):0.1:max(Delta_logL21),min(Delta_logL21):0.1:max(Delta_logL21),'k--');
 hold off;
 xlabel('LSER-Estimated \Deltalog L_2_1 Values');
 ylabel('Elution Model Fitted \Deltalog L_2_1 Values');
 box on;
end

r2_DeltaL21 = corrcoef(Delta_logL21, Delta_logL21_pred).^2;
rmse_DeltaL21 = sqrt(sum((Delta_logL21_pred-Delta_logL21).^2)./length(Delta_logL21));

if strcmp(prompt_output,'verbose')
 disp('The r^2 and RMSE values of fitted Delta logL_21 values by eq 4 are:'); 
 disp([r2_DeltaL21(1,2) rmse_DeltaL21]);
end

% Save the LSER-determined logL1 and DeltalogL21 values of the calibration analytes to file
u12_calib = [logL1 Delta_logL21];
savefilename = strcat('../',output_path,'logL1_DeltalogL21_calibration.dat');
save(savefilename,'u12_calib','-ascii');

% Save the model-fitted logL1 and DeltalogL21 values of the calibration analytes to file
u12_pred = [logL1_pred Delta_logL21_pred];
savefilename = strcat('../',output_path,'logL1_DeltalogL21_fitted.dat');
save(savefilename,'u12_pred','-ascii');

