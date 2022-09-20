function [t1_pred t2_pred] = convert_logL1_DeltalogL21_to_t1_t2(logL1_DeltalogL21,alpha1_alpha2,alpha3,beta,calibration_analytes); 

logL1 = logL1_DeltalogL21(:,1);
logL2 = logL1_DeltalogL21(:,2)+beta*logL1;

% Use given logL1 values to recover Ni* values, by eq S3.
N_star = (logL1-alpha1_alpha2(2))/alpha1_alpha2(1);

% Use Ni* values to predict t1 values by way of eq 7.
t1_pred = interp1(calibration_analytes.alkanes(:,1),calibration_analytes.alkanes(:,2),N_star,'linear','extrap');

% Use predicted t1 values and alkanes data to obtain t2* values by eq S2.
t2_star = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,3),t1_pred,'linear','extrap');

% Use the t2* values and Ni* values determined above together with given log2 values to predict t2, by eq 8.
t2_pred = alpha3 + (t2_star-alpha3).*10.^(logL2-calibration_analytes.logL2_alkane_Nstar_parms(1)*N_star-calibration_analytes.logL2_alkane_Nstar_parms(2));

