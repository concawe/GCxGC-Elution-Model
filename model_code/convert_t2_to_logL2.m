function [logL2_To] = convert_t2_to_logL2(x,calibration_analytes);

% Determine the value of t_2,i_star given by eq S2 in Arey et al (2022)
t2_star = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,3),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

% Determine the value of N_i_star given by eq S1 in Arey et al (2022)
N_star = interp1(calibration_analytes.alkanes(:,2),calibration_analytes.alkanes(:,1),calibration_analytes.numdata.data(:,1),'linear','extrap'); 

% Determine logL2 at 121 C for n-alkanes by eq 5 of Arey et al (2022)
% Coefficients are set based on logL values predicted by Abraham eq given in
% Abraham 1999 (review) for OV 17 (Table 13) at 120 C.
logL2ref_To = calibration_analytes.logL2_alkane_Nstar_parms(1)*N_star + calibration_analytes.logL2_alkane_Nstar_parms(2);

% Determine logL2 at 121 C according to eq 4 in Arey et al (2022)
logL2_To = logL2ref_To + log10((calibration_analytes.numdata.data(:,2) - x(1))./(t2_star - x(1)));
