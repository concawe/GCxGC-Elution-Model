function [logL1_Delta_logL21, beta] = determine_logL1_Delta_logL21(ESABVL, beta_input, output_path, prompt_output);

% Determine the number of structures in the library
N_structures = length(ESABVL);

% Set up an array of the LSER parameters library parameters which can be multiplied by the LSER coefficients
LSERparms = [ESABVL ones(N_structures,1)];

% Set the LSER coefficients for column stationary phase 1, in the order [e s a b v l c]
% The default values are SE-30 coefficients taken from the 1999 review by Abraham
LSERcoeffs_1 = [0.024 0.190 0.125 0.0 0.0 0.498 -0.194];

% Set the LSER coefficients for column stationary phase 2, in the order [e s a b v l c]
% The default values are OV-17 coefficients taken from the 1999 review by Abraham 
LSERcoeffs_2 = [0.071 0.653 0.263 0.0 0.0 0.518 -0.372];

% Determine the LSER model gas-stationary partition constant values for columns 1 and 2
logL1 = sum(LSERparms.*(ones(N_structures,1)*LSERcoeffs_1),2);
logL2 = sum(LSERparms.*(ones(N_structures,1)*LSERcoeffs_2),2);

% Save logL1 and logL2 values of the library to file
logL1_logL2_library = [logL1 logL2];
savefilename = strcat('../',output_path,'logL1_logL2_library.dat');
save(savefilename,'logL1_logL2_library','-ascii','-tabs');

% If the user has set beta = 0, then a new value of beta will be determined by finding
% the constant which produces an orthogonal vector from logL1 and logL2 by
% Schmidt orthogonalization, where <logL1*(logL1-beta*logL2)> = 0:
if beta_input == 0
 beta = sum(logL1.*logL2)/sum(logL1.^2);
 if strcmp(prompt_output,'verbose')
  disp('The beta value determined by Gram-Schmidt orthogonalization for the Constituent Library is:');
  disp(beta);
 end
else
 beta = beta_input;
end

% The orthogonal vectors logL1 and Delta_logL21 satisfy the relationship: <logL1*Delta_logL21> = 0
Delta_logL21 = logL2-beta*logL1;

% Determine the correlation of logL1 with Delta_logL21 
r2_logL1_Delta_logL21_training_set = corrcoef(logL1-mean(logL1),Delta_logL21-mean(Delta_logL21)).^2;

if strcmp(prompt_output,'verbose')
 disp('The squared correlation (r^2) of logL1 and Delta_logL21 of the Constituent Library is:');
 disp(r2_logL1_Delta_logL21_training_set(1,2))
end

% Save the LSER-predicted logL1 and Delta_logL21 values of the library to file
logL1_Delta_logL21 = [logL1 Delta_logL21];
savefilename = strcat('../',output_path,'logL1_Delta_logL21_library.dat');
save(savefilename,'logL1_Delta_logL21','-ascii');

