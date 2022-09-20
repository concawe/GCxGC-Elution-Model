clear;
clear functions;
close all;
clc;

% -----------------------------------------------------------------------------
% ALL PROGRAM PARAMETERS THAT ARE TYPICALLY SET BY THE USER ARE SHOWN HERE.

% MODEL INPUT AND CALIBRATION DATA

% Set the input file path
input_path = 'user/input/';

% Indicate the name of the input file that contains quantitative data about the n-alkanes.
% This file should be located in the folder indicated by input_path.
alkanes_data_file = 'alkane_data_progA.dat'; 

% Indicate the name of the input file that contains quantitative data about the calibration analytes.
% This file should be located in the folder indicated by input_path.
calibration_data_file = 'calibration_data_progA.dat';

% Indicate the name of the input file that contains text data which describe the classes of the 
% calibration analytes.
% This file should be located in the folder indicated by input_path.
calibration_classes_file = 'calibration_classes_progA.dat';

% Set the modulation period (units of seconds)
modulation_period = 12.5;

% The default value of beta is 1.07, based on the parameterization described by
% Arey et al. (2022). If the user wishes to re-determine the value of beta from 
% a user-provided constituent library, set beta = 0 here.
beta_input = 1.07;

% HYDROCARBON LIBRARY PARAMETERS

% Indicate the file name for the UFZ-LSER-estimated values of E, S, A, B, V, and L (floating point) 
% for the Hydrocarbon Library. This library file should be located in the model_parameters folder.
% The default path points to an abridged library (971 structures) which is provided on the github 
% repository. The corresponding files for the complete library can be obtained from Concawe upon 
% request at the contact address provided in the user documentation.
library_ESABVL_file = 'mini_constituent_library_ESABVL.dat';

% Indicate the file name for the assigned values of class (text) and carbon number (integers) 
% for the Hydrocarbon Library. This library file should be located in the model_parameters folder.
% The default path points to an abridged library (971 structures) which is provided on the github 
% repository. The corresponding files for the complete library can be obtained from Concawe upon 
% request at the contact address provided in the user documentation.
library_class_CN_file = 'mini_constituent_library_carbon_number_class.dat';

% MODEL OUTPUT PARAMETERS

% Set the output file path
output_path = 'user/output/';

% Set generic_plot_flag to a value of 1 to see plots of basic output from the GCxGC elution model
% Set to a value of 0 to suppress plots
generic_plot_flag = 1;

% Visualize the simulated retention times of only a subset of constituent library members having 
% a single carbon number
% Set to a value of 0 to suppress the plot
single_CN_plot_flag = 1;

% If the user has set single_CN_plot_flag = 1, then the carbon number value must be indicated here
single_CN_value = 20;

% Set Matlab console output level. Choose: 'minimal' or 'verbose'.
% The default value is 'verbose'
prompt_output = 'verbose';

% END OF USER INPUT SECTION
% -----------------------------------------------------------------------------

cd('../model_code');

addpath .
main;

cd('../user');
