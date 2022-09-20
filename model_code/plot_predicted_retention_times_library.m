function plot_predicted_retention_times_library(calibration_analytes, t1_calib, t2_calib, t1_lib, t2_lib, modulation_period, prompt_output)

% Plot simulated retention time data of library constituents having predicted first-dimension retention times >2 min
t1_lib_gt2min = t1_lib(find(t1_lib > 2));
t2_lib_gt2min = t2_lib(find(t1_lib > 2));

figure;
plot(t1_lib_gt2min,t2_lib_gt2min,'k.','MarkerSize',4);
hold on;
plot(t1_calib,t2_calib,'ro');
hold off;
xlabel('1^{st} Dimension Retention Time (min)')
ylabel('2^{nd} Dimension Retention Time (s)')
title('Simulated retention times of the constituent library');
ylim([0 modulation_period]);
box on;
