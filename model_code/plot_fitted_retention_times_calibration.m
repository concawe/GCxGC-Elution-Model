function plot_fitted_retention_times_calibration(calibration_analytes, t1_calib, t2_calib, modulation_period, prompt_output)


% The vectors t1_meas and t2_meas represent the measured retention time values for the set of calibration analytes.
t1_meas = calibration_analytes.numdata.data(:,1);
t2_meas = calibration_analytes.numdata.data(:,2);

figure;
ax = axes;

% To visualize the colors that correspond to the integers in the colorlist, search for "Why are plot lines different colors"
% in the Matlab documentation. The default color order is:
% 1 dark blue; 2 orange; 3 gold; 4 purple; 5 green; 6 cyan; 7 dark red
colorlist = load('../model_parameters/colorlist_class.dat');

% The class list ordering is provided in the function create_class_index.
class_number_calib = create_class_index(calibration_analytes.class);
for class_ind = unique(class_number_calib')
 colorindex = colorlist(class_ind);
 hold on;
 ax.ColorOrderIndex = colorindex;
 plot(t1_meas(find(class_number_calib==class_ind)),t2_meas(find(class_number_calib==class_ind)),'^');
 hold on;
 ax.ColorOrderIndex = colorindex;
 plot(t1_calib(find(class_number_calib==class_ind)),t2_calib(find(class_number_calib==class_ind)),'o');
 analyte_ind = find(class_number_calib==class_ind);
 for count = 1:length(analyte_ind)
  hold on;
  ax.ColorOrderIndex = colorindex;
  plot([t1_meas(analyte_ind(count)) t1_calib(analyte_ind(count))],[t2_meas(analyte_ind(count)) t2_calib(analyte_ind(count))],'-');
 end
 clear analyte_ind;
end

plot(t1_meas(find(class_number_calib==11)),t2_meas(find(class_number_calib==11)),'k^');
plot(t1_calib(find(class_number_calib==11)),t2_calib(find(class_number_calib==11)),'ko');
alkane_ind = find(class_number_calib==11);

for count = 1:length(alkane_ind)
 plot([t1_meas(alkane_ind(count)) t1_calib(alkane_ind(count))],[t2_meas(alkane_ind(count)) t2_calib(alkane_ind(count))],'k-');
end
hold off;

box on;
xlabel('1^{st} Dimension Retention Time (min)')
ylabel('2^{nd} Dimension Retention Time (s)')
title('Comparison of model-fitted retention times with measured values');
ylim([0 modulation_period]);

rmse_t1 = sqrt(mean((t1_meas-t1_calib).^2));
rmse_t2 = sqrt(mean((t2_meas-t2_calib).^2));
maxdev_t1 = max(abs(t1_meas-t1_calib));
maxdev_t2 = max(abs(t2_meas-t2_calib));
r_t1 = corrcoef(t1_meas,t1_calib);
r_t2 = corrcoef(t2_meas,t2_calib);
r2_t1 = r_t1(1,2)^2;
r2_t2 = r_t2(1,2)^2;

if or(strcmp(prompt_output,'normal'), strcmp(prompt_output,'verbose'))
 disp('The r^2 and RMSE values (min) of fitted retention time 1 values for calibration analytes are:'); 
 disp([r2_t1 rmse_t1]);
 disp('The r^2 and RMSE values (s) of fitted retention time 2 values for calibration analytes are:'); 
 disp([r2_t2 rmse_t2]);
end

