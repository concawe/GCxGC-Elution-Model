function plot_predicted_retention_times_library_class_coded(calibration_analytes, library, t1_lib, t2_lib, modulation_period, prompt_output)

% Plot simulated retention time data of library constituents having predicted first-dimension retention times >2 min
t1_lib_gt2min = t1_lib(find(t1_lib > 2));
t2_lib_gt2min = t2_lib(find(t1_lib > 2));
class_number_lib_gt2min = create_class_index(library.class(find(t1_lib > 2)));

colorlist = load('../model_parameters/colorlist_class.dat');

figure;
ax = axes;

colorindex = 0;
for class_index = 1:length(unique(library.class))
 colorindex = colorindex + 1;
 t1_class_gt2min = t1_lib_gt2min(class_number_lib_gt2min==class_index);
 t2_class_gt2min = t2_lib_gt2min(class_number_lib_gt2min==class_index);
 hold on;
 ax.ColorOrderIndex = colorlist(colorindex);
 plot(t1_class_gt2min,t2_class_gt2min,'.','MarkerSize',4);
 if class_index == 11
  plot(t1_class_gt2min,t2_class_gt2min,'k.','MarkerSize',4);  
 end
end
hold off;
xlabel('1^{st} Dimension Retention Time (min)')
ylabel('2^{nd} Dimension Retention Time (s)')
ylim([0 modulation_period]);
title('Simulated retention times of the constituent library','Color-coded by class');
box on;
