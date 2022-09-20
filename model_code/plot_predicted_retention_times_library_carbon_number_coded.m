function plot_predicted_retention_times_library_carbon_number_coded(calibration_analytes, library, t1_lib, t2_lib, modulation_period, prompt_output)

% Plot simulated retention time data of library constituents having predicted first-dimension retention times >2 min
CN_lib_gt2min = library.CN(find(t1_lib > 2));
t1_lib_gt2min = t1_lib(find(t1_lib > 2));
t2_lib_gt2min = t2_lib(find(t1_lib > 2));

f_CN = figure;
cm = colormap(f_CN, 'jet');
custom_colormap = downsample(cm,12);

colorlist_CN = [1:length(unique(CN_lib_gt2min))+1];

colororder(custom_colormap);
ax = axes;
for CN_index = unique(CN_lib_gt2min')
 if CN_index > 0
  colorindex_CN = colorlist_CN(CN_index-min(unique(CN_lib_gt2min))+1);
  hold on;
  ax.ColorOrderIndex = colorindex_CN;
  plot(t1_lib_gt2min(find(CN_lib_gt2min==CN_index)),t2_lib_gt2min(find(CN_lib_gt2min==CN_index)),'.','MarkerSize',4);
 end
end
xlabel('1^{st} Dimension Retention Time (min)');
ylabel('2^{nd} Dimension Retention Time (s)');
title('Simulated retention times of the constituent library','Color-coded by carbon number');
ylim([0 modulation_period]);
box on;

