function plot_predicted_retention_times_library_class_single_CN(calibration_analytes, library, t1_lib, t2_lib, modulation_period, CN_to_visualize, prompt_output);

% Plot simulated retention time data of library constituents having predicted first-dimension retention times >2 min
t1_lib_gt2min = t1_lib(find(t1_lib > 2));
t2_lib_gt2min = t2_lib(find(t1_lib > 2));
class_number_lib_gt2min = create_class_index(library.class(find(t1_lib > 2)));
CN_lib_gt2min = library.CN(find(t1_lib > 2));

array_of_inpolygon_booleans = zeros(length(t1_lib_gt2min),length(library.HCBs));

colorlist = load('../model_parameters/colorlist_class.dat');

figure;
ax = axes;
polygon_index = 1;
colorindex = 0;
for class_index = 1:length(unique(library.class))
 colorindex = colorindex + 1;
 for CN_index = CN_to_visualize
  t1_class_CN_gt2min = t1_lib_gt2min(class_number_lib_gt2min==class_index & CN_lib_gt2min==CN_index);
  t2_class_CN_gt2min = t2_lib_gt2min(class_number_lib_gt2min==class_index & CN_lib_gt2min==CN_index);
  if length(t1_class_CN_gt2min) < 3
   hold on;
   ax.ColorOrderIndex = colorlist(colorindex);
   plot(t1_class_CN_gt2min,t2_class_CN_gt2min,'.','MarkerSize',4);
   if class_index == 11
    plot(t1_class_CN_gt2min,t2_class_CN_gt2min,'k.','MarkerSize',4);  
   end
  else
   t1t2 = unique([t1_class_CN_gt2min t2_class_CN_gt2min],'rows');
   DT = delaunayTriangulation(t1t2(:,1),t1t2(:,2));
   k = convexHull(DT);
   array_of_inpolygon_booleans(:,polygon_index) = inpolygon(t1_lib_gt2min,t2_lib_gt2min,DT.Points(k,1),DT.Points(k,2));
   polygon_index = polygon_index + 1;
   hold on;
   ax.ColorOrderIndex = colorlist(colorindex);
   plot(t1t2(:,1),t1t2(:,2),'.','MarkerSize',4);
   hold on;
   ax.ColorOrderIndex = colorlist(colorindex);
   plot(DT.Points(k,1),DT.Points(k,2));
  end
 end
end
hold off;
xlabel('1^{st} Dimension Retention Time (min)')
ylabel('2^{nd} Dimension Retention Time (s)')
ylim([0 modulation_period]);
title('Simulated retention times of the constituent library','for a single carbon number, color-coded by class');
box on;
