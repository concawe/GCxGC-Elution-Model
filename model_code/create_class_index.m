function class_number = create_class_index(class_data)

class_number = zeros(length(class_data),1);

class_number(find(class_data=="iP")) = 1;
class_number(find(class_data=="mN")) = 2;
class_number(find(class_data=="dN")) = 3;
class_number(find(class_data=="polyN")) = 4;
class_number(find(class_data=="mAr")) = 5;
class_number(find(class_data=="dAr")) = 6;
class_number(find(class_data=="polyAr")) = 7;
class_number(find(class_data=="NmAr")) = 8;
class_number(find(class_data=="NdAr")) = 9;
class_number(find(class_data=="NpolyAr")) = 10;
class_number(find(class_data=="nP")) = 11;

