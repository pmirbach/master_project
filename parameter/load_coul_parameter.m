function [coul_parameter, parameter_names] = load_coul_parameter( Ctrl )

RawData = importdata('Coul_parameters.txt','\t');

material_list = RawData.textdata(2:end,1);

parameter_names{1} = RawData.textdata(2:7,2);
parameter_names{2} = RawData.textdata(1,3:end);

material_index = ismember(material_list,Ctrl.material);

coul_parameter = RawData.data(material_index,:);
