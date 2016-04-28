function [tb_parameter, parameter_names] = load_TB_parameter( Ctrl )

RawData = importdata('TB_Liu_parameter.txt');

material_list = RawData.textdata(1,2:end);
parameter_names = RawData.textdata(2:end,1);

material_index = ismember(material_list,Ctrl.material);

tb_parameter = RawData.data(:,material_index);
