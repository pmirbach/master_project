function [names, values] = load_TB_parameter_old(Ctrl)

% Empirische TB-Parameter für zu berechnendes Single-Layer Dichalcogenide

load TB_Liu_parameter
ind_mat = ismember(TB_Liu_parameter(:,1),Ctrl.material);

if strcmp(Ctrl.method,'NN')
    ind_met = 2;
elseif strcmp(Ctrl.method,'TNN')
    ind_met = 3;
end

names = TB_Liu_parameter{1,ind_met};
values = TB_Liu_parameter{ind_mat,ind_met};
