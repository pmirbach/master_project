function [Para, Data] = load_daniel_kpts(Ctrl, Para, Data)

if any( Ctrl.comp_Daniel )
    
    do
    
end

% % load kpts_11x11.mat
% % k1 = permute(k_11x11,[2,1,3]);
% % k1 = k1(:,1:66);
% 
% load kpts_35x35.mat
% k1 = permute(kpts_35x35,[2,1,3]);
% 
% k2 = k1(1:2,:);
% Data.wk = round( k1(3,:) / min(k1(3,:)) );
% Para.BZsmall.area = 1;
% Para.symm_indices = find( Data.wk == 1 );
% 
% [Data.k] = red_to_BZ(k2);
% Para.nr.k = size(Data.k,2);
% 
% % load kpts_11x11(2).mat
% % k3 = permute(kpts,[2,1,3]);
% % Data.k = k3(1:2,:,:);
% % Para.nr.k = size(Data.k,2);
% 
% Para.BZsmall.area = min(k1(3,:));
% % Para.BZsmall.area = min(k3(3,:));
% Para.coul.pol = 1.27287195103197  / Para.BZsmall.area;        % 35 x 35
% % Para.coul.pol = 4.32776463350871  / Para.BZsmall.area;          % 11 x 11
%  
% Para.k.qmin = 0.193102996717152;                                % 35 x 35
% % Para.k.qmin = 0.656550188837845;                                % 11 x 11








%% Daniels Eigenvektoren
% load('CVec (2).mat')
% Data.Ev = CVec;
% Data.Ev = abs(real(Data.Ev)) + 1i * abs(imag(Data.Ev));

% load('CVec_35.mat')
% Data.Ev = CVec;

% load test_eig_chol.mat
% Data.Ek = ENERGY;
% Data.Ev = coeff;

% load('CVec_35_noSOC.mat')
% 
% compa = Ev_comp(Prep.Ev_noSOC, CVec(1:3,1:3,:,:));
% 
% figure
% set(gcf, 'Color', 'w');
% for ii = 1:3
%     subplot(1,3,ii)
%     imagesc(squeeze(sum(real(compa(ii,:,:)),3)))
%     colorbar
% end

% figure
% set(gcf, 'Color', 'w');
% for ii = 1:6
%     subplot(2,3,ii)
%     scatter3(kx,ky,imag(Data.Ev(1,2,:,ii)))
% end



% load('CVec_35_noSOC.mat')

% Prep.Ev_noSOC(:,1,:,:) = - Prep.Ev_noSOC(:,1,:,:);