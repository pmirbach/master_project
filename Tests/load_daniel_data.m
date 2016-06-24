function [Para, Data, MoS2] = load_daniel_data(Ctrl, Para, Data, file)

if Ctrl.cmp.load
    
    load(file)
    
    if Ctrl.cmp.use_k
        
        k1 = permute(MoS2.kpts,[2,1,3]);
        Data.k = k1(1:2,:,:);
        
        Data.wk = round( k1(3,:) / min(k1(3,:)) );
       
        Para.symm_indices = find( Data.wk == 1 );
        Para.BZsmall.area = min(k1(3,:));
        Para.nr.k = MoS2.npts;
        Para.coul.pol = MoS2.coulpol;
        Para.k.qmin = norm( Data.k(:,1) - Data.k(:,2) ) / 2;
        
        if round ( sum(Data.wk) * Para.BZsmall.area * 1e3 ) ~= round( Para.BZ.area * 1e3 )
            warning('Integrated weights do not agree with area of BZ!')
        end
        
    end
    
else
    MoS2 = [];
end





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