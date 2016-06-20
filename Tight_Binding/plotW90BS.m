clear all; 
clc;

addpath('./Wannier90Interface');

    % Settings ------------------------------------------------------------
    
    %   Hamiltonian
    %seed        = '/obelix/TMDC/02_Results/02_Materials/WSe2/a0_3.325A/01_WannierTB/01_DFT/02_W3d/wannier90';
    %seed        = '/obelix/TMDC/02_Results/02_Materials/WS2/a0_3.191A/01_WannierTB/02_G0W0/02_W3d/wannier90';
    %seed        = '/obelix/TMDC/02_Results/02_Materials/MoSe2/a0_3.320A/01_WannierTB/02_G0W0/02_Mo3d/wannier90';
    
%     seed        = '/obelix/TMDC/02_Results/02_Materials/MoS2_beta/a0_3.160A/01_WannierTB/02_G0W0/02_Mo3d/wannier90'
    seed        = './ab_initio/02_Materials/MoS2/a0_3.180A/01_WannierTB/02_G0W0/02_Mo3d/wannier90';
    
    %   SOC Settings
    SOCSettings.type  = 'none'; % 'first' or 'second' or 'none'
    
    SOCSettings.dampedLambda0 = 1;
    
    % WSe2
    %SOCSettings.Delta =  0.464;
    %SOCSettings.E0p1  = -3.279;
    %SOCSettings.E0m1  = -1.292;
    %SOCSettings.Em2m1 =  0.607;
    %SOCSettings.BG    =  1.5344;
    
    % WS2
    SOCSettings.Delta =  0.464;
    SOCSettings.E0p1  = -3.667;
    SOCSettings.E0m1  = -1.544;
    SOCSettings.Em2m1 =  0.655;
    SOCSettings.BG    =  3.0037;
    
    % MoSe2
    %SOCSettings.Delta =  0.200;
    %SOCSettings.E0p1  = -2.862;
    %SOCSettings.E0m1  = -1.125;
    %SOCSettings.Em2m1 =  0.462;
    %SOCSettings.BG    =  2.2715;
    
    %   choose bands num in dependence of SOC settings
    if(~strcmp(SOCSettings.type, 'none'))
        bandsNum = 6;
    else
        bandsNum = 3;
    end
    
    %   comparison BS
    %BSCompFile  = '/obelix/TMDC/01_Calculations/WSe2/02_DFT/wSOC/bands/fbands.dat';
    %BSCompFile  = '/obelix/TMDC/01_Calculations/MoSe2/02_DFT/bands/fbands.dat';
%     BSCompFile  = '/obelix/TMDC/01_Calculations/WS2/02_DFT/wSOC/bands/fbands.dat';
%     BSComp      = load(BSCompFile);
    
    %   BS Path
    stepsPerPath    = 30;
    
    line(1, 1, 1:3) = [0.0, 0.0, 0.0];  % G
    line(1, 2, 1:3) = [1/2, 0.0, 0.0];  % M
 
    line(2, 1, 1:3) = [0.5, 0.0, 0.0];  % M
    line(2, 2, 1:3) = [1/3, 1/3, 0.0];  % K
    
    line(3, 1, 1:3) = [1/3, 1/3, 0.0];  % K
    line(3, 2, 1:3) = [0, 0, 0];        % Gamma
    
    % =====================================================================
    
    % load wannier90 Data
    W90Data = loadW90Data(seed, SOCSettings);
    
    disp('Bandstructure'); % ----------------------------------------------
    
    a0ToA = 0.529177249;           % 1a0  = 0.529177249Ang
    
    % If Data is already in Ang and not in a0
    a0ToA = 1;
    
    kToEv = 1.973269630783068e+03; % 1/Ang = 1.97... c^-1 hq^-1 eV
    m0 = 510.99906 * 10^3;         % electron mass in eV

    figure(1);
    clf;
    hold on;

    pp = 1;
    steps = stepsPerPath;
    kpNorm(1) = 0;
    for ll = 1 : length(line(:,1,1))
        
        disp('new line ...');
        for jj = 1 : steps
            
            if(pp > 1)
                kpNorm(pp) = kpNorm(pp-1) + norm(stepPlot);
            end
            
            a(1:3) = line(ll, 1, 1:3);
            b(1:3) = line(ll, 2, 1:3);
            
            a(1:3) = W90Data(1).kLat' * a(1:3)';
            b(1:3) = W90Data(1).kLat' * b(1:3)';
            
            stepPlot = (b - a) / steps;
            
            step(1:3) = (line(ll, 2, 1:3) - line(ll, 1, 1:3)) / steps;

            tmp2(1:3) = line(ll, 1, 1:3);
            kp(pp, 1:3) = tmp2 + (jj-1) .* step;
            
            pp = pp + 1;
            
        end
        
    end
    
    %get hamiltonian for every kk point
    HH = getW90Hamiltonian(W90Data, kp);
    
    for ii = 1 : length(kp)
       
        tmp(1:bandsNum, 1:bandsNum) = HH(ii, 1:bandsNum, 1:bandsNum);
        
        [ef, ev] = eig(tmp);

        [ev, I] = sort(diag(real(ev)));
        ef = ef(:, I);
        
        values(ii, 1:bandsNum) = sort(ev);
        
    end
    
%     plot(BSComp(:,1), BSComp(:,2), 'linestyle', 'none', 'marker', '.', 'markersize', 10);
    plot(kpNorm, values(:, 1:bandsNum), 'color', 'red', 'linewidth', 2);
    
    ylim([ min(min(values(:, :)))-0.2, max(max(values(:, :)))+0.2 ]); 
    
    %fid = fopen('BS.dat', 'w');
    %for ii = 1:length(values(1,:))
    %    for kk = 1:length(values(:,1))
    %        fprintf(fid, '%f\t\t%f \n', kpNorm(kk), values(kk, ii));
    %    end
    %    fprintf(fid, ' \n');
    %end
    %fclose(fid);


