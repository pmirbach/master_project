clear variables
close all

konvergenz = cell(1,4);
konvergenz{1} = importdata('Ergebnisse/Konvergenz_MoS2.txt');
konvergenz{2} = importdata('Ergebnisse/Konvergenz_MoSe2.txt');
konvergenz{3} = importdata('Ergebnisse/Konvergenz_WS2.txt');
konvergenz{4} = importdata('Ergebnisse/Konvergenz_WSe2.txt');

materials = {'MoS_2','MoSe_2','WS_2','WSe_2'};

for jj = 1:2
    figure
    set(gcf,'color','w');
    for ii = 1:2
        
        NN = ii + (jj-1) * 2;
        
        subplot(2,1,ii)
        
        title(materials{NN})
        title('test')
        
        plot(konvergenz{NN}(:,2),konvergenz{NN}(:,3),'b-x')
        hold on
        plot(konvergenz{NN}(:,2),konvergenz{NN}(:,4),'r-x')
        
        set(gca,'fontsize',18)
        xlabel('Anzahl k-Punkte in BZ')
        ylabel('Energie E-E_{Gap} in meV')
        
        legend('A-Exziton','B-Exziton','location','se')
        
    end
end


% 
% plot(N_k,E_alpha,'b-x')
% hold on
% plot(N_k,E_beta,'r-x')
% 
% set(gca,'fontsize',18)
% 
% legend('A-Exziton','B-Exziton','location','se')
% xlabel('Anzahl k-Punkte in BZ')
% ylabel('Energie E-E_{Gap} in meV')