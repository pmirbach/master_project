clear variables
close all


load Final_MoS2_4band.mat
load Final_MoSe2_4band.mat
load Final_WS2_4band.mat
load Final_WSe2_4band.mat

figure
set(gcf,'color','w');

title('MoS_2 - 4 Band vs 6 Band')

materials = {'MoS_2','MoSe_2','WS_2','WSe_2'};

subplot(4,1,1)
plot(Ergebnis_MoS2.E,Ergebnis_MoS2.alpha)
subplot(4,1,2)
plot(Ergebnis_MoSe2.E,Ergebnis_MoSe2.alpha)
subplot(4,1,3)
plot(Ergebnis_WS2.E,Ergebnis_WS2.alpha)
subplot(4,1,4)
plot(Ergebnis_WSe2.E,Ergebnis_WSe2.alpha)


for ii = 1:4
    subplot(4,1,ii)
    legendstr = [materials{ii}, ' - 4 Band'];
    set(gca,'fontsize',18)
    legend(legendstr,'location','northwest')
    ylabel('\alpha')
end

xlabel('Energie E - Egap in meV')

