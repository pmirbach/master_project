clear variables
close all


load Final_MoS2_4band_up.mat
load Final_MoS2_4band_down.mat
load Final_WSe2_4band_up.mat
load Final_WSe2_4band_down.mat

figure
set(gcf,'color','w');


materials = {'MoS_2','WSe_2'};

subplot(2,1,1)
plot(Ergebnis_MoS2_up.E,Ergebnis_MoS2_up.alpha)
hold on
plot(Ergebnis_MoS2_down.E,Ergebnis_MoS2_down.alpha)

subplot(2,1,2)
plot(Ergebnis_WSe2_up.E,Ergebnis_WSe2_up.alpha)
hold on
plot(Ergebnis_WSe2_down.E,Ergebnis_WSe2_down.alpha)


for ii = 1:2
    subplot(2,1,ii)
    legendstr = [materials{ii}];
    set(gca,'fontsize',18)
    legend([materials{ii},' \uparrow'],[materials{ii},' \downarrow'],'location','northwest')
    ylabel('\alpha')
end

xlabel('Energie E - Egap in meV')

