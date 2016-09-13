clear variables
close all


load Final_MoS2_4band.mat
load Final_MoS2_6band.mat

figure
set(gcf,'color','w');

plot(Ergebnis_MoS2.E,Ergebnis_MoS2.alpha)
hold on
plot(Ergebnis_MoS2_6B.E,Ergebnis_MoS2_6B.alpha,'r--')

set(gca,'fontsize',18)

title('MoS_2 - 4 Band vs 6 Band')
xlabel('Energie E - Egap in meV')
ylabel('Absorptionskoeffizient \alpha')
legend('4 Band','6 Band','location','northwest')