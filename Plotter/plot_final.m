clear variables
% close all

% load Final_MoS2_q150.mat
% load Final_MoSe2_q150.mat
% load Final_WS2_q150.mat
load Final_WSe2_q150.mat

% load Final_MoS2_q150.mat
% load Final_MoS2_q150_up.mat
% load Final_MoS2_q150_dwn.mat

%%

figure(10)
hold on
set(gcf,'color','w')
set(gcf,'units','normalized','position',[.1 .1 .8 .45])

plot( (Ergebnis.energy + Ergebnis.EGap)/1000 + Ergebnis.EGapCorr  , Ergebnis.alpha )
set(gca,'fontsize',16)
box on

xlabel('Energie in eV')
ylabel('Absoprtionskoeff. \alpha')

%%
xlim([min((Ergebnis.energy + Ergebnis.EGap)/1000 + Ergebnis.EGapCorr) max((Ergebnis.energy + Ergebnis.EGap)/1000 + Ergebnis.EGapCorr)])
% legend('\uparrow & \downarrow', '\uparrow', '\downarrow','location','northwest')
linie = plot([Ergebnis.EGap/1000 + Ergebnis.EGapCorr , Ergebnis.EGap/1000 + Ergebnis.EGapCorr ],get(gca,'YLim'),'color',0.4*[1 1 1]);