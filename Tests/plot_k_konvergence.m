clear variables
close all
clc



N_k = [ 109 166 235 316 409 514 631 760 901 1054 1219 ...
    1396 1585 1786 1999 2224 2461 2710 2971 3244 3529 3826 ];

E_alpha = -[ 867 781 730 697 676 661 649 641 635 630 626 ...
    623 621 619 617 616 615 614 613 613 612 612 ];

E_beta = -[ 731 645 595 562 541 525 514 506 499 495 491 ...
    488 485 483 482 481 479 479 478 477 477 477 ];

figure
set(gcf,'color','w');

plot(N_k,E_alpha,'b-x')
hold on
plot(N_k,E_beta,'r-x')

set(gca,'fontsize',18)

legend('A-Exziton','B-Exziton','location','se')
xlabel('Anzahl k-Punkte in BZ')
ylabel('Energie E-E_{Gap} in meV')