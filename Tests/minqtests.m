function minqtests(Parameter, Prep, k)

corners = [Parameter.symmpts{2} ...
    * [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0], [0; 0]];
[corners] = red_to_BZ(corners);

figure
hold on
set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);

colors = {'b','r','y','g','c','m'};
marker = {'h','p','x','+','^','v'};

for ii = 1:6
    in = k(3,:,ii) == 6 * Parameter.area_sBZ;
    on = k(3,:,ii) == 3 * Parameter.area_sBZ;
    symm = k(3,:,ii) == 1 * Parameter.area_sBZ;
    
    instr = strcat(colors{ii},marker{1+mod(ii,2)});
    onstr = strcat(colors{ii},marker{3+mod(ii,2)});
    symmstr = strcat(colors{ii},marker{5+mod(ii,2)});
    
    plot(corners(1,:,ii),corners(2,:,ii),'color',0.9 * [1 1 1])
    
    plot(k(1,in,ii),k(2,in,ii),instr)
    plot(k(1,on,ii),k(2,on,ii),onstr)
    plot(k(1,symm,ii),k(2,symm,ii),symmstr)
    
end
axis([min(min(k(1,:,:))) - 0.2, max(max(k(1,:,:))) + 0.2, ...
    min(min(k(2,:,:))) - 0.2, max(max(k(2,:,:))) + 0.2])
axis equal
set(gca,'FontSize',16)
title('k-mesh in BZ mit markierten äquivalenten k-Punkten')
xlabel('k_x')
ylabel('k_y')

set(gcf, 'Color', 'w');


ind = 2;
plot(k(1,ind,1),k(2,ind,1),'ko','markers',8)


% for ii = 1:6
%     plot(k(1,ind,ii),k(2,ind,ii),'ko','markers',8)
% end

q = zeros(size(k,2),6);
for ii = 1:6
    for jj = 1:size(k,2)
        q(jj,ii) = norm(k(1:2,ind,1)-k(1:2,jj,ii));
    end
end






1