function [handle_coulomb_up, handle_coulomb_down] = plot_coulomb(Ctrl,Para,k,V,index)

if Ctrl.plot.coul == 1
    
    d = 0;
    
    handle_coulomb_up = figure;
    set(handle_coulomb_up, ...
        'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    set(handle_coulomb_up, 'Color', 'w');
    for ii = 1:9
        subplot(3,3,ii)
        
        color_sc = ( V(index,:,ii+d).' / max( V(index,:,ii+d) ) );
        
        C = [ zeros(Para.nr.k,1) , color_sc , zeros(Para.nr.k,1) ];     % green
        
        scatter3(k(1,:,1),k(2,:,1),V(index,:,ii+d)',12,C)
        
        title(num2str(Para.coul_indices(ii+d,1:2)))
    end
    
    
    d = 9;
    
    handle_coulomb_down = figure;
    set(handle_coulomb_down, ...
        'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    set(handle_coulomb_down, 'Color', 'w');
    for ii = 1:9
        subplot(3,3,ii)
        
        color_sc = ( V(index,:,ii+d).' / max( V(index,:,ii+d) ) );
        
        C = [ zeros(Para.nr.k,1) , color_sc , zeros(Para.nr.k,1) ];     % green
        
        scatter3(k(1,:,1),k(2,:,1),V(index,:,ii+d)',12,C)
        
        title(num2str(Para.coul_indices(ii+d,1:2)))
    end
    
else
    handle_coulomb_up = [];
    handle_coulomb_down = []; 
end