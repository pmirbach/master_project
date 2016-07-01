function [] = compare_alex( k_P, k_A, Data_P, Data_A , art )

% k size in 1st dim, kx ky in columns
if size( k_P , 1 ) < size( k_P , 2 )
    k_P = permute( k_P , [2 1] );
end

if size( k_A , 1 ) < size( k_A , 2 )
    k_A = permute( k_A , [2 1] );
end

if size( Data_P , 1 ) < size( Data_P , 2 )
    Data_P = permute( Data_P , [2 1] );
end

if size( Data_A , 1 ) < size( Data_A , 2 )
    Data_A = permute( Data_A , [2 1] );
end



k_P = round( k_P , 6 ) ;
k_A = round( k_A , 6 ) ;

[k_P_sort, index_P] = sortrows( k_P );
[k_A_sort, index_A] = sortrows( k_A );

if any( any( k_P_sort - k_A_sort ) )
    error('k points do not match!')
end

Data_P_sort = Data_P( index_P );
Data_A_sort = Data_A( index_A );

if strcmp(art,'rel')
    plotdata = Data_P_sort ./ Data_A_sort - 1;   
elseif strcmp(art,'abs')
    plotdata = ( Data_P_sort - Data_A_sort );
else
    error('Use abs or rel!')
end

color_sc =  ( (plotdata - min(plotdata)) / max(plotdata - min(plotdata) )  ) ;
C = [ zeros( length(color_sc) ,1) , color_sc , zeros(length(color_sc),1) ];                 % green

scatter3( k_P_sort(:,1),k_P_sort(:,2), plotdata , 80, C , 'filled' )