function [ sample ] = get_sample( parameter , quantity )

nr_parameter = size( parameter , 2 );
sample = zeros( quantity , nr_parameter );


for jj = 1:nr_parameter
    
    sample(:,jj) = round( rand( quantity , 1 ) * ( parameter( jj ) - 1 ) ) + 1;
    
end
