function W90Data = loadW90Data(seed, trans_me , SOCSettings)

% settings ================================================================

w90path = fullfile( seed , '01_WannierTB','02_G0W0',['02_' , trans_me , '3d'] );

% wannier90 input file (to get the basis labeling and lattice vectors)
wannier90File = fullfile( w90path , 'wannier90.win' );

% filename of Hamiltonian in Wannierbasis (wannier90 output)
hamiltonianFile = fullfile( w90path , 'wannier90_hr.dat' );

% filename of details txt file for orbital order and band gap correction
detailsFile = fullfile( w90path , 'Details.txt' );


% basis / EGap ============================================================

fid = fopen(detailsFile, 'r');

tline = fgetl(fid);
while ischar( tline )
    
    test = strfind(tline, 'orbital labels');
    if ~isempty( test )
        
        fgetl(fid);
        goon = true;
        bb = 1;
        while goon
            
            tline = fgetl(fid);
            
            if ~all(isspace(tline)) && ischar( tline )
                
                tline_split = strsplit(tline,' ');
                
                W90Data(1).basis{bb} = tline_split{end};
                bb = bb + 1;
            else
                goon = false;
            end
        end
    end
    
    test = strfind(tline, 'EGapAct');
    if ~isempty( test )
        
        tline_split = strsplit(tline,' ');
        W90Data(1).EGapAct = str2double( tline_split{end} );
        
    end

    test = strfind(tline, 'EGapCorr');
    if ~isempty( test )
        
        tline_split = strsplit(tline,' ');
        W90Data(1).EGapCorr = str2double( tline_split{end} );
        
    end
    
    tline = fgetl(fid);
end

fclose(fid);


% preojection / lattice ===================================================

% open file
fid = fopen(wannier90File, 'r');

%   read each line untile
%       - 'begin projections'
%       - 'begin unit_cell_cart'
%   are found
tline = fgetl(fid);
while ischar(tline)
    
    % get basis
    test = strfind(tline, 'begin projections');
    if(~isempty(test))
        
        goon = true;
        bb = 1;
        while(goon == true)
            
            tline = fgetl(fid);
            test = strfind(tline, 'end projections');
            if(~isempty(test))
                goon = false;
            else
                W90Data(1).prejection{bb} = strtrim(tline);
                bb = bb + 1;
            end
            
        end
        
    end
    
    % get lattice
    test = strfind(tline, 'begin unit_cell_cart');
    if(~isempty(test))
        
        % maybe there is an "bohr" statement following
        tline = fgetl(fid);
        test = strfind(tline, 'bohr');
        if(~isempty(test))
            nextLine = fgetl(fid);
        else
            nextLine = tline;
        end
        
        % the next three lines will include the lattice vectors
        W90Data(1).rLat(1,1:3) = sscanf(nextLine, '%f %f %f');
        W90Data(1).rLat(2,1:3) = sscanf(fgetl(fid), '%f %f %f');
        W90Data(1).rLat(3,1:3) = sscanf(fgetl(fid), '%f %f %f');
        
        % calculate k lattice vectors
        W90Data(1).kLat(1,1:3) = ...
            2 * pi ...
            * cross(W90Data(1).rLat(2,1:3), W90Data(1).rLat(3,1:3)) ...
            / ( W90Data(1).rLat(1,1:3) ...
            * (cross(W90Data(1).rLat(2,1:3), W90Data(1).rLat(3,1:3) ))');
        W90Data(1).kLat(2,1:3) = ...
            2 * pi ...
            * cross(W90Data(1).rLat(3,1:3), W90Data(1).rLat(1,1:3)) ...
            / ( W90Data(1).rLat(1,1:3) ...
            * (cross(W90Data(1).rLat(2,1:3), W90Data(1).rLat(3,1:3) ))');
        W90Data(1).kLat(3,1:3) = ...
            2 * pi ...
            * cross(W90Data(1).rLat(1,1:3), W90Data(1).rLat(2,1:3)) ...
            / ( W90Data(1).rLat(1,1:3) ...
            * (cross(W90Data(1).rLat(2,1:3), W90Data(1).rLat(3,1:3) ))');
        
    end
    
    tline = fgetl(fid);
    
end

fclose(fid);


% details =================================================================


% hamilatonian ============================================================

% open and read file
disp(['reading ', hamiltonianFile, ' ...']);
fid = fopen(hamiltonianFile, 'r');

%   1. line: date and time at which the file was created
W90Data.dateTime = fgetl(fid);
%   2. line: number of Wannier functions
W90Data.numWann = sscanf(fgetl(fid), '%f');
%   3. line: the number of Wigner-Seitz grid-points
W90Data.numWSGridPoints = sscanf(fgetl(fid), '%f');
%   The next block of nrpts integers gives the degeneracy of each
%   Wigner-Seitz grid point, with 15 entries per line.
numLines = ceil(W90Data.numWSGridPoints / 15);
for nn = 1 : numLines
    tmp = sscanf(fgetl(fid), '%i');
    W90Data.degeneracy((nn*15)-14:(nn*15)-15+length(tmp)) = tmp;
end

fclose(fid);

%   the remaining num wann^2×nrpts lines each contain, respectively, the
%   components of the vector R in terms of the lattice vectors {Ai}, the
%   indices m and n, and the real and imaginary parts of the Hamiltonian
%   matrix element H(R)_mn in the WF basis
hamiltonianRaw = importdata(hamiltonianFile, ' ', 3+numLines);
W90Data.hamiltonianRaw = hamiltonianRaw.data;

W90Data.SOCSettings = SOCSettings;


W90Data = orderfields(W90Data, [ 7 4 1 5 6 8 9 10 11 2 3 12 ] );

end
