function [ Coul_ME ] = read_Coulomb_data( seed )

coulomb_path = fullfile( seed , '02_CoulombME' );
coulombFile = fullfile( coulomb_path , 'FitParameters.txt' );

disp(['reading ', coulombFile, ' ...']);

fid = fopen(coulombFile, 'r');

tline = fgetl(fid);      % Read first line for while loop
while ischar( tline )
    
    test = strfind(tline, '(1.1.1.1)');
    if ~isempty( test )
        fgetl(fid);
        tline = fgetl(fid);
        para_value = textscan(tline,'%s %f');
        Coul_ME.U.gamma_quad = para_value{2};
    end
    
    test = strfind(tline, '(1.1.1.2)');
    if ~isempty( test )
        fgetl(fid);
        tline = fgetl(fid);
        para_value = textscan(tline,'%s %f');
        Coul_ME.U.gamma = para_value{2};
        
        tline = fgetl(fid);
        para_value = textscan(tline,'%s %f');
        Coul_ME.U.delta = para_value{2};
    end
    
    test = strfind(tline, '(1.1.2)');
    if ~isempty( test )
        fgetl(fid);
        for ii = 1:2
            tline = fgetl(fid);
            para_value = textscan(tline,'%s %f');
            Coul_ME.U.micro(ii) = para_value{2};
        end
    end
    
    tests = { '(1.2.1)' , '(1.2.2)' , '(1.2.3)' };
    for jj = 1:3
        test = strfind(tline, tests{jj});
        if ~isempty( test )
            fgetl(fid);
            for ii = 1:3
                tline = fgetl(fid);
                para_value = textscan(tline,'%s %f');
                Coul_ME.U_Ev(ii,jj) = para_value{2};
            end
        end
    end
    
    test = strfind(tline, '(2.1.1.1)');
    if ~isempty( test )
        fgetl(fid);
        tline = fgetl(fid);
        para_value = textscan(tline,'%s %f');
        Coul_ME.eps.inf = para_value{2};
        
        tline = fgetl(fid);
        para_value = textscan(tline,'%s %f');
        Coul_ME.eps.L_2d_macro = para_value{2};
    end
    
    test = strfind(tline, '(2.1.1.2)');
    if ~isempty( test )
        fgetl(fid);
        for ii = 1:5
            tline = fgetl(fid);
            para_value = textscan(tline,'%s %f');
            Coul_ME.eps.Resta(ii) = para_value{2};
        end
    end
    
    test = strfind(tline, '(2.1.2)');
    if ~isempty( test )
        fgetl(fid);
        for ii = 1:2
            tline = fgetl(fid);
            para_value = textscan(tline,'%s %f %f');
            Coul_ME.eps.micro(ii) = para_value{end};
        end
    end    
    
    tline = fgetl(fid);
end

fclose(fid);


Coul_ME.U.gamma_quad = Coul_ME.U.gamma_quad * 0.1;      % to nm
Coul_ME.U.gamma = Coul_ME.U.gamma * 0.1;                % to nm
Coul_ME.U.delta = Coul_ME.U.delta * 0.01;               % to nm^2

Coul_ME.U.micro = Coul_ME.U.micro * 1000;               % to meV

Coul_ME.eps.L_2d_macro = Coul_ME.eps.L_2d_macro * 0.1;  % to nm
Coul_ME.eps.Resta(1) = Coul_ME.eps.Resta(1) * 100;      % to 1/nm^2
Coul_ME.eps.Resta(3) = Coul_ME.eps.Resta(3) * 0.1;      % to nm
Coul_ME.eps.Resta(4) = Coul_ME.eps.Resta(4) * 0.1;      % to nm

Coul_ME.U_Ev = Coul_ME.U_Ev( [1 3 2], : );              % same basis order as TB
