clear variables
close all
clc

Ctrl.material = 'MoS2';
Ctrl.lattice_constant = 0.318;
a0_form = num2str( 10 * Ctrl.lattice_constant , '%.3f' );
seed = fullfile('Tight_Binding','ab_initio','02_Materials',Ctrl.material,['a0_',a0_form, 'A']);

coulomb_path = fullfile( seed , '02_CoulombME' );

details_file = fullfile( coulomb_path , 'FitParameters.txt' );



fid = fopen(details_file, 'r');


tline = fgetl(fid);      % Read first line for while loop

while ischar( tline )
    %     disp( tline )
    
%     test = strfind(tline, 'orbital labels');
%     if ~isempty( test )
%         
%         fgetl(fid);
%         goon = true;
%         bb = 1;
%         while goon
%             
%             tline = fgetl(fid);
%             
%             if ~all(isspace(tline)) && ischar( tline )
%                 
%                 tline_split = strsplit(tline,' ');
%                 
%                 W90Data(1).basis{bb} = tline_split{end};
%                 bb = bb + 1;
%             else
%                 goon = false;
%             end
%         end
%     end
%     
%     test = strfind(tline, 'EGapAct');
%     if ~isempty( test )
%         
%         tline_split = strsplit(tline,' ');
%         W90Data(1).EGapAct = str2double( tline_split{end} );
%         
%     end
%     
%     
%     test = strfind(tline, 'EGapCorr');
%     if ~isempty( test )
%         
%         tline_split = strsplit(tline,' ');
%         W90Data(1).EGapCorr = str2double( tline_split{end} );
%         
%     end
    
    %     test = strfind(tline, 'band gap');
%     if ~isempty( test )
%         
%         fgetl(fid);
%         goon = true;
%         while goon
%             
%             tline = fgetl(fid);
%             if ~all(isspace(tline)) && ischar( tline )
%                 disp( tline )
%                 
%                 1
%             elseif tline == -1
%                 goon = false;
%             else
%                 goon = false;
%             end
%         end
%     end
    
    tline = fgetl(fid);
end




fclose(fid);