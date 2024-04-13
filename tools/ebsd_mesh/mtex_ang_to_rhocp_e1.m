clc;clear;
%% Conversion of .ang file to rhocp-mesh
% example 1: dummy_scan with Ni - phase 1 (p1) and AlNi3 - phase (p2)

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [3.48 3.48 3.48], 'mineral', 'Ni', 'color', 'red'),...
  crystalSymmetry('m-3m', [3.52 3.52 3.52], 'mineral', 'Al Ni3', 'color', 'green')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

% which files to be imported
fname = 'example_1_file.ang';

ebsd = EBSD.load(fname,CS,'interface','ang','convertEuler2SpatialReferenceFrame','setting 1');
[grains,ebsd.grainId] = calcGrains(ebsd);
figure(1)
plot(grains); % phase map
exportgraphics(gcf,'phase_map_example_1.png','Resolution',300);

ipfKey = ipfColorKey(grains('Ni'));
ipfKey.inversePoleFigureDirection = zvector;
color = ipfKey.orientation2color(grains('Ni').meanOrientation);
figure(2)
plot(grains('Ni'),color);
figure(3)
plot(ipfHSVKey(ebsd('Ni').CS))

for ia = 1:size(ebsd.rotations.phi1,1)
    if (ebsd.phaseId(ia)==2) 
        sym_id(ia,1) = 225;
    else
        sym_id(ia,1) = 221;
    end
end

x_min = ebsd.prop.x(1); % start X co-ordinate
y_min = ebsd.prop.y(1); % start Y co-ordinate
z_min = 0.0; % start Z co-ordinate

x_step = 0.5 ; % step size in X direction
y_step = 0.5 ; % step size in Y direction
z_step = 0.5 ; % one element along Z direction

x_max = ebsd.prop.x(size(ebsd.rotations.phi1,1)); % end X co-ordinate
y_max = ebsd.prop.y(size(ebsd.rotations.phi1,1)); % end Y co-ordinate
z_max = 0.0; % one element along Z-direction

x_dim = (x_max - x_min + x_step) / x_step; % dimension - X
y_dim = (y_max - y_min + y_step) / y_step; % dimension - Y
z_dim = (z_max - z_min + z_step) / z_step; % dimension - Z

array = [ebsd.rotations.phi1 * 180/pi  ebsd.rotations.Phi * 180/pi ebsd.rotations.phi2 * 180/pi ... 
    ebsd.prop.x   abs(max(ebsd.prop.y)-ebsd.prop.y)  zeros(size(ebsd.rotations.phi1,1),1) ...
    ebsd.grainId  ebsd.phaseId-1  sym_id];

phase_1_name = 'p1';
phase_1_sym = 225;
phase_2_name = 'p2';
phase_2_sym = 221;
total_grains = max(array(:,7));
    
% write data into a text file
file = fopen('example_1_file_rhocp.txt','wt');
fprintf(file,'%s \n','# Header:    Marmot Input File');
fprintf(file,'%s \n','# Date:      11-April-2024 18:09:16');
fprintf(file,'%s \n','#');
fprintf(file,'%s \n','# Column 1:  Euler angle "phi1" ');
fprintf(file,'%s \n','# Column 2:  Euler angle "PHI" ');
fprintf(file,'%s \n','# Column 3:  Euler angle "phi2" ');
fprintf(file,'%s \n','# Column 4:  x-coordinate (in microns)');
fprintf(file,'%s \n','# Column 5:  y-coordinate (in microns)');
fprintf(file,'%s \n','# Column 6:  z-coordinate (in microns)');
fprintf(file,'%s \n','# Column 7:  grain number (integer)');
fprintf(file,'%s \n','# Column 8:  phase number (integer)');
fprintf(file,'%s \n','# Column 9:  Symmetry class (from TSL)');
fprintf(file,'%s \n','#');
fprintf(file,'%s %s\n','# Phase 1: ', phase_1_name);
fprintf(file,'%s %d \n','# Phase 1 Crystal Class: ', phase_1_sym);
fprintf(file,'%s %s\n','# Phase 2: ', phase_2_name);
fprintf(file,'%s %d \n','# Phase 2 Crystal Class: ', phase_2_sym);
fprintf(file,'%s %d \n','# Total Number of Grains: ', total_grains);
fprintf(file,'%s \n','#');
fprintf(file,'%s %f \n','# X_Min:  ',x_min);
fprintf(file,'%s %f \n','# X_Max:  ',x_max);
fprintf(file,'%s %f \n','# X_step:  ',x_step);
fprintf(file,'%s %d \n','# X_Dim:   ',x_dim);
fprintf(file,'%s \n','#');
fprintf(file,'%s %f \n','# Y_Min:  ',y_min);
fprintf(file,'%s %f \n','# Y_Max:  ',y_max);
fprintf(file,'%s %f \n','# Y_step:  ',y_step);
fprintf(file,'%s %d \n','# Y_Dim:   ',y_dim);
fprintf(file,'%s \n','#');
fprintf(file,'%s %f \n','# Z_Min:  ',z_min);
fprintf(file,'%s %f \n','# Z_Max:  ',z_max);
fprintf(file,'%s %f \n','# Z_step:  ',z_step);
fprintf(file,'%s %d \n','# Z_Dim:   ',z_dim);
fprintf(file,'%s \n','#');
fprintf(file,'%s \n','# Note:  This is a 3D dataset with one layer in Z!');
fprintf(file,'%s \n','#');
for ia = 1:size(array,1)
    fprintf(file,'%4.2f%s %4.2f%s %4.2f%s %4.2f%s %4.2f%s %4.2f%s %d%s %d%s %d\n',array(ia,1),' ',array(ia,2),' ', ... 
        array(ia,3),' ',array(ia,4),' ',array(ia,5),' ',array(ia,6),' ', ...
        array(ia,7),' ', array(ia,8),' ',array(ia,9));
end
fclose(file);




