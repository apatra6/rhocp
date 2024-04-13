function paraview_to_ang()
%% MATLAB script to generate .ang file from MOOSE Exodus output
    clc; clear;
    
    % original input microstructure
    fname = 'example_2_file.ang';
    ebsd_data = importdata(fname,' ',218);
    inp_data = ebsd_data.data;
    
    % csv file exported from Paraview
    imp_data = importdata('point_to_cell_data.csv');
    data = imp_data.data;
    ori = [data(:,11)  data(:,5)  data(:,12)] * pi/180;
    phase_id = data(:,3);
    
    array = [ori inp_data(:,4) abs( max(inp_data(:,5)) - inp_data(:,5) ) ... 
        ones(size(inp_data,1),1) ones(size(inp_data,1),1) phase_id ... 
        ones(size(inp_data,1),1) ones(size(inp_data,1),1)];
    
    
    [pathstr,name,ext] = fileparts(fname);
    outfile = fullfile(pathstr,['example_2_out',ext]);
    write_OIM_data(array,fname,outfile);
end

