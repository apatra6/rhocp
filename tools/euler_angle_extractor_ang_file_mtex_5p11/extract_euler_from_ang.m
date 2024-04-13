function extract_euler_from_ang()
%% Extracting euler angles from EBSD (.ang) format data 
    clc;clear;
    % the header end line should be checked by opening the ang file ...
    % in a text editor before running this script
    header_line_end = 126;
    imp_data = importdata('input_ang_file.ang',' ',header_line_end);
    data = imp_data.data;
    
    % convert to degrees
    phi1 = data(:,1) * 180/pi;
    Phi = data(:,2) * 180/pi;
    phi2 = data(:,3) * 180/pi;
    
    % write data into a text file
    file = fopen('full_data_euler.txt','wt');
    for ia = 1:size(data,1)
        fprintf(file,'%e%s%e%s%e\n',phi1(ia),' ', ...
            Phi(ia),' ',phi2(ia));
    end
    fclose(file);
end
