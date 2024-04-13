function neper_to_MOOSE()
%% Convert neper generated .inp file to .inp file that can be read by MOOSE/rho-CP
    clc; clear;
    filename = 'cube2_meshed.inp';
    id = readlines(filename);
    count = 1;
    for ia = 3:size(id,1)-2
        clear line_data
        line_data = id(ia);
        if (line_data=="")
            continue;
        end 
        new_id(count,:) = line_data;
        count = count + 1;
    end
    
    writelines(new_id,'cube2_meshed_MOOSE.inp');
end