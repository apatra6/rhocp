function random_ori_generator()
%% Generate random orientations for n_total_grains
    n_total_grains = 512;
    ori = orientation.rand(n_total_grains);
    
    % convert to degrees
    phi1_val = ori.phi1 * 180/pi; 
    phi_val = ori.Phi * 180/pi; 
    phi2_val = ori.phi2 * 180/pi;

    % write data into a text file
    file = fopen('orientations.in','wt');
    fprintf(file,'%d\n',n_total_grains);
    for ia = 1:n_total_grains
        fprintf(file,'%e%s%e%s%e\n',phi1_val(ia),' ', ...
            phi_val(ia),' ',phi2_val(ia));
    end
    fclose(file);
end