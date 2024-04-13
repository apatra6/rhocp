function write_inp_cubic_grains()

% This code generate a mesh file in the ABAQUS .inp file format with cube-shaped grains
% The .inp file can be read by MOOSE/rho-CP

    % input file name
    file = fopen(strcat(int2str(total_grains),'grains_', ...
        int2str(total_elems),'elements.inp'),'wt');

    %%  Define the size of the model in terms of unit cells and grain size
    n_x_grain = 4; % no. of grains along x-direction
    n_y_grain = 4; % no. of grains along y-direction
    n_z_grain = 4; % no. of grains along z-direction
    grain_size = 0.1; % (mm)
    grn_xdir = 0.1; % (mm)
    grn_ydir = 0.1; % (mm)
    grn_zdir = 0.1; % (mm)
    %  Define the number of elements of each grain
    n_x_ele = 2; % no. of elements along x-direction in each grain
    n_y_ele = 2; % no. of elements along y-direction in each grain
    n_z_ele = 2; % no. of elements along z-direction in each grain
	
    el_size = grain_size/n_x_ele;
    total_grains = n_x_grain * n_y_grain * n_z_grain
    total_elems = n_x_grain * n_x_ele * n_y_grain * n_y_ele * ...
        n_z_grain * n_z_ele
    
    %  Define the size of model
    x_size = n_x_grain*grn_xdir;
    y_size = n_y_grain*grn_ydir;
    z_size = n_z_grain*grn_zdir;
    %  Define the local origin
    x0 = 0.0;
    y0 = 0.0;
    %  Calculate Loop variables
    n_node_x = n_x_grain * n_x_ele + 1;
    n_node_y = n_y_grain * n_y_ele + 1;
    n_node_z = n_z_grain * n_z_ele + 1;
    n_total_grain=n_x_grain*n_y_grain*n_z_grain;
    
    %% Writing into ABAQUS input file
    
    fprintf(file,'%s\n','*NODE');
    nctr = 0;
    for iz = 1:n_node_z
        for iy = 1:n_node_y
            for ix = 1:n_node_x
                nctr = ix + (iy - 1)*1000 + (iz -1)*1000000;
                fprintf(file,'%d%s%4.2f%s%4.2f%s%4.2f\n',nctr,',',(ix - 1)*el_size,',',(iy - 1)*el_size,',',(iz - 1)*el_size);
            end
        end
    end
    
    %  Generate vertex node sets
    fprintf(file,'%s\n','** Vertex node sets:8');
    %    generate vertex (0,0,0)
    fprintf(file,'%s\n','*NSET, nset=V000');
    fprintf(file,'%s\n','1');
    %    generate vertex (1,0,0)
    fprintf(file,'%s\n','*NSET, nset=V100');
    fprintf(file,'%6d\n',n_node_x);
    %    generate vertex (0,1,0)
    fprintf(file,'%s\n','*NSET, nset=V010');
    fprintf(file,'%6d\n',(n_node_y-1)*1000+1);
    %    generate vertex (1,1,0)
    fprintf(file,'%s\n','*NSET, nset=V110');
    fprintf(file,'%6d\n',(n_node_y-1)*1000+n_node_x);
    %    generate vertex (0,0,1)
    fprintf(file,'%s\n','*NSET, nset=V001');
    fprintf(file,'%6d\n',(n_node_z-1)*1000000+1);
    %    generate vertex (1,0,1)
    fprintf(file,'%s\n','*NSET, nset=V101');
    fprintf(file,'%6d\n',(n_node_z-1)*1000000+n_node_x);
    %    generate vertex (0,1,1)
    fprintf(file,'%s\n','*NSET, nset=V011');
    fprintf(file,'%6d\n',(n_node_z-1)*1000000+(n_node_y-1)*1000+1);
    %    generate vertex (1,1,1)
    fprintf(file,'%s\n','*NSET, nset=V111');
    fprintf(file,'%6d\n',(n_node_z-1)*1000000+(n_node_y-1)*1000+n_node_x);
    %  generate point (1/2,1/2,0)
    fprintf(file,'%s\n','*NSET, nset=Vh0');
    fprintf(file,'%6d\n',floor(1+((n_node_y-1)*1000+n_node_x)/2));
    
    %  Generate all elements
    fprintf(file,'%s\n','*ELEMENT, type=C3D8R');
    
    for iz = 1:n_z_ele*n_z_grain
        for iy = 1:n_y_ele*n_y_grain
            for ix = 1:n_x_ele*n_x_grain
                el_no = ix + 1000*(iy - 1) + 1000000*(iz - 1);
                node1 = ix + 1000*(iy - 1) + 1000000*(iz - 1);
                node2 = ix + 1 + 1000*(iy - 1) + 1000000*(iz - 1);
                node3 = ix + 1 + 1000*iy + 1000000*(iz - 1);
                node4 = ix + 1000*iy + 1000000*(iz - 1);
                node5 = ix + 1000*(iy - 1) + 1000000*iz;
                node6 = ix + 1 + 1000*(iy - 1) + 1000000*iz;
                node7 = ix + 1 + 1000*iy + 1000000*iz;
                node8 = ix + 1000*iy + 1000000*iz;
    
                fprintf(file,'%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d\n',el_no,',',node1,',',node2,',',node3,',',node4,',',node5,',',node6,',',node7,',',node8);
            end
        end
    end
    
    
    %  Generate element sets for each grain
    for i=1:n_z_grain
        for j=1:n_y_grain
            for k=1:n_x_grain
                n_grain=(i-1)*n_y_grain*n_x_grain+(j-1)*n_x_grain+k;
                fprintf(file,'%s%d\n','*ELSET, ELSET=GRAIN',n_grain);
                for i1=1:n_z_ele
                for j1=1:n_y_ele
                    fprintf(file,'%d%s%d%s',((i-1)*n_z_ele+i1-1)*1000000+((j-1)*n_y_ele+j1-1)*1000+((k-1)*n_x_ele)+1,',',((i-1)*n_z_ele+i1-1)*1000000+((j-1)*n_y_ele+j1-1)*1000+((k-1)*n_x_ele)+n_x_ele,',');
                end
                end
                fprintf(file,'\n');
            end
        end
    end
	
fclose(file);
