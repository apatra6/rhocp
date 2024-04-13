function write_OIM_data(OIM_data,header,outfile)
fid_outfile = fopen(outfile,'wb');
fid_header = fopen(header,'rb');
        while 1
            tline = fgetl(fid_header);
            if ~ischar(tline), break, end
            
            if(tline(1)~='#')
%                 disp('Running');
                break;
            end
            fprintf(fid_outfile,'%s \r\n',tline);
%             disp(tline)
        end
fclose(fid_header);
fprintf(fid_outfile,'%4.2f   %4.2f   %4.2f   %4.2f   %4.2f   %4.2f   %4.2f   %d   %4.2f   %4.2f \r\n',OIM_data');
fclose(fid_outfile);
