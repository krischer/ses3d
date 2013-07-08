%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- write stf to file ------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if (accept>0)
               
    fn=[output_path num2str(event_list(idx_eq)) '/ad_src_' int2str(idx)];
    fid=fopen(fn,'w');
    
    %- write header -------------------------------------------------------
    
    fprintf(fid,'-- adjoint source ------------------\n');
    fprintf(fid,'-- source coordinates (colat,lon,depth)\n');
    fprintf(fid,'%f %f %f\n',90-locations.stla,locations.stlo,1e3);
    fprintf(fid,'-- source time function (x, y, z) --\n');

    %- write adjoint source time function ---------------------------------
    
    for m=1:nt
        fprintf(fid,'%f %f %f\n',weight_total*ad_src_x(m),weight_total*ad_src_y(m),weight_total*ad_src_z(m));
    end
                
    fclose(fid);
        
    %- note adjoint source location ---------------------------------------

    adloc.x(idx)=90-locations.stla;
    adloc.y(idx)=locations.stlo;
    adloc.z(idx)=1e3;
    
    %- rotate adjoint source location if necessary ------------------------
    
    if (strcmp(rot_flag,'yes'))
        
       [adloc.x(idx) adloc.y(idx)]=rotate_vector(rot_axis,rot_angle,adloc.x(idx),adloc.y(idx),'silent');
        
    end
       
end