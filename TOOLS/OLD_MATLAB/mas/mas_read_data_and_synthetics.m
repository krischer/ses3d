%==================================================================
% read data and synthetics
%==================================================================
                
fprintf(1,'\n*******************************************************\n');
            
[t data synth locations name comp]=read_data_and_synthetics(path_data,path_synth,[num2str(event_list(idx_eq)) period_tag],stations(irec,:),components(irec,1:3));
    
if (strcmp(old,'yes'))
    [t_old data_old synth_old locations name comp]=read_data_and_synthetics(path_data,path_synth_old,[num2str(event_list(idx_eq)) period_tag],stations(irec,:),components(irec,1:3));
end
                       
if (irec==1)
    fprintf(fid_err,'source coordinates (lat,lon): %g %g\n',locations.evla,locations.evlo);
end

%- apply station correction -----------------------------------------------

if (station_correction==1000)
    
   %- open station file ---------------------------------------------------
   
   fn_station=[path_data '/STATIONFILES/' name '.' comp];
   
   if (exist(fn_station,'file'))
   
        fid_station=fopen(fn_station,'r');
   
        fgetl(fid_station);
        fgetl(fid_station);
        fgetl(fid_station);

        ev=0;
    
        while (isempty(ev)==0)
        
            ev=fscanf(fid_station,'%d',1);
            dummy=fscanf(fid_station,'%g',1);
     
        end
   
        fgetl(fid_station);
        fgetl(fid_station);
        fgetl(fid_station);
   
        %- apply amplitude correction -------------------------------------
   
        amplitude_correction=fscanf(fid_station,'%g');
        fprintf(1,'apply amplitude correction: %g\n',amplitude_correction);
   
        synth=amplitude_correction*synth;
   
        if (strcmp(old,'yes'))
            synth_old=amplitude_correction*synth_old;
        end
   
        %- clean up -------------------------------------------------------
   
        fclose(fid_station);
        
   else
      
       fprintf(1,'no station file available\n');
       
   end
   
end