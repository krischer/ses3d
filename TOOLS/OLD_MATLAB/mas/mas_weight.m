
%------------------------------------------------------------------
%- read or generate subjective station weight ---------------------
%------------------------------------------------------------------    

if (strcmp(mode,'auto'))
    
    weight_sub=fscanf(fid_logold,'%g',1);
    fgets(fid_logold);
    
elseif (strcmp(mode,'manual'))

    weight_sub=input('weight: '); 
  
end

%------------------------------------------------------------------
%- read and generate geographic weights ---------------------------
%------------------------------------------------------------------

%- epicentral distance weight -------------------------------------
    
weight_geo=1-exp(-Delta^2/5^2);

%- geographic station weight from station file --------------------
   
fn_station=[path_data 'STATIONFILES/' name '.' comp];
   
if (exist(fn_station,'file'))
    
    fid_station=fopen(fn_station,'r');

    weight_response=fscanf(fid_station,'%g',1);
    
    fclose(fid_station);
    
    fprintf(1,'response weight: %g\n',weight_response);
    
else
    
    weight_response=0.2;
    
    fprintf(1,'no station file, response weight: %g\n',weight_response);
    
end

weight_geo=weight_geo*weight_response;

%------------------------------------------------------------------
%- make total weight ----------------------------------------------
%------------------------------------------------------------------

weight_total=weight_sub*weight_geo*weight_event;

%------------------------------------------------------------------
%- write output ---------------------------------------------------
%------------------------------------------------------------------

fprintf(1,'weights: subjective: %g, geographic: %g, total: %g\n',weight_sub,weight_geo,weight_total);
   