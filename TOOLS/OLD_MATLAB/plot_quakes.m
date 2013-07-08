%--------------------------------------------------------------------------
% open earthquake catalogue
%--------------------------------------------------------------------------

    
if ((value>=0) && (value<30))
       
    fid=fopen('QUAKES/quakes_europe_0_30','r');
    opened=1;
       
elseif ((value>=30) && (value<60))
       
    fid=fopen('QUAKES/quakes_europe_30_60','r');
    opened=1;
    
elseif ((value>=60) && (value<100))
       
    fid=fopen('QUAKES/quakes_europe_60_100','r');
    opened=1;
    
elseif ((value>=100) && (value<150))
       
    fid=fopen('QUAKES/quakes_europe_100_150','r');
    opened=1;
       
elseif ((value>=150) && (value<200))
       
    fid=fopen('QUAKES/quakes_europe_150_200','r');
    opened=1;
    
elseif ((value>=200) && (value<300))
       
    fid=fopen('QUAKES/quakes_europe_200_300','r');
    opened=1;
    
elseif ((value>=300) && (value<400))
       
    fid=fopen('QUAKES/quakes_europe_300_400','r');
    opened=1;    
    
elseif ((value>=400) && (value<600))
       
    fid=fopen('QUAKES/quakes_europe_400_600','r');
    opened=1;    
       
end

%--------------------------------------------------------------------------
% plot earthquakes
%--------------------------------------------------------------------------

while (feof(fid)==0)

    fscanf(fid,'%s',1);
    fscanf(fid,'%d',3);
    fscanf(fid,'%g',1);
    lat_quake=fscanf(fid,'%g',1);
    lon_quake=fscanf(fid,'%g',1);
    m_plot(lon_quake,lat_quake,'k.','MarkerSize',5);
    fgetl(fid);
    
end
