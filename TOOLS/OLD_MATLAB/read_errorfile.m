
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads the error-file and plots the results
%
% last modfied: 29 Januar, 2009 by Andreas Fichtner
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'/home/fichtner/mmap/m_map');

path_out='/nfs/stig/fichtner/SES3D_3.MER.T1/ADJOINT/42.tf/'

event_list=[1002 1005 1009 1010 1012 1026 1027 1029 1030 1014 1017 1018 1019 1020 1021 1035 1037 1039 1045 1048 1051 1054 1062 1064 1065 1066 1067 1068 1069 1071 1072 1073 1074 1075 1076];
event_list=[1002 1005 1009 1010 1012 1027 1029 1030 1014 1017 1018 1019 1020 1021 1037 1039 1048 1054 1062 1065 1066 1068 1069 1071 1072 1073 1074 1075 1076];

cmap='mono';        % rgb, gray, mono

bin_limits=0:0.05:3.5;
%bin_limits=-10:0.5:10;

bin=zeros(1,length(bin_limits));

plot_station_position=1;
plot_station_name=1;
plot_ray_paths=0;

%upper_bound=3.5;        %- upper misfit bound 
upper_bound=30.5;        %- upper misfit bound 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plot_ray_paths)

    figure
    title('Misfits and ray paths');
    m_proj('lambert','lon',[22 46],'lat',[33 44]);
    m_coast('patch',[0 1 0],'edgecolor','k','FaceAlpha',0.1,'LineWidth',2);
    m_grid('linest','-','xticklabels',[],'yticklabels',[])
    hold on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over events, collect misfits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=0;

E_total=0;
E_mean=0;

window_length=0;
window_length_total=0;
window_max=0;
n_window_total=0;
n_window_event=0;

E_min=inf;
E_max=-inf;

num_events=length(event_list);

E_ev=zeros(1,num_events);

E_vector=zeros(1,400);

N_old=0;

nnz=0;

for ev=1:num_events
    
    %- open error file and read event information -------------------------
    
    fn=[path_out num2str(event_list(ev)) '/errorfile'];
    fprintf(1,'reading errorfile %s\n',fn);
    fid=fopen(fn,'r');
    
    fscanf(fid,'%c',29);
    evla=fscanf(fid,'%g',1);
    evlo=fscanf(fid,'%g',1);
    fgetl(fid);
    
    n_window_event=0;
    E_mean_event=0;
    n_event=0;
    
    %- loop through individual stations -----------------------------------
    
    while (feof(fid)==0)
    
        N=N+1;
        n_event=n_event+1;
        
        %- read station information ---------------------------------------

        fgetl(fid);
        name=fscanf(fid,'%c',4);
        fscanf(fid,'%c',1);
        comp=fscanf(fid,'%c',1);
        fscanf(fid,'%c',24);
        stla=fscanf(fid,'%g',1);
        stlo=fscanf(fid,'%g',1);
        fgetl(fid);
        
        epi=epicentral_distance(evla,evlo,stla,stlo);
        
        %- read misfits ---------------------------------------------------
        
        window_length=0;
        
        while (strcmp(fscanf(fid,'%c',5),'taper'));
            
            fscanf(fid,'%c',31);
            taper_left=fscanf(fid,'%g',1);
            taper_right=fscanf(fid,'%g',1);
            fgetl(fid);
            window_length=window_length+(taper_right-taper_left);
            n_window_event=n_window_event+1;
            n_window_total=n_window_total+1;
            
            if (taper_right>window_max)
                window_max=taper_right;
            end
            
        end
                
        window_length_total=window_length_total+window_length;
        
        fscanf(fid,'%c',9);
        E=fscanf(fid,'%g',1);
        w_sub=fscanf(fid,'%g',1);
        w_geo=fscanf(fid,'%g',1);
        w_event=fscanf(fid,'%g',1);
        fgetl(fid);

        E=E*exp(-E^4/upper_bound^4);
        E=w_event*w_sub*w_geo*E;
        
        %- fill bins ------------------------------------------------------
        
        for i=length(bin_limits):-1:1
            if (E>bin_limits(i))
                bin(i)=bin(i)+1;
                break;
            end
        end
        
        %- sum misfits ----------------------------------------------------
        
        E_ev(ev)=E_ev(ev)+abs(E);

        E_total=E_total+abs(E);
        E_mean=E_mean+E;
        E_mean_event=E_mean_event+E;
        
        fscanf(fid,'%c',1);

        %- minima, maxima, total ------------------------------------------
       
        E_min=min(E_min,E);
        E_max=max(E_max,E);
        
    end
    
    E_ev(ev)=E_ev(ev)/n_window_event;
    
    fprintf(1,'windows: %d, misfit per window: %g, mean misfit: %g\n',n_window_event, E_ev(ev),E_mean_event/n_event);
    fclose(fid);
    
end

E_max=1.1*max(E_max,abs(E_min));
E_max=max(E_max,upper_bound);
E_min=-E_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over events, generate map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if (plot_ray_paths)
if (1==1)

for ev=1:num_events
    
    %- open error file and read event information -------------------------
    
    fn=[path_out num2str(event_list(ev)) '/errorfile'];
    fid=fopen(fn,'r');
    
    fscanf(fid,'%c',29);
    evla=fscanf(fid,'%g',1);
    evlo=fscanf(fid,'%g',1);
    fgetl(fid);
    
    %m_plot(evlo,evla,'bx','MarkerSize',20,'LineWidth',3);
    %m_text(1.01*evlo,1.01*evla,num2str(event_list(ev)),'Color','k','FontSize',15,'FontWeight','bold');
    plot(evlo,evla,'ro','MarkerSize',10,'LineWidth',3);

    %- loop through individual stations -----------------------------------
    
    while (feof(fid)==0)
        
        %- read station information ---------------------------------------

        fgetl(fid);
        name=fscanf(fid,'%c',4);
        fscanf(fid,'%c',1);
        comp=fscanf(fid,'%c',1);
        fscanf(fid,'%c',24);
        stla=fscanf(fid,'%g',1);
        stlo=fscanf(fid,'%g',1);
        fgetl(fid);
        
        epi=epicentral_distance(evla,evlo,stla,stlo);
        
        %- read misfits ---------------------------------------------------
        
        window_length=0;
        
        while (strcmp(fscanf(fid,'%c',5),'taper'))
            fscanf(fid,'%c',31);
            taper_left=fscanf(fid,'%g',1);
            taper_right=fscanf(fid,'%g',1);
            fgetl(fid);
            window_length=window_length+(taper_right-taper_left);
        end
        
        window_length_total=window_length_total+window_length;
        
        fscanf(fid,'%c',9);
       
        E=fscanf(fid,'%g',1);
        w_sub=fscanf(fid,'%g',1);
        w_geo=fscanf(fid,'%g',1);
        w_event=fscanf(fid,'%g',1);
        
        E=E*exp(-E^4/upper_bound^4);
        E=w_event*w_sub*w_geo*E;
        
        fgetl(fid);

        %- generate ray path + error --------------------------------------
     
        if (strcmp(comp,'x') || strcmp(comp,'y') || strcmp(comp,'z'))
        
            %- plot station position and station names --------------------
            
            if (plot_station_position==1)
                %m_plot(stlo,stla,'r.','MarkerSize',20,'LineWidth',2);
                plot(stlo,stla,'k.','MarkerSize',20,'LineWidth',2);
            end
            
            if (plot_station_name==1)
                %m_text(1.01*stlo,1.01*stla,name,'Color','k','FontSize',10,'FontWeight','bold');
                text(1.001*stlo,1.001*stla,name,'Color','k','FontSize',12,'FontWeight','bold');
            end
            
            %- compute ray path -------------------------------------------

        	[range,ln,lt]=m_lldist([stlo evlo],[stla evla],50);
                
            for k=1:length(ln)
                if (ln(k)<0)
                    ln(k)=360+ln(k);
                end
            end
            
            %- plot ray path ----------------------------------------------

            if (strcmp(cmap,'rgb'))
                if (E>0)
                    red=1;          green=1-E/E_max;    blue=1-E/E_max;
                elseif (E<0)
                    red=1+E/E_max;  green=1+E/E_max;    blue=1;
                elseif (E==0)
                    red=1;          green=1;            blue=1;
                end
            end
        
            th=0.85;
            if (strcmp(cmap,'gray'))
                if (E>th)
                    red=0;  blue=0; green=0;
                else
                    red=1-E/th;    blue=1-E/th;   green=1-E/th;
                end
            end
        
            if (strcmp(cmap,'mono'))
                    red=0;  blue=0; green=0;
            end
           
            if (abs(E)<=upper_bound)
                %m_line(ln,lt,'color',[red green blue],'LineWidth',1);
                %m_line(ln(5:45),lt(5:45),'color',[red green blue],'LineWidth',2);
                %m_line(ln(24:26),lt(24:26),'color',[red green blue],'LineWidth',2);
            end
            
        end
 
    end
        
    fclose(fid);
  
    pause(0.1)
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'misfit (total, mean, min, max): %g, %g, %g, %g\n',E_total,E_mean/N,E_min,E_max);
fprintf(1,'number of recordings: %d\n',N);
fprintf(1,'number of windows: %d\n',n_window_total);
fprintf(1,'total window length: %g hours\n',window_length_total/(60*60));

fprintf(1,'window max=%g\n',window_max);

figure
bar(bin_limits,bin,'k');
title('misfit distribution');
