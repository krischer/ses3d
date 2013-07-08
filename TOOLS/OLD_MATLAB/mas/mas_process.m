%==============================================================
%- taper data and synthetics
%==============================================================
            
%- check whether data are acceptable --------------------------
            
if (strcmp(mode,'auto'))
                
    s=fscanf(fid_logold,'%c',4);        %- station name
    accept=fscanf(fid_logold,'%d',1);   %- station accepted or not
             
    if ((strcmp(s,name)==1) & (accept>0))
               
        fprintf(1,'accepted\n');
                    
    elseif (accept==0)
                    
        fprintf(1,'not accepted\n');
        fscanf(fid_logold,'%d',5);
        fgets(fid_logold);
                    
    else
                    
        fprintf(1,'name inconsistency!\n');
        break;
                    
    end

elseif (strcmp(mode,'manual'))

    accept=[];
    while (isempty(accept) | not(isfloat(accept)))
        accept=input('accept=1, accept&flip=2, reject=0: ');
    end
    
end
            
%- taper if acceptable ----------------------------------------
                        
if (accept>0)
                
    fprintf(fid_err,'%d\n',idx);
    fprintf(fid_err,'%c',name(1:4));
    fprintf(fid_err,'.%c\n',comp);
    fprintf(fid_err,'coordinates (lat,lon): %f %f\n',locations.stla,locations.stlo);
                
    fprintf(fid_log,'%c',name(1:4));
    fprintf(fid_log,' %d ',accept);
                
    if (accept==2)
                    
        data=-data;
                    
        if (strcmp(plot_mode,'plot'))
                    
            plot(t,synth/ms,'r','LineWidth',1);
            hold on
            plot(t,data/md,'k','LineWidth',1);
            title(['station: ' name ', \Delta=' num2str(Delta) ' deg, black: data, red: synthetics']);
            xlabel('t [s]')
            axis([0 max(t) -1.2 1.2])
            hold off
            fprintf(1,'data flipped!\n');
                        
        end
                    
    end
                
    %- read from logfile if mode=='auto' ----------------------
                
    if (strcmp(mode,'auto'))    
                    
        num_win=fscanf(fid_logold,'%d',1);
        fprintf(1,'number of windows: %d\n',num_win);
        fprintf(fid_log,' %d ',num_win);
        
        taper_time=0*t;
        taper_time_weighted=0*t;
        
        t_width=taper_width;

        %- generate weighted multi-window ---------------------
                    
        for idx_win=1:num_win
                    
            x_left=fscanf(fid_logold,'%g',1);
            fprintf(1,'left taper boundary: %g\n',x_left);
                    
            x_right=fscanf(fid_logold,'%g',1);
            fprintf(1,'right taper boundary: %g\n',x_right);
                       
            weight_win=fscanf(fid_logold,'%g',1);
            fprintf(1,'window weight: %g\n',weight_win);
                     
            taper_time=taper_time+taper(0*t+1,t,x_left,x_right,t_width);
            taper_time_weighted=taper_time_weighted+weight_win*taper(0*t+1,t,x_left,x_right,t_width);
                        
            fprintf(fid_log,'%g %g %g ',x_left, x_right, weight_win);
            fprintf(fid_err,'taper (left, right, width, weight)=(%f %f %f %f)\n',x_left,x_right,t_width,weight_win);
                        
        end
                                     
    elseif (strcmp(mode,'manual'))
                    
        num_win=[];
        while (isempty(num_win) | not(isfloat(num_win)) | (num_win>5))
            num_win=input('number of windows: ');
        end

        fprintf(fid_log,' %d ', num_win);
        
        taper_time=0*t;
        taper_time_weighted=0*t;
        
        t_width=taper_width;
                    
        %- generate weighted multi-window ---------------------
                    
        for idx_win=1:num_win
                        
            fprintf(1,'left taper boundary: ');
            x_left=[];
            while (isempty(x_left))
                [x_left dummy]=ginput(1);
            end
            fprintf(1,'%f s',x_left);
                    
            fprintf(1,'\nright taper boundary: ');
            x_right=[];
            while (isempty(x_right))
                [x_right dummy]=ginput(1);
            end
            fprintf(1,'%f s\n',x_right);
                        
            weight_win=[];
            while (isempty(weight_win))
                weight_win=input('window weight: ');
            end
                        
            taper_time=taper_time+taper(0*t+1,t,x_left,x_right,t_width);
            taper_time_weighted=taper_time_weighted+weight_win*taper(0*t+1,t,x_left,x_right,t_width);
                        
            fprintf(fid_log,'%g %g %g ',x_left, x_right, weight_win);
            fprintf(fid_err,'taper (left, right, width, weight)=(%f %f %f %f)\n',x_left,x_right,t_width,weight_win);
                        
        end
            
    end

    idx=idx+1;
    
    data=data.*taper_time_weighted;
    synth_tap=synth.*taper_time;
    synth_tap_weighted=synth.*taper_time_weighted;
            
    if (strcmp(plot_mode,'plot'))
        
        plot(t,data/md,'k');
        hold on
        plot(t,synth_tap_weighted/ms,'r');

        if (strcmp(old,'yes'))
            plot(t,synth_old(1:length(t)).*taper_time_weighted/ms,'r:')
        end
                    
        title('black: tapered data, red: tapered synthetics, dashed: tapered synthetic-data');
        xlabel('t [s]')
        grid on
        hold off
        
        pause(0.5)
                    
    end
    
end