function mas(mode,old,plot_mode)

%==========================================================================
% generate adjoint sources
% last modified 29 March, 2010. Andreas Fichtner
%==========================================================================

% mode = ('auto','manual')                      % make everything manually or automatically based on old logfile
% old = ('yes','no')                            % plot data from previous iteration
% plot_mode = ('silent','plot')                 % make figures ('plot') or not ('silent')

%- input ------------------------------------------------------------------

mas_input

num_events=length(event_list);
        
%==========================================================================
% loop over earthquakes
%==========================================================================

for idx_eq=1:num_events
    
    %- initialisation -----------------------------------------------------
    
    mas_init;
    
    %- loop over recordings -----------------------------------------------
    
    idx=0;
    for irec=1:nrec
        
        mas_read_data_and_synthetics
        mas_map_epi_arrivals
        mas_plot_data_and_synthetics
        mas_process

        
        %- make adjoint source time function ------------------------------
                
        if (accept>0)
            
            mas_weight
            
            if (strcmp(misfit,'cc'))
                adsrc_cc  
            elseif (strcmp(misfit,'tf'))
                adsrc_tf
            elseif (strcmp(misfit,'amp'))
                adsrc_amp
            end
            
        end
            
        %- write adjoint source time function to file ---------------------
        
        mas_write_adstf

        %- write logfile and close figures --------------------------------
        
        if (accept>0)
            fprintf(fid_log,'%g\n',weight_sub);     % subjective weight
        end
            
        if (accept==0)
            fprintf(fid_log,'%c',name(1:4));
            fprintf(fid_log,' 0 0 0 0 0 0\n');    
        end
            
    end     % loop over recordings

    %- write adjoint source locations -------------------------------------
    
    fprintf(fid_src,'%d\n',idx);
    
    for i=1:idx
        
        fprintf(fid_src,'%f %f %f\n',adloc.x(i),adloc.y(i),adloc.z(i));
        
    end
    
    %- clean up -----------------------------------------------------------
    
    fclose(fid_src);
    fclose(fid_err);
    fclose(fid_log);
        
    if (strcmp(mode,'auto') & strcmp(old,'yes'))    
        fclose(fid_logold);
    end
        
end             % loop over earthquakes
