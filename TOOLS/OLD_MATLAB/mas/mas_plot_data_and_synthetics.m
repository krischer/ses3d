            %==============================================================
            %- plot data and synthetics
            %==============================================================
         
            if (strcmp(plot_mode,'plot'))
            
                %- plot data and synthetics -----------------------------------
            
                figure(h_seis);

                md=max(abs(data));
                ms=max(abs(synth));
                
                if (strcmp(old,'yes'))
                    plot(t_old,synth_old/ms,'r:','LineWidth',3);
                    hold on
                end
                
                plot(t,synth/ms,'r','LineWidth',2);
                hold on
                plot(t,data/md,'k','LineWidth',2);
                
                axis([0 max(t) -1.1 1.1]);
                
                %- plot arrival times of principal phases ---------------------
            
                if (1==0)
                
                    scal=1;
                    
                    plot([tS tS],[-0.3*scal -0.2*scal],'k');
                    text(tS,-0.4*scal,'S');
                    plot([tSS tSS],[-0.4*scal -0.5*scal],'k');
                    text(tSS,-0.6*scal,'SS');
                    plot([tP tP],[-0.3*scal -0.2*scal],'k');
                    text(tP,-0.4*scal,'P');
                    plot([tPP tPP],[-0.4*scal -0.5*scal],'k');
                    text(tPP,-0.6*scal,'PP');
                    plot([tScS tScS],[-0.3*scal -0.2*scal],'k');
                    text(tScS,-0.4*scal,'ScS');
                    plot([tPcP tPcP],[-0.4*scal -0.5*scal],'k');
                    text(tPcP,-0.6*scal,'PcP');
                    plot([tScP tScP],[-0.3*scal -0.2*scal],'k');
                    text(tScP,-0.4*scal,'ScP');
                
                    plot([Delta*111/5 Delta*111/5],[-0.15*scal 0.15*scal],'b');
                    text(Delta*111/5,0.2*scal,'5 km/s');

                    plot([Delta*111/4 Delta*111/4],[-0.15*scal 0.15*scal],'b');
                    text(Delta*111/4,0.3*scal,'4 km/s');
                
                    plot([Delta*111/3.5 Delta*111/3.5],[-0.15*scal 0.15*scal],'b');
                    text(Delta*111/3.5,0.2*scal,'3.5 km/s');
                
                    plot([Delta*111/3 Delta*111/3],[-0.15*scal 0.15*scal],'b');
                    text(Delta*111/3,0.3*scal,'3 km/s');

                end
            
                title(['station: ' name '.' comp ', \Delta=' num2str(Delta) ' deg, black: data, red: synthetics'],'FontSize',14);
                xlabel('t [s]')
                %grid on
                hold off
        
                pause(1.5)
                
            end
        