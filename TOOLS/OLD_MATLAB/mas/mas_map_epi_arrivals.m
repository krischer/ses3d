            %==============================================================
            % generate map
            %==============================================================

            dt=t(2)-t(1);
            nt=length(t);
      
            if (strcmp(plot_mode,'plot'))
                
                %tic
                %h_map=figure('Color',[1 1 1],'Position',[10 40 S_height/2.5 S_height/3]);
                %h_seis=figure('Position',[10 S_height-10 0.98*S_width S_height/2.2]);
                %toc
                
                %- make a dot in the map --------------------------------------
	
                figure(h_map)
                
                m_plot(locations.stlo,locations.stla,'xg');
                m_plot(locations.evlo,locations.evla,'or');
                
            	[range,ln,lt]=m_lldist([locations.stlo locations.evlo],[locations.stla locations.evla],10);
                
                for k=1:length(ln)
                    if (ln(k)<0)
                        ln(k)=360+ln(k);
                    end
                end
                
                if (exist('ln_old'))
                    m_line(ln_old,lt_old,'color','r','LineWidth',3);
                end
                
            	m_line(ln,lt,'color','b','LineWidth',3);
 
                ln_old=ln;
                lt_old=lt;
                
            end
	
            %==============================================================
            %- compute epicentral distance and output information
            %==============================================================
            
            Delta=epicentral_distance(locations.evla,locations.evlo,locations.stla,locations.stlo);
            fprintf(1,'station: %s: lat=%g °, lon=%g °, delta=%f °\n',name,locations.stla,locations.stlo,Delta);
            
            %==============================================================
            %- compute arrival times of major phases
            %==============================================================
            
            if (strcmp(plot_mode,'plot'))
            
                tS=interp1(epi,S,Delta);
                tSS=interp1(epi,SS,Delta);
                tP=interp1(epi,P,Delta);
                tPP=interp1(epi,PP,Delta);
                tPcP=interp1(epi,PcP,Delta);
                tScS=interp1(epi,ScS,Delta);
                tScP=interp1(epi,ScP,Delta);
                
            end