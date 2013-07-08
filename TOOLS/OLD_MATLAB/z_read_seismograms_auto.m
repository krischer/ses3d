function p=z_read_seismograms_auto(directory,component,station,mode)

% function p=z_read_seismograms(array,spacing)
%
% reads seismograms out of the output files and plots them
%
% INPUT:
%
% array:    array of the indeces of the seismograms that will be plotted,
%           enter [] if all seismograms are to be plotted
%
% spacing:  spacing of seismograms in the plot as multiple of the maximum
%           amplitude. A reasonable value is 2.0 .
%
% mode:     'plot' for plotting, otherwise 'silent'
%
% OUTPUT:   a seismogram structure
%
% last modified: 19 October, 2009

filename=[directory station '_' component];
%disp(filename);
%disp(directory);

if (exist(filename,'file'))

    fid=fopen(filename,'r');

    fgets(fid);
    fscanf(fid,'%c',4);
    p.nt=fscanf(fid,'%d',1);
    fscanf(fid,'%c',5);
    p.dt=fscanf(fid,'%g',1);

    p.seismograms=zeros(1,p.nt);
    p.t0=0.0;

    fgets(fid);
    fgets(fid);
    fscanf(fid,'%c',4);
    p.r_x=fscanf(fid,'%g',1);
    fscanf(fid,'%s',1);
    p.r_y=fscanf(fid,'%g',1);
    fscanf(fid,'%s',1);
    p.r_z=fscanf(fid,'%g',1);

    fgets(fid);
    fgets(fid);
    fscanf(fid,'%s',1);
    p.s_x=fscanf(fid,'%g',1);
    fscanf(fid,'%s',1);
    p.s_y=fscanf(fid,'%g',1);
    fscanf(fid,'%s',1);
    p.s_z=fscanf(fid,'%g',1);
    
    p.s_x=90-p.s_x;
    p.r_x=90-p.r_x;
    Delta=epicentral_distance(p.s_x,p.s_y,p.r_x,p.r_y);

    p.seismograms(1:p.nt)=fscanf(fid,'%g',p.nt);
   
    fclose(fid);
    
else
    
    p.nt=1000;
    p.dt=0.2;
    p.seismograms(1:1000)=0.0;
    p.t0=0;
    p.r_x=0.0;
    p.r_y=0.0;
    p.r_z=0.0;
    p.s_x=0.0;
    p.s_y=0.0;
    p.s_z=0.0;
    
end

%==========================================================================
% plot seismograms
%==========================================================================

if (strcmp(mode,'plot'))

    m=max(abs(p.seismograms));
    t=0:p.dt:((p.nt-1)*p.dt);

    figure
    hold on

    plot(t,p.seismograms,'k');
    axis([t(1)-0.1*(t(end)-t(1)) t(end)+0.1*(t(end)-t(1)) -1.1*m 1.1*m]);

    text(t(1), m, [station ' (\Delta=' num2str(Delta) '°)']);
    text(t(1),0.9*m, ['source: lat=' num2str(p.s_x) '°, lon=' num2str(p.s_y) '°, depth=' num2str(p.s_z/1000) 'km']);
    text(t(1),0.8*m, ['receiver: lat=' num2str(p.r_x) '°, lon=' num2str(p.r_y) '°, depth=' num2str(p.r_z/1000) 'km']);
    
    xlabel('time [s]');

    if (strcmp(component,'x'))

        ylabel('\theta displacement [m]');

    elseif (strcmp(component,'y'))

        ylabel('\phi displacement [m]');

    elseif (strcmp(component,'z'))

        ylabel('z displacement [m]');

    end

end