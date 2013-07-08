function p=read_seismograms(directory,station,rot_mode,plot_mode)

% function p=read_seismograms(directory,station,mode)
%
% reads seismograms out of the output files and plots them
%
% INPUT:
%
% array:    directory where the SES3D seismograms are located
%
% station:  station name
%
% mode:     'plot' for plotting, otherwise 'silent'
%
% OUTPUT:   a seismogram structure
%
% last modified: 29 March, 2010

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- set paths --------------------------------------------------------------

path(path,'/home/fichtner/Matlab/Rotation/');

%- rotation axis unit vector and rotation angle ---------------------------

n=[0 1 0]';
phi=57.5;

%- filenames --------------------------------------------------------------

filename_x=[directory station '_x'];
filename_y=[directory station '_y'];
filename_z=[directory station '_z'];

%- read synthetic data ----------------------------------------------------

if (exist(filename_x,'file') && exist(filename_y,'file') && exist(filename_z,'file'))

    fid_x=fopen(filename_x,'r');
    fid_y=fopen(filename_y,'r');
    fid_z=fopen(filename_z,'r');

    fgets(fid_x);
    fscanf(fid_x,'%c',4);
    p.nt=fscanf(fid_x,'%d',1);
    fscanf(fid_x,'%c',5);
    p.dt=fscanf(fid_x,'%g',1);
    
    p.seismograms_x=zeros(1,p.nt);
    p.seismograms_y=zeros(1,p.nt);
    p.seismograms_z=zeros(1,p.nt);
    p.t0=0.0;

    fgets(fid_x);
    fgets(fid_x);
    
    tmp=regexp(fgetl(fid_x),'[xyz]=','split');
    p.r_x=sscanf(tmp{2},'%g');
    p.r_y=sscanf(tmp{3},'%g');
    p.r_z=sscanf(tmp{4},'%g');
    
    fgets(fid_x);
    tmp=regexp(fgetl(fid_x),'[xyz]=','split');
    p.s_x=sscanf(tmp{2},'%g');
    p.s_y=sscanf(tmp{3},'%g');
    p.s_z=sscanf(tmp{4},'%g');
    
    p.s_x=90-p.s_x;
    p.r_x=90-p.r_x;
    Delta=epicentral_distance(p.s_x,p.s_y,p.r_x,p.r_y);

    for i=1:7
        fgetl(fid_y);
        fgetl(fid_z);
    end
    
    p.seismograms_x(1:p.nt)=fscanf(fid_x,'%g',p.nt);
    p.seismograms_y(1:p.nt)=fscanf(fid_y,'%g',p.nt);
    p.seismograms_z(1:p.nt)=fscanf(fid_z,'%g',p.nt);
    
    fclose(fid_x);
    fclose(fid_y);
    fclose(fid_z);
    
else
    
    p.nt=1000;
    p.dt=0.2;
    p.seismograms_x(1:1000)=0.0;
    p.seismograms_y(1:1000)=0.0;
    p.seismograms_z(1:1000)=0.0;
    p.t0=0;
    p.r_x=0.0;
    p.r_y=0.0;
    p.r_z=0.0;
    p.s_x=0.0;
    p.s_y=0.0;
    p.s_z=0.0;
    
end

%==========================================================================
% rotate seismograms (preparatory steps)
%==========================================================================

if (strcmp(rot_mode,'yes'))

%- compute source and receiver coordinates in original system -------------

colat_new=90-p.r_x;    %- colatitude in rotated system
lon_new=p.r_y;         %- longitude in rotated system

[colat lon]=rotate_vector(n,-phi,90-p.r_x,p.r_y,'silent');
p.r_x=90-colat;         %- latitude in original system
p.r_y=lon;              %- longitude in original system

[p.s_x p.s_y]=rotate_vector(n,-phi,90-p.s_x,p.s_y,'silent');
p.s_x=90-p.s_x;         %- colatitude of source in original system

%- transform to radians ---------------------------------------------------

phi=phi*pi/180;

colat=colat*pi/180;
lon=lon*pi/180;

colat_new=colat_new*pi/180;
lon_new=lon_new*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotate seismograms (main calculations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=Rot(n,phi);   % rotation matrix that transforms the original system into the transformed system

%- compute e_theta, e_phi, e_r in the original coordinate system ----------
%- These are the local basis vectors in the original coordinate system ----

e1=[cos(lon)*cos(colat); sin(lon)*cos(colat); -sin(colat)];     % e_theta
e2=[-sin(lon); cos(lon); 0];                                    % e_phi
e3=[cos(lon)*sin(colat); sin(lon)*sin(colat); cos(colat)];      % e_r

%- compute e_theta, e_phi, e_r in the transformed coordinate system -------
%- These are the local basis vectors in the rotated coordinate system -----

e1_new=[cos(lon_new)*cos(colat_new); sin(lon_new)*cos(colat_new); -sin(colat_new)];         % e_theta
e2_new=[-sin(lon_new); cos(lon_new); 0];                                                    % e_phi
e3_new=[cos(lon_new)*sin(colat_new); sin(lon_new)*sin(colat_new); cos(colat_new)];          % e_r

%- express the rotated local basis vectors in terms of the original unit vector basis

e1_new=R'*e1_new;
e2_new=R'*e2_new;
e3_new=R'*e3_new;

%- compute transformation matrix ------------------------------------------

A=zeros(3);

A(1,1)=e1_new'*e1;
A(1,2)=e2_new'*e1;
A(1,3)=e3_new'*e1;

A(2,1)=e1_new'*e2;
A(2,2)=e2_new'*e2;
A(2,3)=e3_new'*e2;

A(3,1)=e1_new'*e3;
A(3,2)=e2_new'*e3;
A(3,3)=e3_new'*e3;

%- rotate seismogram ------------------------------------------------------

p_new=p;

p_new.seismograms_x=A(1,1)*p.seismograms_x+A(1,2)*p.seismograms_y+A(1,3)*p.seismograms_z;
p_new.seismograms_y=A(2,1)*p.seismograms_x+A(2,2)*p.seismograms_y+A(2,3)*p.seismograms_z;
p_new.seismograms_z=A(3,1)*p.seismograms_x+A(3,2)*p.seismograms_y+A(3,3)*p.seismograms_z;

p=p_new;

end

%==========================================================================
% plot seismograms
%==========================================================================

if (strcmp(plot_mode,'plot'))

    m=max(abs([p.seismograms_x p.seismograms_y p.seismograms_z]));
    t=0:p.dt:((p.nt-1)*p.dt);

    figure
    hold on

    subplot(3,1,1)
    plot(t,p.seismograms_x,'k');
    axis([t(1)-0.1*(t(end)-t(1)) t(end)+0.1*(t(end)-t(1)) -1.1*m 1.1*m]);
    text(t(1),0.9*m, [station ' (\Delta=' num2str(Delta) '°)']);
    text(t(1),0.7*m, ['source: lat=' num2str(p.s_x) '°, lon=' num2str(p.s_y) '°, depth=' num2str(p.s_z/1000) 'km']);
    text(t(1),0.5*m, ['receiver: lat=' num2str(p.r_x) '°, lon=' num2str(p.r_y) '°, depth=' num2str(p.r_z/1000) 'km']);
    xlabel('time [s]');
    ylabel('\theta velocity [m/s]');
    
    subplot(3,1,2)
    plot(t,p.seismograms_y,'k');
    axis([t(1)-0.1*(t(end)-t(1)) t(end)+0.1*(t(end)-t(1)) -1.1*m 1.1*m]);
    text(t(1),0.9*m, [station ' (\Delta=' num2str(Delta) '°)']);
    text(t(1),0.7*m, ['source: lat=' num2str(p.s_x) '°, lon=' num2str(p.s_y) '°, depth=' num2str(p.s_z/1000) 'km']);
    text(t(1),0.5*m, ['receiver: lat=' num2str(p.r_x) '°, lon=' num2str(p.r_y) '°, depth=' num2str(p.r_z/1000) 'km']);
    xlabel('time [s]');
    ylabel('\phi velocity [m/s]');
    
    subplot(3,1,3)
    plot(t,p.seismograms_z,'k');
    axis([t(1)-0.1*(t(end)-t(1)) t(end)+0.1*(t(end)-t(1)) -1.1*m 1.1*m]);
    text(t(1),0.9*m, [station ' (\Delta=' num2str(Delta) '°)']);
    text(t(1),0.7*m, ['source: lat=' num2str(p.s_x) '°, lon=' num2str(p.s_y) '°, depth=' num2str(p.s_z/1000) 'km']);
    text(t(1),0.5*m, ['receiver: lat=' num2str(p.r_x) '°, lon=' num2str(p.r_y) '°, depth=' num2str(p.r_z/1000) 'km']);
    xlabel('time [s]');
    ylabel('r velocity [m/s]');

end