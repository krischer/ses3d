%==========================================================================
% initialisation
%==========================================================================

x_min_global=inf;
y_min_global=inf;
z_min_global=inf;

x_max_global=-inf;
y_max_global=-inf;
z_max_global=-inf;

%==========================================================================
% read boxfile
%==========================================================================

h=figure('Color','w');
maxval=0.0;

fpt=fopen([path_boxfile 'boxfile'],'r');

for i=1:14
    fgets(fpt);
end

p.pp=fscanf(fpt,'%i',1);
p.px=fscanf(fpt,'%i',1);
p.py=fscanf(fpt,'%i',1);
p.pz=fscanf(fpt,'%i',1);

fgets(fpt);
fgets(fpt);

p.mi=[1:p.pp];
p.mi_theta=[1:p.pp];
p.mi_phi=[1:p.pp];
p.mi_r=[1:p.pp];

p.itheta_min=[1:p.pp];
p.itheta_max=[1:p.pp];
p.iphi_min=[1:p.pp];
p.iphi_max=[1:p.pp];
p.ir_min=[1:p.pp];
p.ir_max=[1:p.pp];

p.theta_min=[1:p.pp];
p.theta_max=[1:p.pp];
p.phi_min=[1:p.pp];
p.phi_max=[1:p.pp];
p.r_min=[1:p.pp];
p.r_max=[1:p.pp];

for i=1:p.pp
    
    p.mi(i)=fscanf(fpt,'%d',1);
    p.mi_theta(i)=fscanf(fpt,'%d',1);
    p.mi_phi(i)=fscanf(fpt,'%d',1);
    p.mi_z(i)=fscanf(fpt,'%d',1);
    p.itheta_min(i)=fscanf(fpt,'%i',1);
    p.itheta_max(i)=fscanf(fpt,'%i',1);
    p.iphi_min(i)=fscanf(fpt,'%i',1);
    p.iphi_max(i)=fscanf(fpt,'%i',1);
    p.ir_min(i)=fscanf(fpt,'%i',1);
    p.ir_max(i)=fscanf(fpt,'%i',1);
    
    p.theta_min(i)=fscanf(fpt,'%lf',1)*180/pi;
    p.theta_max(i)=fscanf(fpt,'%lf',1)*180/pi;
    p.phi_min(i)=fscanf(fpt,'%lf',1)*180/pi;
    p.phi_max(i)=fscanf(fpt,'%lf',1)*180/pi;
    p.r_min(i)=fscanf(fpt,'%lf',1);
    p.r_max(i)=fscanf(fpt,'%lf',1);
    
    fgets(fpt);
    fgets(fpt);
    
end

fclose(fpt);

%==========================================================================
% determine filename
%==========================================================================

plane=input('plane (x (theta), y (phi), z (r)): ','s');
comp=input('component (x (theta), y (phi), z (r)): ', 's');
it=input('iteration: ');

if (strcmp(plane,'x'))
   theta=input('theta [deg]: ');
   theta=theta*pi/180;
end

%==========================================================================
% read slice log file
%==========================================================================

fid_log=fopen([path_slices_displacement '/slice_log_' plane],'r');

    fgetl(fid_log);
    lat_max=90-fscanf(fid_log,'%g',1);
    lat_min=90-fscanf(fid_log,'%g',1);
    lon_min=fscanf(fid_log,'%g',1);
    lon_max=fscanf(fid_log,'%g',1);
    fgetl(fid_log);
    fgetl(fid_log);
    is_diss=fscanf(fid_log,'%d',1);
    
    fgetl(fid_log);
    fgetl(fid_log);
    fgetl(fid_log);
    
    nproc=fscanf(fid_log,'%d',1);
    fgetl(fid_log);
    fgetl(fid_log);
    
    n=fscanf(fid_log,'%d',nproc);

fclose(fid_log);

%==========================================================================
% successively open files
%==========================================================================

for index=1:length(n)

   filename=[path_slices_displacement 'slice_' plane '_u' comp '_' num2str(n(index)) '_' num2str(it)];
   fpt=fopen(filename,'r');

   fprintf(1,'opening file %s\n',filename);

   %==========================================================================
   % load coordinate lines
   %==========================================================================

   filename_xco=[path_coordinates 'xco_' num2str(n(index))];
   filename_yco=[path_coordinates 'yco_' num2str(n(index))];
   filename_zco=[path_coordinates 'zco_' num2str(n(index))];

   fid_xco=fopen(filename_xco,'r');
   fid_yco=fopen(filename_yco,'r');
   fid_zco=fopen(filename_zco,'r');

   xco=fscanf(fid_xco,'%f',inf);
   yco=fscanf(fid_yco,'%f',inf);
   zco=fscanf(fid_zco,'%f',inf);
   
   fclose(fid_xco);
   fclose(fid_yco);
   fclose(fid_zco);
   
   %==========================================================================
   % read and plot fields
   %==========================================================================

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %- slice in theta=const plane ------------------------------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (strcmp(plane,'x'))
       
       %- read field values -----------------------------------------------
       
       field=zeros(length(yco),length(zco));
    
       for k=1:length(zco)
           for j=1:length(yco)
               field(j,k)=fscanf(fpt,'%f',1);
           end
       end
    
       if (maxval<max(max(abs(field))))
           
           maxval=max(max(abs(field)));
           
       end
       
       [YCO,ZCO]=meshgrid(yco,zco);
				
       y=ZCO.*sin(YCO).*sin(theta);
       z=ZCO.*cos(YCO);
       
       %- determine local and global maxima and minima --------------------
       
       y_min_local=min(min(y));
       y_max_local=max(max(y));
       z_min_local=min(min(z));
       z_max_local=max(max(z));
       
       if (y_min_local<y_min_global) y_min_global=y_min_local; end
       if (y_max_local>y_max_global) y_max_global=y_max_local; end
       if (z_min_local<z_min_global) z_min_global=z_min_local; end
       if (z_max_local>z_max_global) z_max_global=z_max_local; end

       %- plot field ------------------------------------------------------
       
       pcolor(y,z,field');
       hold on
       shading flat
       
       axis([y_min_global y_max_global z_min_global z_max_global]);
       axis off
       pause(0.2)
       
       %- plot coordinate lines -------------------------------------------
       
       if (plot_coordinates==1)
       
           plot(y(1,:),z(1,:),'k');
           plot(y(end,:),z(end,:),'k');
           plot(y(:,1),z(:,1),'k');
           plot(y(:,end),z(:,end),'k');       
           
           if (max(zco)==6371e3)
               
               text(6500e3*sin(min(yco))*sin(theta),6500e3*cos(min(yco)),[num2str(round(min(yco*180/pi))) '°'],'FontSize',20);
               text(6500e3*sin(max(yco))*sin(theta),6500e3*cos(max(yco)),[num2str(round(max(yco*180/pi))) '°'],'FontSize',20);
               
           end
           
       end
       
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %- slice in phi=const plane --------------------------------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (strcmp(plane,'y'))
    
       %- read field values -----------------------------------------------
       
       field=zeros(length(xco),length(zco));
    
       for k=1:length(zco)
           for i=1:length(xco)
               field(i,k)=fscanf(fpt,'%f',1);
           end
       end
    
       
       if (maxval<max(max(abs(field))))
           
           maxval=max(max(abs(field)));
           
       end
       
       [XCO,ZCO]=meshgrid(xco,zco);

       x=ZCO.*cos(XCO);
       z=ZCO.*sin(XCO);
       
       %- determine local and global maxima and minima --------------------
       
       x_min_local=min(min(x));
       x_max_local=max(max(x));
       z_min_local=min(min(z));
       z_max_local=max(max(z));
       
       if (x_min_local<x_min_global) x_min_global=x_min_local; end
       if (x_max_local>x_max_global) x_max_global=x_max_local; end
       if (z_min_local<z_min_global) z_min_global=z_min_local; end
       if (z_max_local>z_max_global) z_max_global=z_max_local; end

       
       %- plot field ------------------------------------------------------
       
       pcolor(x,z,field');
       hold on
       shading flat
       
       axis([x_min_global x_max_global z_min_global z_max_global]);
       axis off
       pause(0.2)
       
       %- plot coordinate lines -------------------------------------------
       
       if (plot_coordinates==1)
       
           plot(x(1,:),z(1,:),'k');
           plot(x(end,:),z(end,:),'k');
           plot(x(:,1),z(:,1),'k');
           plot(x(:,end),z(:,end),'k');       
           
           if (max(zco)==6371e3)
               
               text(6500e3*cos(min(xco)),6500e3*sin(min(xco)),[num2str(round(min(xco*180/pi))) '°'],'FontSize',20);
               text(6500e3*cos(max(xco)),6500e3*sin(max(xco)),[num2str(round(max(xco*180/pi))) '°'],'FontSize',20);
               
           end
           
       end
       
      
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %- slice in z=const plane ---------------------------------------------
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (strcmp(plane,'z'))
    
       field=zeros(length(xco),length(yco));
    
       for j=1:length(yco)
           for i=1:length(xco)
               field(i,j)=fscanf(fpt,'%f',1);
           end
       end
    
       
       if (maxval<max(max(abs(field))))
           
           maxval=max(max(abs(field)));
           
       end
       
       [lat,lon]=meshgrid(90-xco*180/pi,yco*180/pi);
       m_pcolor(lon,lat,field');
       shading flat
       hold on
     
   end
   
   fclose(fpt);

end

%=========================================================================
% adjust colors and style
%=========================================================================

load cm_slice
colormap(cm_slice);

if (strcmp(plane,'z'))

    m_coast('patch',[0.7 0.75 0.7],'edgecolor','k','FaceAlpha',0.4);
    m_grid('box','fancy','linestyle','none','FontSize',20);
    caxis([-maxval maxval]);

else
    
    axis equal;
    caxis([-maxval maxval]);
    
end