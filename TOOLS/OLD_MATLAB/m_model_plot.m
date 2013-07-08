function m_model_plot(model,coordinate,value,quakes)

% coordinate='depth' [km], 'latitude' [deg], 'longitue' [deg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'/home/fichtner/mmap/m_map');
path(path,'/home/fichtner/Matlab/Rotation/');

n_rot=[0 1 0];                      %- coordinate rotation vector
phi_rot=57.5;                       %- rotation angle

phi_rot=0;

figure('Color',[1 1 1]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% depth slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(coordinate,'depth')==1)
    
    %- initialisations ----------------------------------------------------
    
    f_max=-1e99;
    f_min=1e99;
    
    %- make map -----------------------------------------------------------
    
    %m_proj('stereographic','lat',60,'lon',30,'radius',55);
    m_proj('lambert','lon',[21 47],'lat',[32 45]);
    set(gca,'color',[1.0 1.0 1.0]);
    hold on
    
    %- loop over subvolumes -----------------------------------------------
    
    for n=1:model.nsubvol
    
        if ((min(model.m(n).r)<=(6371-value)) && (max(model.m(n).r)>=(6371-value)))
        
            %- coordinate lines -------------------------------------------
        
            d_min=min(abs(6371-value-model.m(n).r));
            idx=min(find((abs(6371-value-model.m(n).r))==d_min));
        
            nx=length(model.m(n).theta);
            ny=length(model.m(n).phi);
        
            [lat lon]=meshgrid(90-model.m(n).theta(1:nx),model.m(n).phi(1:ny));
    
            %- transform coordinate system ------------------------------------
       
            lon_t=zeros(size(lat));
            colat_t=zeros(size(lat));
       
            for ilat=1:nx
                for ilon=1:ny
                    [colat_t(ilon,ilat) lon_t(ilon,ilat)]=rotate_vector(n_rot,-phi_rot,90-lat(ilon,ilat),lon(ilon,ilat),'silent');
                end
            end
       
            lat=90-colat_t;
            lon=lon_t;
    
            %- find extremal values ---------------------------------------
            
            f_max=max(f_max,max(max(abs(model.m(n).v(:,:,idx)))));
            f_min=min(f_min,min(min(abs(model.m(n).v(:,:,idx)))));
            
            %- plot -----------------------------------------------------------
    
            M=zeros(nx,ny);
            M(1:nx-1,1:ny-1)=model.m(n).v(:,:,idx);
            
            m_pcolor(lon,lat,M');
    
            if (strcmp(quakes,'yes'))
                plot_quakes;
            end
    
            shading flat
        
        end
    
    end
    
    m_coast('line','Color','k');
    m_grid('linestyle','none','tickdir','out','linewidth',3);
    %m_grid('linest','-','xticklabels',[],'yticklabels',[])
    title([num2str(value) ' km'],'FontSize',24);
    
    f_mean=(f_max+f_min)/2;
    df=f_max-f_mean;
    
    %f_mean=0;
    caxis([f_mean-df f_mean+df]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% latitude slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
elseif (strcmp(coordinate,'latitude')==1)
    
    min_lon=1000;
    max_lon=-1000;
    min_d=10000;
    max_d=-10000;
    
    %- loop over subvolumes -----------------------------------------------
    
    for n=1:model.nsubvol
    
        if ((max(90-model.m(n).theta)>=value) && (min(90-model.m(n).theta)<=value))

            min_lon=min(min_lon,min(model.m(n).phi));
            max_lon=max(max_lon,max(model.m(n).phi));
            min_d=min(min_d,min(6371-model.m(n).r));
            max_d=max(max_d,max(6371-model.m(n).r));
        
            ny=length(model.m(n).phi);
            nz=length(model.m(n).r);
        
            Mplot=zeros(ny,nz);
    
            idx=min(find((abs(90-model.m(n).theta-value))==min(abs(90-model.m(n).theta-value))));
    
            for i=1:ny-1
                for j=1:nz-1
            
                    Mplot(i,j)=model.m(n).v(idx,i,j);
            
                end
            end
    
            [lon r]=meshgrid(model.m(n).phi,model.m(n).r);
            pcolor(lon,6371-r,Mplot');
            hold on
            axis ij
            %shading flat
            
        end
    
    end
        
    %- arrange figure -----------------------------------------------------
    
    xlabel('longitude [�]');
    ylabel('depth [km]');
        
    axis([min_lon max_lon min_d max_d]);
    
    title(['latitude=' num2str(value) ' �']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% longitude slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
elseif (strcmp(coordinate,'longitude')==1)
    
    min_lat=1000;
    max_lat=-1000;
    min_d=10000;
    max_d=-10000;
    
    %- loop over subvolumes -----------------------------------------------
    
    for n=1:model.nsubvol
    
        if ((max(model.m(n).phi)>=value) && (min(model.m(n).phi)<=value))

            min_lat=min(min_lat,min(90-model.m(n).theta));
            max_lat=max(max_lat,max(90-model.m(n).theta));
            min_d=min(min_d,min(6371-model.m(n).r));
            max_d=max(max_d,max(6371-model.m(n).r));
        
            nx=length(model.m(n).theta);
            nz=length(model.m(n).r);
        
            Mplot=zeros(nx,nz);
    
            idx=min(find((abs(model.m(n).phi-value))==min(abs(model.m(n).phi-value))));
    
            for i=1:nx-1
                for j=1:nz-1
            
                    Mplot(i,j)=model.m(n).v(i,idx,j);
            
                end
            end
    
            [lat r]=meshgrid(90-model.m(n).theta,model.m(n).r);
            pcolor(lat,6371-r,Mplot');
            hold on
            axis ij
            %shading flat
            
        end
    
    end
        
    %- arrange figure -----------------------------------------------------
    
    xlabel('latitude [�]');
    ylabel('depth [km]');
        
    axis([min_lat max_lat min_d max_d]);
    
    title(['latitude=' num2str(value) ' �']);

end
  
load colormap_maps_3
colormap(colormap_maps_3);
colorbar