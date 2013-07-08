function s=m_model_make_spline(model,h,depth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s.colat=[];
s.lon=[];
s.h=h;

eta1=[];
eta2=[];
eta3=[];

hw=180*acos(2^(2/3)-(0.5/h)*(2^(2/3)-1)*(1+h^2))/pi;    %- halfwidth in deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read model slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npt=0;

for n=1:model.nsubvol
    
        if ((min(model.m(n).r)<=(6371-depth)) && (max(model.m(n).r)>=(6371-depth)))
            
            %- coordinate lines -------------------------------------------
        
            d_min=min(abs(6371-depth-model.m(n).r));
            idx=min(find((abs(6371-depth-model.m(n).r))==d_min));
            
            nx=length(model.m(n).theta)-1;
            ny=length(model.m(n).phi)-1;

            colat=(model.m(n).theta(1:nx)+model.m(n).theta(2:nx+1))/2;
            lon=(model.m(n).phi(1:ny)+model.m(n).phi(2:ny+1))/2;
    
            %- append colat and lon vectors -------------------------------
        
            for ilat=1:nx
                for ilon=1:ny
                    
                    npt=npt+1;
                    
                    s.colat(npt)=colat(ilat);
                    s.lon(npt)=lon(ilon);
                    y(npt)=model.m(n).v(ilat,ilon,idx);
                       
                    eta1(npt)=cos(s.lon(npt)*pi/180).*sin(s.colat(npt)*pi/180);
                    eta2(npt)=sin(s.lon(npt)*pi/180).*sin(s.colat(npt)*pi/180);
                    eta3(npt)=cos(s.colat(npt)*pi/180);
            
                end
            end
            
        end
    
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(npt);
m=(1/(4*pi))*(1-h^2);

tic

for i=1:npt
    
    for k=1:npt
                
        if ((abs(s.lon(i)-s.lon(k))<5*hw) && (abs(s.colat(i)-s.colat(k))<5*hw))
        
            %- scalar product ---------------------------------------------
        
            p=eta1(i)*eta1(k)+eta2(i)*eta2(k)+eta3(i)*eta3(k);
        
            %- evaluate Abel-Poisson kernel -------------------------------
        
            A(i,k)=m./(1+h^2-2*h*p).^(3/2);
            
        end
        
    end
    
end
   
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
s.v=cgs(A,y',[],10);
toc