function m_model_2vtk_point(model,filename)

%--------------------------------------------------------------------------
% initialisation
%--------------------------------------------------------------------------

path(path,'/home/fichtner/Matlab/Rotation/');

n_rot=[0 1 0];                      %- coordinate rotation vector
phi_rot=57.5;                       %- rotation angle

%--------------------------------------------------------------------------
%- number of grid points --------------------------------------------------
%--------------------------------------------------------------------------

nx=zeros(1,model.nsubvol);
ny=zeros(1,model.nsubvol);
nz=zeros(1,model.nsubvol);

for i=1:model.nsubvol

    nx(i)=length(model.m(i).theta)-1;
    ny(i)=length(model.m(i).phi)-1;
    nz(i)=length(model.m(i).r)-1;
    
end

N=0;

for i=1:model.nsubvol
    N=N+nx(i)*ny(i)*nz(i);
end

%--------------------------------------------------------------------------
% open file and write header
%--------------------------------------------------------------------------

fid=fopen(filename,'w');

fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'vtk output\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');

%--------------------------------------------------------------------------
% write grid points
%--------------------------------------------------------------------------

fprintf(fid,'POINTS %d float\n',N);

for n=1:model.nsubvol

    disp(n)
    
    for i=1:nx(n)
        for j=1:ny(n)
            for k=1:nz(n)
            
                %- rotate coordinate system
                
                theta=(model.m(n).theta(i)+model.m(n).theta(i+1))/2;
                phi=(model.m(n).phi(j)+model.m(n).phi(j+1))/2;
                r=(model.m(n).r(k)+model.m(n).r(k+1))/2;
                
                [theta_t phi_t]=rotate_vector(n_rot,-phi_rot,model.m(n).theta(i),model.m(n).phi(j),'silent');
                theta=theta_t*pi/180;
                phi=phi_t*pi/180;
                
                %- transform to cartesian coordinates
                
                x=r*sin(theta)*cos(phi);
                y=r*sin(theta)*sin(phi);
                z=r*cos(theta);
            
                fprintf(fid,'%f %f %f\n',x,y,z);
            
            end
        end
    end
end
  

%--------------------------------------------------------------------------
% write connectivity
%--------------------------------------------------------------------------

%- total number of cells --------------------------------------------------

n_cells=0;

for n=1:model.nsubvol

    n_cells=n_cells+(nx(n)-1)*(ny(n)-1)*(nz(n)-1);
    
end

%- write cell connectivity ------------------------------------------------

fprintf(fid,'\n');
fprintf(fid,'CELLS %d %d\n',n_cells,9*n_cells);

count=0;

for n=1:model.nsubvol

    for i=1:nx(n)-1
        for j=1:ny(n)-1
            for k=1:nz(n)-1
                                                               % i j k
                a=count+k+(j-1)*nz(n)+(i-1)*ny(n)*nz(n)-1;     % 0 0 0
                b=count+k+(j-1)*nz(n)+(i-1)*ny(n)*nz(n);       % 0 0 1
                c=count+k+(j)*nz(n)+(i-1)*ny(n)*nz(n)-1;       % 0 1 0
                d=count+k+(j)*nz(n)+(i-1)*ny(n)*nz(n);         % 0 1 1
                e=count+k+(j-1)*nz(n)+(i)*ny(n)*nz(n)-1;       % 1 0 0
                f=count+k+(j-1)*nz(n)+(i)*ny(n)*nz(n);         % 1 0 1
                g=count+k+(j)*nz(n)+(i)*ny(n)*nz(n)-1;         % 1 1 0
                h=count+k+(j)*nz(n)+(i)*ny(n)*nz(n);           % 1 1 1
            
                fprintf(fid,'8 %d %d %d %d %d %d %d %d\n',a,b,c,d,e,f,g,h);
            
            end
        end
    end
    
    count=count+nx(n)*ny(n)*nz(n);
    
end
 
%--------------------------------------------------------------------------
% write cell types
%--------------------------------------------------------------------------

fprintf(fid,'\n');
fprintf(fid,'CELL_TYPES %d\n',n_cells);

for n=1:model.nsubvol

    for i=1:nx(n)-1
        for j=1:ny(n)-1
            for k=1:nz(n)-1
            
                fprintf(fid,'11\n');
            
            end
        end
    end
end

%--------------------------------------------------------------------------
% write data
%--------------------------------------------------------------------------

fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %d\n',N);
fprintf(fid,'SCALARS scalars float\n');
fprintf(fid,'LOOKUP_TABLE mytable\n');

for n=1:model.nsubvol
    
    for i=1:nx(n)
        for j=1:ny(n)
            for k=1:nz(n)
            
                fprintf(fid,'%f\n',model.m(n).v(i,j,k));
            
            end
        end
    end
end

fclose(fid);















