function model=m_model_read(directory,filename)

%- open files -------------------------------------------------------------

disp([directory filename])

fid=fopen([directory filename],'r');

fid_x=fopen([directory 'block_x'],'r');
fid_y=fopen([directory 'block_y'],'r');
fid_z=fopen([directory 'block_z'],'r');

%- number of subvolumes ---------------------------------------------------

model.nsubvol=fscanf(fid_x,'%d',1);
model.nsubvol=fscanf(fid_y,'%d',1);
model.nsubvol=fscanf(fid_z,'%d',1);
model.nsubvol=fscanf(fid,'%d',1);

%- loop over subvolumes ---------------------------------------------------

for n=1:model.nsubvol

    %- read coordinates ---------------------------------------------------
    
    nx=fscanf(fid_x,'%d',1);
    ny=fscanf(fid_y,'%d',1);
    nz=fscanf(fid_z,'%d',1);
    
    model.m(n).theta=fscanf(fid_x,'%g',nx);
    model.m(n).phi=fscanf(fid_y,'%g',ny);
    model.m(n).r=fscanf(fid_z,'%g',nz);
    
    %- read model ---------------------------------------------------------
    
    num=fscanf(fid,'%d',1);
    
    if (num==(nx-1)*(ny-1)*(nz-1))
        
        model.m(n).v=zeros(nx-1,ny-1,nz-1);
        
        for i=1:nx-1
            for j=1:ny-1
                for k=1:nz-1
                    
                    model.m(n).v(i,j,k)=fscanf(fid,'%g',1);
                    
                end
            end
        end
        
    else
        
        fprintf(1,'Error. Dimension mismatch.\n');
        break;
        
    end
    
end

%- clean up ---------------------------------------------------------------

fclose(fid_x);
fclose(fid_y);
fclose(fid_z);

fclose(fid);
