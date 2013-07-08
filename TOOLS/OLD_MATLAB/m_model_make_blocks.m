%==========================================================================
% generate and plot model blocks
% write blockfile
%==========================================================================

%==========================================================================
% input
%==========================================================================

model.nsubvol=7;

model.m(1).r=6021:10:6371;
model.m(1).phi=-15:0.5:32;
model.m(1).theta=87:0.5:115;

model.m(2).r=6021:10:6371;
model.m(2).phi=-15:1:32;
model.m(2).theta=115:1:120;

model.m(3).r=6021:10:6371;
model.m(3).phi=-15:1:32;
model.m(3).theta=50:1:87;

model.m(4).r=6021:10:6371;
model.m(4).phi=-30:1:-15;
model.m(4).theta=50:1:120;

model.m(5).r=6021:10:6371;
model.m(5).phi=32:1:60;
model.m(5).theta=50:1:120;

model.m(6).r=5701:20:6021;
model.m(6).phi=-30:1:60;
model.m(6).theta=50:1:120;

model.m(7).r=4951:50:5701;
model.m(7).phi=-30:2:60;
model.m(7).theta=50:2:120;

%==========================================================================
% output
%==========================================================================

n_blocks=0;

fidx=fopen('block_m_x','w');
fidy=fopen('block_m_y','w');
fidz=fopen('block_m_z','w');

%- number of subvolumes

fprintf(fidx,'%d\n',model.nsubvol);
fprintf(fidy,'%d\n',model.nsubvol);
fprintf(fidz,'%d\n',model.nsubvol);

for n=1:model.nsubvol

    fprintf(fidx,'%d\n',length(model.m(n).theta));
    fprintf(fidy,'%d\n',length(model.m(n).phi));
    fprintf(fidz,'%d\n',length(model.m(n).r));

    n_blocks=n_blocks+length(model.m(n).theta)*length(model.m(n).phi)*length(model.m(n).r);
  
    for k=1:length(model.m(n).r)
        fprintf(fidz,'%f\n',model.m(n).r(k));
    end

    for k=1:length(model.m(n).phi)
        fprintf(fidy,'%f\n',model.m(n).phi(k));
    end

    for k=1:length(model.m(n).theta)
        fprintf(fidx,'%f\n',model.m(n).theta(k));
    end
    
end

fclose(fidx);
fclose(fidy);
fclose(fidz);

fprintf(1,'number of model blocks: %d\n',n_blocks);

%==========================================================================
% generate blank model
%==========================================================================

for n=1:model.nsubvol
    
    nx=length(model.m(n).theta)-1;
    ny=length(model.m(n).phi)-1;
    nz=length(model.m(n).r)-1;
    
    model.m(n).v=rand(nx,ny,nz);
    
end