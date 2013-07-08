function m_model_save(model,filename)

fid=fopen(filename,'w');

fprintf(fid,'%d\n',model.nsubvol);

for n=1:model.nsubvol
    
    fprintf(fid,'%d\n',(length(model.m(n).theta)-1)*(length(model.m(n).phi)-1)*(length(model.m(n).r)-1));
    
    for i=1:length(model.m(n).theta)-1
        for j=1:length(model.m(n).phi)-1
            for k=1:length(model.m(n).r)-1
                
                fprintf(fid,'%g\n',model.m(n).v(i,j,k));
                
            end
        end
    end
end

fclose(fid);