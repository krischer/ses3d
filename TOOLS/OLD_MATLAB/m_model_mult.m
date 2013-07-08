function model_new=m_model_add(factor,model)

model_new=model;

for k=1:model.nsubvol
    
    model_new.m(k).v=factor*model.m(k).v;
    
end