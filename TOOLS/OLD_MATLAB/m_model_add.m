function model_c=m_model_add(model_a,model_b)

model_c=model_a;

for k=1:model_a.nsubvol
    
    model_c.m(k).v=model_a.m(k).v+model_b.m(k).v;
    
end