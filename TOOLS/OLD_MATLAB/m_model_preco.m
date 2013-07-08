function model_preco=m_model_preco(model)

model_preco=model;

for k=1:model.nsubvol
    
    model_preco.m(k).v=sign(model.m(k).v).*sqrt(abs(model.m(k).v));
    
end