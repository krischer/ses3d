function n=m_model_norm(model)

n=0;

for i=1:model.nsubvol
    
    V=(model.m(i).theta(2)-model.m(i).theta(1))*(model.m(i).phi(2)-model.m(i).phi(1))*(model.m(i).r(2)-model.m(i).r(1));
    
    n=n+V*sum(sum(sum((model.m(i).v).^2)));
    
end

n=sqrt(n);