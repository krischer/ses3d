function model_smooth=m_model_smooth(model,iterations)

model_smooth=model;

%--------------------------------------------------------------------------
% determine lateral extent of the model
%--------------------------------------------------------------------------

theta_min=1e6;
theta_max=-1e6;

phi_min=1e6;
phi_max=-1e6;

d_theta=zeros(1,model.nsubvol);
d_phi=zeros(1,model.nsubvol);

for i=1:model.nsubvol
    
    d_theta(i)=model.m(i).theta(2)-model.m(i).theta(1);
    d_phi(i)=model.m(i).phi(2)-model.m(i).phi(1);
    
    if (model.m(i).theta(1)<theta_min)
        theta_min=model.m(i).theta(1);
    end
    
    if (model.m(i).theta(end)>theta_max)
        theta_max=model.m(i).theta(end);
    end
    
    if (model.m(i).phi(1)<phi_min)
        phi_min=model.m(i).phi(1);
    end
    
    if (model.m(i).phi(end)>phi_max)
        phi_max=model.m(i).phi(end);
    end

end

%--------------------------------------------------------------------------
% determine finest depth discretisation 
%--------------------------------------------------------------------------

r=model.m(1).r(1);

for i=1:model.nsubvol
    r=[r model.m(i).r'];
end

r=sort(r);

k=1;
for i=2:length(r)
    if (r(i)>r(i-1))
        r_new(k)=r(i);
        k=k+1;
    end
end

r=r_new;

%--------------------------------------------------------------------------
% loop through depth levels
%--------------------------------------------------------------------------

for n=1:length(r)-1
    
    disp(r(n))
    
    %----------------------------------------------------------------------
    % generate the finest lateral grid
    %----------------------------------------------------------------------

    d_theta_min=1e6;
    d_phi_min=1e6;
    
    for k=1:model.nsubvol
        
        if ((model.m(k).r(1)<=r(n)) && (model.m(k).r(end)>r(n)))
             
            if (d_theta(k)<d_theta_min)
                d_theta_min=d_theta(k);
            end
            
            if (d_phi(k)<d_phi_min)
                d_phi_min=d_phi(k);
            end
           
        end
        
    end
    
    [theta,phi]=meshgrid(theta_min:d_theta_min:theta_max,phi_min:d_phi_min:phi_max);
    
    F=0*theta;
    
    %----------------------------------------------------------------------
    % transcribe to the finest grid
    %----------------------------------------------------------------------
    
    for k=1:model.nsubvol
        
        %------------------------------------------------------------------
        % find depth index and transcribe
        %------------------------------------------------------------------
        
        if ((model.m(k).r(1)<=r(n)) && (model.m(k).r(end)>r(n)))
            
            idx=find([abs(model.m(k).r(1:end-1)-r(n))]==[min(abs(model.m(k).r(1:end-1)-r(n)))]);
            
            [theta_old,phi_old]=meshgrid(model.m(k).theta(1:end-1),model.m(k).phi(1:end-1));
            
            F=F+interp2(theta_old,phi_old,model.m(k).v(:,:,idx)',theta,phi,'nearest',0);
           
        end
    end
    
    %----------------------------------------------------------------------
    % smooth
    %----------------------------------------------------------------------
        
    F_old=F;  % redundant, nur zum plotten
    
    for it=1:iterations
        
        F_new=F;
        
        ni=length(F(:,1));
        nj=length(F(1,:));
        
        for i=2:ni-1
            for j=2:nj-1
                
                F_new(i,j)=(2*F(i,j)+F(i+1,j)+F(i-1,j)+F(i,j+1)+F(i,j-1))/6;
                
            end
        end
        
        F(2:ni-1,2:nj-1)=F_new(2:ni-1,2:nj-1);
        
    end
    
    %subplot(1,2,1)
    %pcolor(theta,phi,F_old);
    
    %subplot(1,2,2)
    %pcolor(theta,phi,F);
    
    %----------------------------------------------------------------------
    % transcribe back to the original grid
    %----------------------------------------------------------------------
    
    %figure
    
    theta=theta_min:d_theta_min:theta_max;
    phi=phi_min:d_phi_min:phi_max;
    
    for k=1:model.nsubvol
        
        if ((model.m(k).r(1)<=r(n)) && (model.m(k).r(end)>r(n)))
        
            idx=find([abs(model.m(k).r(1:end-1)-r(n))]==[min(abs(model.m(k).r(1:end-1)-r(n)))]);
            
            [theta_old,phi_old]=meshgrid(model.m(k).theta(1:end-1),model.m(k).phi(1:end-1));
        
            %--------------------------------------------------------------
            % loop over blocks
            %--------------------------------------------------------------
        
            for i=1:length(model.m(k).theta)-1
           
                idx_i=find([theta<model.m(k).theta(i+1)].*[theta>=model.m(k).theta(i)]);
            
                for j=1:length(model.m(k).phi)-1
                
                    idx_j=find([phi<model.m(k).phi(j+1)].*[phi>=model.m(k).phi(j)]);
                
                    model_smooth.m(k).v(i,j,idx)=mean(mean(F(idx_j,idx_i)));
                
                end
            
            end
        
            %m_model_plot(model_smooth,'depth',60,'no');
            
            %disp(k)
            %pcolor(theta_old,phi_old,model.m(k).v(:,:,n)');
            %axis([40 130 -40 70])
            %caxis([-2e-18 2e-18])
            %hold on
            %pause(1.1)
            
        end
        
    end
    
end



