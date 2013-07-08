sigma=-4e6;  % step length

%--------------------------------------------------------------------------
%- load descent directions
%--------------------------------------------------------------------------

if (1==1)

    fprintf(1,'load descent directions ...\n');

    grad_csv=m_model_read('GRADIENTS/GRADIENTS.7/grad_csv');
    grad_csh=m_model_read('GRADIENTS/GRADIENTS.7/grad_csh');
    grad_cp=m_model_read('GRADIENTS/GRADIENTS.7/grad_cp');
    grad_rho=m_model_read('GRADIENTS/GRADIENTS.7/grad_rho');
    
end

%--------------------------------------------------------------------------
%- load current model
%--------------------------------------------------------------------------

if (1==1)

    fprintf(1,'load current model ...\n');
    
    d_csv=m_model_read('d_beta_sv_7');
    d_csh=m_model_read('d_beta_sh_7');
    d_cp=m_model_read('d_alpha_7');
    d_rho=m_model_read('d_rho_7');
    
end

%--------------------------------------------------------------------------
%- compute update
%--------------------------------------------------------------------------

d_csv_new=m_model_add(m_model_mult(sigma,grad_csv),d_csv);
d_csh_new=m_model_add(m_model_mult(sigma,grad_csh),d_csh);
d_cp_new=m_model_add(m_model_mult(sigma,grad_cp),d_cp);
d_rho_new=m_model_add(m_model_mult(sigma,grad_rho),d_rho);

%--------------------------------------------------------------------------
%- plot
%--------------------------------------------------------------------------

m_model_plot(d_csv,'depth',100,'no');
m_model_plot(d_csv_new,'depth',100,'no');