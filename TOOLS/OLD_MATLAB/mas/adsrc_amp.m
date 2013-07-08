%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% misfit for rms amplitude ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DAS IST NUR EIN DUMMY. DIE BERECHNUNG DER ADJOINT SOURCE IST NOCH NICHT
%%% KORREKT

%--------------------------------------------------------------------------
%- compute l2 amplitude misfit
%--------------------------------------------------------------------------

Ep=sum(abs(synth_tap_weighted))/sum(abs(data));

%--------------------------------------------------------------------------
%- plot and write results
%--------------------------------------------------------------------------

fprintf(fid_err,'phase misfit: %g %g %g\n',Ep, weight_sub, weight_geo);
fprintf(1,'phase misfit (unweighted): %g\n',Ep);

%--------------------------------------------------------------------------
%- compute adjoint source time function
%--------------------------------------------------------------------------

ad_src_x=0*t;
ad_src_y=0*t;
ad_src_z=0*t;

ad_src=synth_tap/(dt*sum(synth_tap.*synth_tap));
iv=length(t):-1:1;

ad_src=ad_src(iv);

if (comp=='x')
    ad_src_x=ad_src;
elseif (comp=='y')
    ad_src_y=ad_src;
elseif (comp=='z')
    ad_src_z=ad_src;
end

%--------------------------------------------------------------------------
%- rotate adjoint source time function if necessary
%--------------------------------------------------------------------------

if (strcmp(rot_flag,'yes'))
    
    colat=90-locations.stla;
    lon=locations.stlo;
    
    [colat_new lon_new]=rotate_vector(rot_axis,rot_angle,colat,lon,'silent');
    
    colat=pi*colat/180;
    lon=pi*lon/180;
    
    colat_new=pi*colat_new/180;
    lon_new=pi*lon_new/180;
    
    %- unit vectors in original coordinate system -------------------------
   
    e1=[cos(lon)*cos(colat); sin(lon)*cos(colat); -sin(colat)];     % e_theta
    e2=[-sin(lon); cos(lon); 0];                                    % e_phi
    e3=[cos(lon)*sin(colat); sin(lon)*sin(colat); cos(colat)];      % e_r
    
    %- transform original unit vectors to rotated basis -------------------
    
    R=Rot(rot_axis,pi*rot_angle/180);
    
    e1=R*e1;
    e2=R*e2;
    e3=R*e3;
    
    %- unit vectors in rotated coordinate system --------------------------
    
    e1_new=[cos(lon_new)*cos(colat_new); sin(lon_new)*cos(colat_new); -sin(colat_new)];         % e_theta
    e2_new=[-sin(lon_new); cos(lon_new); 0];                                                    % e_phi
    e3_new=[cos(lon_new)*sin(colat_new); sin(lon_new)*sin(colat_new); cos(colat_new)];          % e_r
    
    %- compute transformation matrix --------------------------------------

    A=zeros(3);

    A(1,1)=e1_new'*e1;
    A(1,2)=e1_new'*e2;
    A(1,3)=e1_new'*e3;

    A(2,1)=e2_new'*e1;
    A(2,2)=e2_new'*e2;
    A(2,3)=e2_new'*e3;

    A(3,1)=e3_new'*e1;
    A(3,2)=e3_new'*e2;
    A(3,3)=e3_new'*e3;
    
    %- transform adjoint source time function -----------------------------
    
    ad_src_x_new=ad_src_x;
    ad_src_y_new=ad_src_y;
    ad_src_z_new=ad_src_z;
    
    ad_src_x_new=A(1,1)*ad_src_x+A(1,2)*ad_src_y+A(1,3)*ad_src_z;
    ad_src_y_new=A(2,1)*ad_src_x+A(2,2)*ad_src_y+A(2,3)*ad_src_z;
    ad_src_z_new=A(3,1)*ad_src_x+A(3,2)*ad_src_y+A(3,3)*ad_src_z;

    ad_src_x=ad_src_x_new;
    ad_src_y=ad_src_y_new;
    ad_src_z=ad_src_z_new;
    
end