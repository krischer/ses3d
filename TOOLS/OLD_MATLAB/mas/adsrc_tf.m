 %- compute TFR via cross-correlation function ---------------------
 
 [tau_cc nu_cc tf_cc]=tfa_cc(t,data,synth_tap_weighted,dt_new,width,'silent');
        
 %- compute TFR of data and synthetics -----------------------------
        
 [tau nu tf_synth_weighted]=tfa(t,synth_tap_weighted,dt_new,width,threshold,'silent');
 [tau nu tf_synth]=tfa(t,synth_tap,dt_new,width,threshold,'silent');
   
 tf_cc_interp=interp2(tau_cc,nu_cc,tf_cc',tau,nu);
 tf_cc=tf_cc_interp';

 %------------------------------------------------------------------
 %- make window function -------------------------------------------
 %------------------------------------------------------------------

 %- noise taper -----------------
 
 m=max(max(abs(tf_cc)))/10;
 W=1-exp(-(abs(tf_cc).^2)/(m^2));
    
 %- high-pass filter ------------
    
 W=W.*(1-exp(-nu'.^2/0.002^2));
    
 %- low-pass filter -------------
    
 W=W.*(exp(-(nu'-0.045).^4/0.045^4).*[nu'>0.045]+[nu'<=0.045]);
 
 %- smoothing -------------------
 
 N=max(size(W));
 for si=1:3
    W(2:N-1,:)=(W(1:N-2,:)+W(3:N,:)+W(2:N-1,:))/3; 
    W(:,2:N-1)=(W(:,1:N-2)+W(:,3:N)+W(:,2:N-1))/3;
 end
 
 %- normalisation ---------------
 
 W=W./max(max(W));
        
 %------------------------------------------------------------------   
 %- compute phase difference ---------------------------------------
 %------------------------------------------------------------------

 DP=imag(log(eps+tf_cc./(eps+abs(tf_cc))));
 
 %------------------------------------------------------------------
 %- detect phase jumps ---------------------------------------------
 %------------------------------------------------------------------
 
 N=max(size(DP));
 
 for si=1:3
    DP(2:N-1,:)=(DP(1:N-2,:)+DP(3:N,:)+DP(2:N-1,:))/3;
    DP(:,2:N-1)=(DP(:,1:N-2)+DP(:,3:N)+DP(:,2:N-1))/3;
 end
 
 test_field=W.*DP./(max(max(abs(W.*DP))));
 criterion_1=max(max(abs(diff(test_field))));
 criterion_2=max(max(abs(diff(test_field'))));
 criterion=max(criterion_1,criterion_2);
 ad_kill=1.0;
 
 %if (criterion>0.7)
 if (criterion>1000.7)
     
     ad_kill=0.0;
     fprintf(1,'possible phase jump\n');
     
 end
 
 %------------------------------------------------------------------
 %- compute phase misfit -------------------------------------------
 %------------------------------------------------------------------

 dnu=nu(2,1)-nu(1,1);
 
 Ep=sqrt(sum(sum(W.*W.*DP.*DP))*dt_new*dnu);
 
 if (isnan(Ep))
     ad_kill=1;
     Ep=0.0;
 end
 
 fprintf(fid_err,'phase misfit: %g %g %g %g\n',Ep, weight_sub, weight_geo, weight_event);
 
 if (Ep>misfit_upper_bound)
     ad_kill=0.0;
     fprintf(1,'misfit exceeds upper bound: set to zero!\n');
 end
 
 fprintf(1,'phase misfit (unweighted): %g\n',Ep);
       
 %------------------------------------------------------------------
 %- make kernel for the inverse tf transform -----------------------
 %------------------------------------------------------------------
 
 IDP=W.*W.*DP.*tf_synth_weighted./(eps+abs(tf_synth).*abs(tf_synth));
 
 %------------------------------------------------------------------
 %- plot results ---------------------------------------------------
 %------------------------------------------------------------------
 
 if (1==0)
        
     %- plot weighted phase misfit ----------------------------------------
     
     load colormap_tf      
     pcolor(tau,nu,W'.*DP');
     colormap(colormap_tf);
     shading flat
     hold on
     m=max(max(abs(W'.*DP')));
     [c h]=contour(tau,nu,W'.*DP',[-1:0.2:1]*m,'LineColor',[0.1 0.1 0.1]);
     colorbar
     caxis([-m m])
     axis([0 max(t) 0 0.12])
     title('weighted phase difference')
     hold off
     
     if (strcmp(mode,'manual'))
        input('press ENTER to continue ');
     else
         pause(1.0)
     end
     
     %- plot window function ----------------------------------------------
     
     load colormap_mono
     pcolor(tau,nu,W');
     colormap(colormap_mono);
     shading flat
     hold on
     [c h]=contour(tau,nu,W',[-1:0.2:1]*max(max(abs(W))),'LineColor',[0.1 0.1 0.1]);
     colorbar
     axis([0 max(t) 0 0.12])
     title('weighting function')
     hold off

     if (strcmp(mode,'manual'))
        input('press ENTER to continue ');
     else
         pause(1.0)
     end
     
      %- plot integral kernel ---------------------------------------------
     
     if (1==0) 
      
     load colormap_tf      
     pcolor(tau,nu,abs(IDP'));
     colormap(colormap_tf);
     shading flat
     hold on
     m=max(max(abs(IDP')));
     colorbar
     caxis([-m m])
     axis([0 max(t) 0 0.12])
     title('integral kernel')
     hold off
     
     if (strcmp(mode,'manual'))
        input('press ENTER to continue ');
     else
         pause(1.5)
     end
     
     end
     
 end
 
 %------------------------------------------------------------------
 %- invert tf transform and make adjoint source --------------------
 %------------------------------------------------------------------
 
 [ad_src,it,I]=itfa(tau,nu,IDP,width,threshold,'silent');

 ad_src=interp1(tau(1,:),imag(ad_src),t,'spline');
       
 ad_src=-ad_src/(Ep+eps);                                     %- divide by misfit
 ad_src=diff(ad_src)/(t(2)-t(1)); ad_src(end+1)=0.0;    %- differentiate
 
 ad_src=ad_kill*ad_src;                                 %- kill if phase jump
 
 if (strcmp(plot_mode,'plot'))
 
     plot(t,ad_src,'k');
     title('adjoint source function')
 
     pause(1.1)
     
 end
 
 iv=length(t):-1:1;                                     %- reverse time
 ad_src=ad_src(iv);
 
 %------------------------------------------------------------------
 %- assign to components -------------------------------------------
 %------------------------------------------------------------------
 
 ad_src_x=0*t;
 ad_src_y=0*t;
 ad_src_z=0*t;
       
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
    
    %- station coordinates in the original coordinate system --------------
    
    colat=90-locations.stla;
    lon=locations.stlo;
    
    %- station coordinates in the rotated coordinate system ---------------
    
    [colat_new lon_new]=rotate_vector(rot_axis,rot_angle,colat,lon,'silent');
    
    %- convert to radians -------------------------------------------------
    
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