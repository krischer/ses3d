%==========================================================================
% input
%==========================================================================

%- set paths --------------------------------------------------------------

path(path,'../');
path(path,'/home/fichtner/Matlab/TFA');
path(path,'/home/fichtner/mmap/m_map');

path_data='/home/fichtner/Data_Europa/';                                      % path of real data

path_synth='/nfs/stig/fichtner/SES3D_3.MER/DATA/OUTPUT/T80.13b/';              % path of new synthetics
output_path='/nfs/stig/fichtner/SES3D_3.MER/ADJOINT/13b.tf/';              % path of new adjoint sources

path_synth_old='/nfs/stig/fichtner/SES3D_3.MER/DATA/OUTPUT/T80.10_new/';           % path of old synthetics
path_logfile='/nfs/stig/fichtner/SES3D_3.MER/ADJOINT/10.tf/';              % path of old logfile

%- parameters for the computation of phase misfits ------------------------

    taper_width=25.0;           % Breite des tapers [s]
    width=45.0;                 % width of the Gaussian for the tf transform [s]
    dt_new=3.0;                 % new time increment for the tf transform [s]
    nu_max=1/(3*dt_new);        % maximum frequency, only for plotting purposes
    threshold=0.0;
    misfit_upper_bound=1.5;    % misfits above this bound are set to zero
    
%- data specifications ----------------------------------------------------
                                          
period_tag='.50s';         % period extension of data directory

event_list=[8 11 24 26 27 33 35 38 39 45 46 47 48];


%- rotation parameters ----------------------------------------------------

rot_flag='yes';
rot_axis=[0 1 0]';
rot_angle=57.5;

%- misfit specification ---------------------------------------------------

misfit='tf';

%- corrections ------------------------------------------------------------

station_correction=1;