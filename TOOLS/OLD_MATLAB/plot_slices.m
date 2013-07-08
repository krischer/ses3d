%==========================================================================
% set paths and small parameters
%==========================================================================

plot_coordinates=1;

path_boxfile='../MODELS/MODELS/';
path_slices_displacement='./';
path_slices_gradient='./';
path_slices_model='../MODELS/MODELS/';
path_coordinates='../DATA/COORDINATES/';

path(path,'/home/fichtner/mmap/m_map');

m_proj('lambert','lon',[21 47],'lat',[27 50]);

%==========================================================================
% call subroutines
%==========================================================================

field_type=input('fields (1=displacement, 2=model, 3=gradients): ');

if (field_type==1)
    
    plot_slices_displacement;
    
elseif (field_type==2)
    
    plot_slices_model;
    
elseif (field_type==3)
    
    plot_slices_gradient
    
end