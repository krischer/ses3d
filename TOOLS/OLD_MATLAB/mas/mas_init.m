
%==========================================================================
% initialise "make adjoint sources (mas)"
%==========================================================================

adloc.x=zeros(1,1000);
adloc.y=zeros(1,1000);
adloc.z=zeros(1,1000);

%==========================================================================
% load traveltime curves
%==========================================================================

cd ../ttcurves

load epi
load S
load SS
load P
load PP
load PcP
load ScS
load ScP

cd ../mas

%==========================================================================
% map setup
%==========================================================================

FS=get(0,'ScreenSize');
S_height=FS(4);
S_width=FS(3);

if (strcmp(plot_mode,'plot'))
                
    h_map=figure('Color',[1 1 1],'Position',[10 40 S_height/2.5 S_height/3]);
    h_seis=figure('Position',[10 S_height-10 0.98*S_width S_height/2.2]);

    figure(h_map)
    m_proj('lambert','lon',[21 47],'lat',[27 50]);
    %m_proj('ortho','lat',52,'lon',5);
    m_coast('patch',[0.9 0.9 0.9],'edgecolor','k','LineWidth',2);
    m_grid('linest','-','xticklabels',[],'yticklabels',[])
    hold on

end

%==================================================================
% open old logfile if auto mode
%==================================================================
        
if (strcmp(mode,'auto'))
            
    fn=[path_logfile num2str(event_list(idx_eq)) '/logfile']

    fid_logold=fopen(fn,'r');
    s=fgets(fid_logold);
            
end
        
%==================================================================
% open ad_srcfile and error file
%==================================================================

fn=[output_path num2str(event_list(idx_eq)) '/ad_srcfile'];
fid_src=fopen(fn,'w');

fn=[output_path num2str(event_list(idx_eq)) '/errorfile'];
fid_err=fopen(fn,'w');

fn=[output_path num2str(event_list(idx_eq)) '/logfile'];
fid_log=fopen(fn,'w');
        
fprintf(fid_log,'station accept num_windows left right weight left right window-weight ... total-weight\n');  

%==================================================================
%- read info file -------------------------------------------------
%==================================================================

fn_temp=[path_data num2str(event_list(idx_eq)) period_tag '/info'];
fid_temp=fopen(fn_temp,'r');

nrec=fscanf(fid_temp,'%d',1);
dummy=fgetl(fid_temp);
dummy=fgetl(fid_temp);
dummy=fgetl(fid_temp);
weight_event=fscanf(fid_temp,'%g',1);
fprintf(1,'event weight = %g\n',weight_event);

fclose(fid_temp);

%==================================================================
%- read station list  ---------------------------------------------
%==================================================================

stations=char(1000,4);
components=char(1000,3);

fn_cont=[path_data num2str(event_list(idx_eq)) period_tag '/content'];
fid_cont=fopen(fn_cont,'r');

for n=1:nrec

    dummy=fscanf(fid_cont,'%c',13);
    stations(n,1:4)=dummy(1:4);
    components(n,1:3)=dummy(6:8);
        
end

fclose(fid_cont);
