function [t data synth_rr synth_tt synth_pp synth_pr synth_tr synth_tp locations name comp]=read_data_and_synthetics_momtens(path_data, path_synth, event_id, station_name, channel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data and synthetics given the event id
%
% input: Event ID, index of recording to be read, component (x,y,z)
%
% output:   time axis and time series for data and synthetics
%           number of recordings
%
% last modified: 23 November, 2010 by Andreas Fichtner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'/home/fichtner/sacmat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_info=[path_data event_id '/info'];

fn_data=[path_data event_id '/' station_name];

fprintf(1,'data: %s\n',fn_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data and synthetics for the given component and station
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- data -------------------------------------------------------------------

p_data=rsac([fn_data '.' channel '.SAC']);
    
data=p_data(:,2);
t_data=p_data(:,1);

%- synthetics -------------------------------------------------------------

p_synth_rr=read_seismograms([path_synth '/Mrr/' event_id '/'],station_name,'yes','silent');
p_synth_tt=read_seismograms([path_synth '/Mtt/' event_id '/'],station_name,'yes','silent');
p_synth_pp=read_seismograms([path_synth '/Mpp/' event_id '/'],station_name,'yes','silent');
p_synth_pr=read_seismograms([path_synth '/Mpr/' event_id '/'],station_name,'yes','silent');
p_synth_tr=read_seismograms([path_synth '/Mtr/' event_id '/'],station_name,'yes','silent');
p_synth_tp=read_seismograms([path_synth '/Mtp/' event_id '/'],station_name,'yes','silent');

if (channel(3)=='Z')
    comp='z';
    synth_rr=p_synth_rr.seismograms_z;
    synth_tt=p_synth_tt.seismograms_z;
    synth_pp=p_synth_pp.seismograms_z;
    synth_pr=p_synth_pr.seismograms_z;
    synth_tr=p_synth_tr.seismograms_z;
    synth_tp=p_synth_tp.seismograms_z;
elseif (channel(3)=='E')
    comp='y';
    synth_rr=p_synth_rr.seismograms_y;
    synth_tt=p_synth_tt.seismograms_y;
    synth_pp=p_synth_pp.seismograms_y;
    synth_pr=p_synth_pr.seismograms_y;
    synth_tr=p_synth_tr.seismograms_y;
    synth_tp=p_synth_tp.seismograms_y;
elseif (channel(3)=='N')
    comp='x';
    synth_rr=p_synth_rr.seismograms_x;
    synth_tt=p_synth_tt.seismograms_x;
    synth_pp=p_synth_pp.seismograms_x;
    synth_pr=p_synth_pr.seismograms_x;
    synth_tr=p_synth_tr.seismograms_x;
    synth_tp=p_synth_tp.seismograms_x;
    data=-data;
end

t_synth=0:p_synth_rr.dt:(p_synth_rr.nt-1)*p_synth_rr.dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get epicentre and station locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locations.evla=lh(p_data,'EVLA');
locations.evlo=lh(p_data,'EVLO');
locations.stla=lh(p_data,'STLA');
locations.stlo=lh(p_data,'STLO');

if (locations.evlo<0)
    
    locations.evlo=360+locations.evlo;
    
end

if (locations.stlo<0)
    
    locations.stlo=360+locations.stlo;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% corrections and equalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- apply amplitude correction ---------------------------------------------

if (channel(3)=='Z')
    comp='z';
elseif (channel(3)=='E')
    comp='y';
elseif (channel(3)=='N')
    comp='x';
end

fn_sw=['/home/fichtner/Data_Europa/DATA_Turkey/STATIONFILES/' station_name '.' comp];

if (exist(fn_sw,'file'))
        
    fid_sw=fopen(fn_sw,'r');
    fgetl(fid_sw);    
    amp_correction=fscanf(fid_sw,'%g',1);
    fclose(fid_sw);
       
    data=amp_correction*data;
    
    fprintf('amplitude correction applied\n');
    
end

%- shift data relative to origin time -------------------------------------

fid_info=fopen(fn_info,'r');
d=fgetl(fid_info);
d=fscanf(fid_info,'%g',1);
fclose(fid_info);

if (abs(lh(p_data,'O'))<12340)

    t_data=t_data-p_data(8,3)+d;

else
    
    t_data=t_data+d;
    
end

%- generate common time axis (time axis of the synthetics) ----------------

t=t_synth;

data=interp1(t_data,data,t,'spline');
data=data.*[t<t_data(end)].*[t>t_data(1)];

%- correct to metres ------------------------------------------------------

data=1e-9*data;

%- assign output (maybe redundant) ----------------------------------------

name=station_name;
