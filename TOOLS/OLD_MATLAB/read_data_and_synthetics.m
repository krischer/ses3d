function [t data synth locations name comp]=read_data_and_synthetics(path_data, path_synth, event_id, station_name, channel)

disp(channel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data and synthetics given the event id
%
% input: Event ID, index of recording to be read, component (x,y,z)
%
% output:   time axis and time series for data and synthetics
%           number of recordings
%
% last modified: 29 March, 2010 by Andreas Fichtner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'/home/fichtner/sacmat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_info=[path_data event_id '/info'];

fn_synth=[path_synth event_id '/' station_name];
fn_data=[path_data event_id '/' station_name];

fprintf(1,'data: %s\n',fn_data);
fprintf(1,'synthetics: %s\n',fn_synth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data and synthetics for the given component and station
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- data -------------------------------------------------------------------

p_data=rsac([fn_data '.' channel '.SAC']);
    
data=p_data(:,2);
t_data=p_data(:,1);

%- synthetics -------------------------------------------------------------

p_synth=read_seismograms([path_synth event_id '/'],station_name,'yes','silent');

if (channel(3)=='Z')
    comp='z';
    synth=p_synth.seismograms_z;
elseif (channel(3)=='E')
    comp='y';
    synth=p_synth.seismograms_y;
elseif (channel(3)=='N')
    comp='x';
    synth=p_synth.seismograms_x;
    data=-data;
end

t_synth=0:p_synth.dt:(p_synth.nt-1)*p_synth.dt;

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

%- shift data relative to origin time -------------------------------------

fid_info=fopen(fn_info,'r');
dummy=fgetl(fid_info);
d=fscanf(fid_info,'%g',1);

if (abs(lh(p_data,'O'))<12340)

    t_data=t_data-lh(p_data,'O')+d;

else
    
    t_data=t_data+d;
    
end

%- correct amplitude of the synthetics ------------------------------------

dummy=fgetl(fid_info);
amp_correction=fscanf(fid_info,'%g',1);
fclose(fid_info);

synth=amp_correction*synth;

%- generate common time axis (time axis of the synthetics) ----------------

t=t_synth;

data=interp1(t_data,data,t,'spline');
data=data.*[t<t_data(end)].*[t>t_data(1)];

%- correct to metres ------------------------------------------------------

data=1e-9*data;

%- assign output (maybe redundant) ----------------------------------------

name=station_name;
