function [t data synth nrec locations name comp]=read_data_and_synthetics(path_data, path_synth, event_id,idx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data synthetics and synthetics given the event id for synthetic
% inversions
%
% input: Event ID, index of recording to be read, component (x,y,z)
%
% output:   time axis and time series for data and synthetics
%           number of recordings
%
% last modified: 13 July, 2009 by Andreas Fichtner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read content files
% get number of recordings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stations=char(400,4);
components=char(400,1);

fn=[path_data event_id '/content'];

fid=fopen(fn,'r');

num=fscanf(fid,'%d',1);
%locations.evla=fscanf(fid,'%g',1);
%locations.evlo=fscanf(fid,'%g',1);

fgetl(fid);

n=0;

while (feof(fid)==0)
    
    dummy=fscanf(fid,'%c',7);
    if (numel(dummy)>1)
        
        stations(n+1,1:4)=dummy(1:4);
        components(n+1)=dummy(6);
        n=n+1;
        
    end
    
end

fclose(fid);

name=stations(idx,:);
comp=components(idx);

fn_synth=[path_synth event_id '/' stations(idx,:)];
fn_data=[path_data event_id '/' stations(idx,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'data: %s\n',fn_data);
fprintf(1,'synthetics: %s\n',fn_synth);

nrec=n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read data and synthetics for the given component and station
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_data=z_read_seismograms_auto([path_data event_id '/'],comp,stations(idx,:),'silent');
data=p_data.seismograms(:);

data=data';

p_synth=z_read_seismograms_auto([path_synth event_id '/'],comp,stations(idx,:),'silent');
synth=p_synth.seismograms(:);

t=0:p_synth.dt:(p_synth.nt-1)*p_synth.dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get epicentre and station locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

locations.stla=p_data.r_x(1);
locations.stlo=p_data.r_y(1);

locations.evla=p_data.s_x(1);
locations.evlo=p_data.s_y(1);