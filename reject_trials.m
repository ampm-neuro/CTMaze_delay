function rejected_trials = reject_trials(type, lost_trials, varargin)
%rejected_trials('stem', lost_trials, out_stem, light_on_sect, accept_light_on)
%rejected_trials('reward', lost_trials, out_rwd)
%rejected_trials('all', lost_trials, out_stem, light_on_sect, accept_light_on, out_rwd)
%
%outputs which trials to reject/exclude

if strcmp(type, 'stem')
    
    out_stem = varargin{1};
    light_on_sect = varargin{2};
    accept_light_on = varargin{3};
    
    rejected_trials = unique([lost_trials; out_stem'; light_on_sect(~ismember(light_on_sect(:,2), accept_light_on), 1)]);
    
elseif strcmp(type, 'reward')
    
    out_rwd = varargin{1};
    if nargin >3
        error('too many input arguments')
    end
    
    rejected_trials = unique([lost_trials; out_rwd']);
    
    
elseif strcmp(type, 'all')
    
    out_stem = varargin{1};
    light_on_sect = varargin{2};
    accept_light_on = varargin{3};
    out_rwd = varargin{4};
    
    rejected_trials = unique([lost_trials; out_stem'; light_on_sect(~ismember(light_on_sect(:,2), accept_light_on), 1); out_rwd']);
    
end
end