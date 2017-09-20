% 2016 08 11  Plot model bp on globe
%             Revised from assemble_pfield_steer_h_20160803

clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\brewermap']);
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Individual pfield path
indiv_src_path = 'calc_indiv_src_pfield_20160811_Ra224-0.6mm';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];

% Src locations path/file
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
nn_all = [SRC.idx_left];

% Set other params
freq = 35*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 180;
num_freq = length(freq);
sm_len = 10;

% Tongue clicking location
tongue_loc = [0.022,0.000,-0.004];

save_fname = sprintf('%s_%s_%s_%02dkHz_x%05.1f_y%05.1f_z%05.1f',...
    script_name,bem_calc_date,model_shape,freq/1e3,...
    tongue_loc(1,1)*1e3,tongue_loc(1,2)*1e3,tongue_loc(1,3)*1e3);
title_text = sprintf('%s %s, %dkHz, tongue (%5.1f,%5.1f,%5.1f)',...
    bem_calc_date,model_shape,freq/1e3,...
    tongue_loc(1,1)*1e3,tongue_loc(1,2)*1e3,tongue_loc(1,3)*1e3);

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results
load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param

% Get shape outline
pad = 0.002;
v_outline = get_head_shape_vertical(shape.nodesb(:,1:3));
v_axlim = [min(v_outline(:,1))-pad,...
    max(v_outline(:,1))+pad,...
    -0.1,0.1,...
    min(v_outline(:,3))-pad,...
    max(v_outline(:,3))+pad];

h_outline = get_head_shape_horizontal(shape.nodesb(:,1:3));
h_axlim = [min(h_outline(:,1))-pad,...
    max(h_outline(:,1))+pad,...
    min(h_outline(:,2))-pad,...
    max(h_outline(:,2))+pad,...
    -0.1,0.1];

%% Plot top and side view of head and tongue locations
shape_title_text = sprintf('Tongue (%5.1f,%5.1f,%5.1f)',...
    tongue_loc(1,1)*1e3,...
    tongue_loc(1,2)*1e3,...
    tongue_loc(1,3)*1e3);

fig_all = figure('position',[200,100,1400,420]);
plot_model_shape_src(subplot(141),...  % shape, horizontal
    h_outline,shape.nodesb(nn_all,1:3),nan(1,3),h_axlim,'h');
plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
    'x','linewidth',2,'markersize',8,'color','r');
plot_model_shape_src(subplot(142),...  % shape, vertical
    v_outline,shape.nodesb(nn_all,1:3),nan(1,3),v_axlim,'v');
plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
    'x','linewidth',2,'markersize',8,'color','r');

%% 
k = 2*pi*freq/param.c;

% Assemble total field
pfield = nan(size(bem_results.pfield));
for iN=1:length(nn_all)
    % load individual src pfield
    indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
        indiv_src_path,nn_all(iN),freq/1e3);
    load(fullfile(base_path,indiv_src_path,indiv_src_file));
    
    % phase delay
    gap_loc = shape.nodesb(nn_all(iN),1:3);
    src_gap_dist = sqrt((tongue_loc-gap_loc)*(tongue_loc-gap_loc)');
    phase_delay = exp(1i*k*src_gap_dist);
    
    pfield(:,iN) = bem_results.pfield*phase_delay;
end

pfield_tot = sum(pfield,2);
pfield_tot_dB = 20*log10(abs(pfield_tot));
pfield_tot_dB = pfield_tot_dB-max(max(pfield_tot_dB));
pp = reshape(pfield_tot_dB,size(bem_results.phi));
pp_plot = pp;
pp_plot(pp_plot<cvec(end)) = cvec(end);

pp = reshape(pp,pfield_size);
pp_plot = reshape(pp_plot,pfield_size);
bem_results.phi = reshape(bem_results.phi,pfield_size);
bem_results.theta = reshape(bem_results.theta,pfield_size);

M.el = bem_results.theta;
M.az = bem_results.phi;
M.pp = pp;  % original assembled pfield
M.pp_plot = pp_plot;  % assembled pfield with values smaller than cvec(end) removed

% summary figure
figure(fig_all)
plot_model_bp_globe_surf(subplot(1,4,3:4),...  % bp globe
    M,cvec,map_proj,az_plot_limit);
suptitle(title_text);
% saveas(fig_all,fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig_all,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');
% epswrite(fullfile(save_path,[save_fname,'.eps']));


