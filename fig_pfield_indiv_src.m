% 2016 08 11  Plot model bp on globe
%             Revised from assemble_pfield_steer_h_20160803

% clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\brewermap']);
end


indiv_src_path = 'calc_indiv_src_pfield_20160803_Ra224-0.6mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160501_bullethead';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];

% Set save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s_%s',script_name,bem_calc_date,model_shape);
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Src locations path/file
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
% nn_all = [SRC.idx_left];
% title_text = sprintf('%s %s, left srcs',bem_calc_date,model_shape);
% save_fname = sprintf('%s_left',script_name);
nn_all = [SRC.idx_right];
title_text = sprintf('%s %s, right srcs',bem_calc_date,model_shape);
save_fname = sprintf('%s_right',script_name);

% Set other params
freq = 55*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 180;
num_freq = length(freq);
sm_len = 10;

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


%% Loop through tongue position

k = 2*pi*freq/param.c;
nn_all_loc = shape.nodesb(nn_all,1:3)*1e2;

% Plot pfield from each individual source
fig = figure('position',[10 50 940 940]);
pfield = nan(size(bem_results.pfield));
num_indiv_src = length(nn_all);
for iN=1:num_indiv_src
    % load individual src pfield
    indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
        indiv_src_path,nn_all(iN),freq/1e3);
    load(fullfile(base_path,indiv_src_path,indiv_src_file));

    pp_plot = 20*log10(abs(bem_results.pfield));
    pp_plot = pp_plot-max(max(pp_plot));
    pp_plot = reshape(pp_plot,pfield_size);
    bem_results.phi = reshape(bem_results.phi,pfield_size);
    bem_results.theta = reshape(bem_results.theta,pfield_size);

    M.el = bem_results.theta;
    M.az = bem_results.phi;
    M.pp_plot = pp_plot;  % assembled pfield with values smaller than cvec(end) removed
   
    figure(fig)
    plot_model_shape_src(subplot(num_indiv_src,3,1+(iN-1)*3),...  % shape, horizontal
        h_outline,shape.nodesb(nn_all,1:3),nan(1,3),h_axlim,'h');
    plot3(nn_all_loc(iN,1),nn_all_loc(iN,2),nn_all_loc(iN,3),...
        'ro','markerfacecolor','r','markersize',4);
    xlabel(''); ylabel('')
    plot_model_shape_src(subplot(num_indiv_src,3,2+(iN-1)*3),...  % shape, vertical
        v_outline,shape.nodesb(nn_all,1:3),nan(1,3),v_axlim,'v');
    plot3(nn_all_loc(iN,1),nn_all_loc(iN,2),nn_all_loc(iN,3),...
        'ro','markerfacecolor','r','markersize',4);
    xlabel(''); ylabel('')
    plot_model_bp_globe(subplot(num_indiv_src,3,3+(iN-1)*3),...  % bp globe
        M,cvec,map_proj,az_plot_limit);
end
suptitle(title_text)

saveas(fig,fullfile(save_path,[save_fname,'.fig']),'fig');
saveSameSize(fig,'file',fullfile(save_path,[save_fname,'.png']),...
    'format','png','renderer','painters');

