% 2016 04 21  Assemble pfield from individual src
% 2016 05 07  Plot for NIFTI poster
% 2016 05 20  Plot for ASA talk
% 2016 07 27  Assemble field for tilted bullethead
% 2016 07 28  Assemble field for piston (on square plate)
% 2016 07 29  Assemble field for piston (circular))

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_plot_opt = 1;
freq = 35*1e3;
cvec = 0:-3:-30;
az_plot_range = 90;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 90;
num_freq = length(freq);
sm_len = 10;

indiv_src_path = 'calc_indiv_src_pfield_20160729_piston';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};

A.base_path = base_path;
A.indiv_src_path = indiv_src_path;
A.model_shape = model_shape;


% Load gap locations
mtx_path = 'calc_bem_20160729_piston';
param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
M = load(fullfile(base_path,mtx_path,param_file));

% Index and combine all files
nn_files = dir(fullfile(base_path,indiv_src_path,'*.mat'));

% Load one sample to get info
load(fullfile(base_path,'pick_src_loc_20160729','nn_all.mat'));

fname = sprintf('%s_nn%04d_freq%02dkHz',indiv_src_path,nn_all(1),freq/1e3);
load(fullfile(base_path,indiv_src_path,fname));

% Assemble all pfield
pfield = nan(size(bem_results.pfield));
for iN=1:length(nn_all)
    % load individual src pfield
    fname = sprintf('%s_nn%04d_freq%02dkHz',indiv_src_path,nn_all(iN),freq/1e3);
    load(fullfile(base_path,indiv_src_path,fname));
    pfield(:,iN) = bem_results.pfield;
end

pfield_tot = sum(pfield,2);
pfield_tot_dB = 20*log10(abs(pfield_tot));
pfield_tot_dB = pfield_tot_dB-max(max(pfield_tot_dB));
pp = reshape(pfield_tot_dB,size(bem_results.phi));
pp_plot = pp;
pp_plot(pp_plot<cvec(end)) = cvec(end);
pp_plot(abs(bem_results.phi)>az_plot_range) = NaN;

% % Get -3dB main contour
% [~,vq_norm,azq,elq] = interp_bp(bem_results.phi/180*pi,bem_results.theta/180*pi,pp,'natural');
% azq = azq/pi*180;
% elq = elq/pi*180;
% idx = abs(azq)>90;
% vq_in = vq_norm;
% vq_in(idx) = NaN;
% [~,c_main_nan] = get_main_contour(vq_in,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
% [c3db_x,c3db_y] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]


M.el = bem_results.theta;
M.az = bem_results.phi;
M.pp_plot = pp_plot;

fig=figure;
plot_model_bp_globe(fig,...  % bp globe
    M,cvec,map_proj,az_plot_limit);



