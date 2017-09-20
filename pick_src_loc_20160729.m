% 2016 04 21  Code to facilitate identifying src locations
% 2016 05 07  Plot for NIFTI poster
% 2016 07 26  Use tilted head
% 2016 07 28  Use piston model

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

mtx_path = 'calc_bem_20160729_piston';
param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
load(fullfile(base_path,mtx_path,param_file));
nodes = shape.nodesb(:,1:3);

nn_all = find(nodes(:,1)>3.5e-5);

save(fullfile(save_path,'nn_all.mat'),'nn_all');



