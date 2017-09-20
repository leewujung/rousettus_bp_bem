% 2016 04 21  Code to facilitate identifying src locations
% 2016 05 07  Plot for NIFTI poster

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

mtx_path = 'calc_bem_20160421_bullethead';
param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
load(fullfile(base_path,mtx_path,param_file));
nodes = shape.nodesb(:,1:3);
[az,el,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
nn_rgape = find( ((el>7/180*pi & el<12/180*pi) |...
                 (el>-12/180*pi & el<-7/180*pi)) &...
                 az>-110/180*pi & az<0/180*pi);
nn_lgape = find( ((el>7/180*pi & el<12/180*pi) |...
                 (el>-12/180*pi & el<-7/180*pi)) &...
                 az>0/180*pi & az<110/180*pi);

figure;
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.','markersize',5);
hold on
plot3(nodes(nn_rgape,1),nodes(nn_rgape,2),nodes(nn_rgape,3),'k.','markersize',10,'linewidth',2);
plot3(nodes(nn_lgape,1),nodes(nn_lgape,2),nodes(nn_lgape,3),'r.','markersize',10,'linewidth',2);


