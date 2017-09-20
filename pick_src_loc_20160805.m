% 2016 04 21  Code to facilitate identifying src locations

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

mtx_path = 'calc_bem_20160803_Ra224-0.6mm';
ss = strsplit(mtx_path,'_');
model_shape = ss{end};
calc_bem_date = ss{3};

param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
load(fullfile(base_path,mtx_path,param_file));
nodes = shape.nodesb(:,1:3);

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,...
    sprintf('%s_%s_%s',script_name,calc_bem_date,model_shape));

if ~exist(save_path,'dir')
    mkdir(save_path);
end

[az,el,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
nn_lhalf = az>0;
nn_rhalf = az<0;


pos = reshape([cursor_left(:).Position],3,[])';
for iD=1:length(pos)
    xyz = nodes(:,1:3)-repmat(pos(iD,:),size(nodes,1),1);
    dist_pos = sqrt(diag(xyz*xyz'));
    [~,idx(iD)] = min(dist_pos);
end



