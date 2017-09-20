% 2016 04 21  Code to facilitate identifying src locations

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\jigsaw']);
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    mesh_path = 'F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input';
else
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\jigsaw']);
    base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling'];
    mesh_path = ['C:\Users\',usrn,'\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input'];
end

% mtx_path = 'calc_bem_20160917_Ra-colony-rotear-0.5mm';
% mtx_path = 'calc_bem_20160817_Ra-colony-0.5mm';
mtx_path = 'calc_bem_20160829_Ra-colony-noear-0.5mm';
ss = strsplit(mtx_path,'_');
model_shape = ss{end};
calc_bem_date = ss{3};

param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
load(fullfile(base_path,mtx_path,param_file));
nodes = shape.nodesb(:,1:3);

% Set save_path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,...
    sprintf('%s_%s_%s',script_name,calc_bem_date,model_shape));
if ~exist(save_path,'dir')
    mkdir(save_path);
end

[az,el,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
nn_lhalf = az>0;
nn_rhalf = az<0;


% % Use jigsaw to draw mesh for identifying mouth location
% % Make MSH file for JIGSAW
% geom = readstl(fullfile(mesh_path,[model_shape,'.stl']));
% makemsh(fullfile(save_path,[model_shape,'.msh']),geom);
% %-- read GEOM file for display
% geom = readmsh(fullfile(save_path,[model_shape,'.msh']));
% %-- draw the output
% ec = [0 0 0];
% fc = [114 213 223]/255;
% view_angle = [-30,25];
% draw_jigsaw_mesh(geom,[1,2,3],fc,ec,view_angle);


% Plot head mesh
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.','markersize',5);
plot3(nodes(nn_rhalf,1),nodes(nn_rhalf,2),nodes(nn_rhalf,3),'.','markersize',5);
plot3(nodes(nn_lhalf,1),nodes(nn_lhalf,2),nodes(nn_lhalf,3),'.','markersize',5);
xlabel('X');ylabel('Y');zlabel('Z');
axis equal; grid;

% Find corresponding src loc
pos = reshape([cursor_left_new(:).Position],3,[])';
% pos = pos/1e2;
for iD=1:size(pos,1)
    xyz = nodes(:,1:3)-repmat(pos(iD,:),size(nodes,1),1);
    dist_pos = sqrt(diag(xyz*xyz'));
    [~,idx(iD)] = min(dist_pos);
end

% Plot selected src loc
figure
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.','markersize',5);
hold on
axis equal; grid;
xlabel('X');ylabel('Y');zlabel('Z');
plot3(nodes(idx_right,1),nodes(idx_right,2),nodes(idx_right,3),'o','linewidth',2);
text(nodes(idx_right,1),nodes(idx_right,2),nodes(idx_right,3),num2str([1:6]'));
plot3(nodes(idx_left,1),nodes(idx_left,2),nodes(idx_left,3),'o','linewidth',2);
text(nodes(idx_left,1),nodes(idx_left,2),nodes(idx_left,3),num2str([1:6]'));
plot3(nodes(idx_right2,1),nodes(idx_right2,2),nodes(idx_right2,3),'o','linewidth',2);
text(nodes(idx_right2,1),nodes(idx_right2,2),nodes(idx_right2,3),num2str([1:6]'));
plot3(nodes(idx_left2,1),nodes(idx_left2,2),nodes(idx_left2,3),'o','linewidth',2);
text(nodes(idx_left2,1),nodes(idx_left2,2),nodes(idx_left2,3),num2str([1:6]'));

