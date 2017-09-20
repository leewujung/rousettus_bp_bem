% 2016 04 21  Code to facilitate identifying src locations

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

% [~,script_name,~] = fileparts(mfilename('fullpath'));
% save_path = fullfile(base_path,script_name);
% if ~exist(save_path,'dir')
%     mkdir(save_path);
% end

mtx_path = '20160317_bem3D_bullethead0.1_surf';
model_shape = 'bullethead';

% Load pre-calculated bem stuff
param_file = sprintf('%s_freq%02dkHz_param.mat',model_shape,35);
load(fullfile(base_path,mtx_path,param_file));

[az,el,r] = cart2sph(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3));
nn_r90 = find(az>-92/180*pi & az<-88/180*pi);  % for bullethead
nn_l90 = find(az>88/180*pi & az<92/180*pi);  % for bullethead

% Plot the model shape
figure;
plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
hold on
plot3(shape.nodesb(nn_r90,1),shape.nodesb(nn_r90,2),shape.nodesb(nn_r90,3),'ro','markersize',5);
plot3(shape.nodesb(nn_l90,1),shape.nodesb(nn_l90,2),shape.nodesb(nn_l90,3),'mo','markersize',5);
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title(model_shape)

R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,0,1;...
          0,1,0;...
          -1,0,0]';
T = R_dist'*R_prox;
shape.nodes_rot = (T*shape.nodes(:,1:3)')';

figure;
plot3(shape.nodesb_rot(:,1),shape.nodesb_rot(:,2),shape.nodesb_rot(:,3),'.','markersize',5);
hold on
plot3(shape.nodesb(nn_r90,1),shape.nodesb(nn_r90,2),shape.nodesb(nn_r90,3),'ro','markersize',5);
plot3(shape.nodesb(nn_l90,1),shape.nodesb(nn_l90,2),shape.nodesb(nn_l90,3),'mo','markersize',5);
plot3(shape.nodesb(idx,1),shape.nodesb(idx,2),shape.nodesb(idx,3),'k*','markersize',10,'linewidth',2);
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title([model_shape,' rotated'])


% Return index at data cursor locations
C = reshape([cursor_info(:).Position],3,length(cursor_info))';
idx = nan(length(C),1);
for iC=1:length(C)
    D = shape.nodesb(:,1:3)-repmat(C(iC,:),length(shape.nodesb),1);
    v = diag(D*D');
    [~,idx(iC)] = min(v);
end

% Plot to check src loc
figure;
plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
hold on
plot3(shape.nodesb(nn_r90,1),shape.nodesb(nn_r90,2),shape.nodesb(nn_r90,3),'ro','markersize',5);
plot3(shape.nodesb(nn_l90,1),shape.nodesb(nn_l90,2),shape.nodesb(nn_l90,3),'mo','markersize',5);
plot3(shape.nodesb(idx,1),shape.nodesb(idx,2),shape.nodesb(idx,3),'k*','markersize',10,'linewidth',2);
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
grid on
title(model_shape)

