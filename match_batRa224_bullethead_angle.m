% 2016 08 12  Scale bullethead mesh according to bat head size

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

%% CT scan bat Ra224
bat_path = 'F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input';
bat_file = 'Ra224_038k_0.5mm.msh';

% Center and scale
[nodes,elements,elementsQUAD] = ...   % Load mesh from GMSH file
    readgeomGMSH(fullfile(bat_path,bat_file)); 
R = 1e-3;
nodes(:,1:3)=nodes(:,1:3)*R;    % scale head size
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % center the mesh location

% Rotate model shape to a different coordinate system
% --> result facing +X-axis
R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,-1,0;...
          1,0,0;...
          0,0,1]';
T = R_dist'*R_prox;
nodes(:,1:3) = (T*nodes(:,1:3)')';

% Apply rotation to compensate for scan orientation
% Rotate around Z-axis
z_rot_rad = -7/180*pi;
Rz = [cos(z_rot_rad), -sin(z_rot_rad), 0;...
      sin(z_rot_rad), cos(z_rot_rad), 0;...
      0, 0, 1];
nodes(:,1:3) = (Rz*nodes(:,1:3)')';
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % re-center

% Rotate around X-axis
x_rot_rad = 2/180*pi;
Rx = [1, 0, 0;...
      0, cos(x_rot_rad), -sin(x_rot_rad);...
      0, sin(x_rot_rad), cos(x_rot_rad)];
nodes(:,1:3) = (Rx*nodes(:,1:3)')';
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % re-center

% Rotate around Y-axis
y_rot_rad = -20/180*pi;
Ry = [cos(y_rot_rad), 0, sin(y_rot_rad);...
      0, 1, 0;...
      -sin(y_rot_rad), 0, cos(y_rot_rad)];
nodes(:,1:3) = (Ry*nodes(:,1:3)')';
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % re-center

nodes_bat = nodes;

%% Scaled bullethead mesh
bul_sc_path = 'F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input';
bul_sc_file = 'bullethead_sc_0.6mm.msh';

% Center and scale
[nodes,elements,elementsQUAD] = ...   % Load mesh from GMSH file
    readgeomGMSH(fullfile(bul_sc_path,bul_sc_file)); 
R = 1e-3;
nodes(:,1:3)=nodes(:,1:3)*R;    % scale head size
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % center the mesh location

% Rotate model shape to a different coordinate system
R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,0,1;...
          0,1,0;...
          -1,0,0]';
T = R_dist'*R_prox;
nodes(:,1:3) = (T*nodes(:,1:3)')';

% Rotate around Y-axis
y_rot_rad = 5/180*pi;
Ry = [cos(y_rot_rad), 0, sin(y_rot_rad);...
      0, 1, 0;...
      -sin(y_rot_rad), 0, cos(y_rot_rad)];
nodes(:,1:3) = (Ry*nodes(:,1:3)')';
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % re-center

nodes_bul = nodes;


%% Plot to compare
figure
plot3(nodes_bat(:,1),nodes_bat(:,2),nodes_bat(:,3),...
    '.','markersize',5);
hold on
plot3(nodes_bul(:,1),nodes_bul(:,2),nodes_bul(:,3),...
    '.','markersize',5);
legend('CT scan bat','bullethead scaled');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal
grid on

