% 2016 08 12  Scale bullethead mesh according to bat head size

clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling'];
end

% [~,script_name,~] = fileparts(mfilename('fullpath'));
% save_path = fullfile(base_path,script_name);
% if ~exist(save_path,'dir')
%     mkdir(save_path);
% end

% CT scan bat Ra224
bat_path = 'calc_bem_20160817_Ra-colony-0.5mm';
bat_file = 'calc_bem_20160817_Ra-colony-0.5mm_freq35kHz_param.mat';
BAT = load(fullfile(base_path,bat_path,bat_file));

% Bullethead_new, original dimension
bul_path = 'calc_bem_20160421_bullethead';
bul_file = 'calc_bem_20160421_bullethead_freq35kHz_param.mat';
BUL = load(fullfile(base_path,bul_path,bul_file));

xscale = range(BAT.shape.nodesb(:,1))/range(BUL.shape.nodesb(:,1));
yscale = range(BAT.shape.nodesb(:,2))/range(BUL.shape.nodesb(:,2));
zscale = range(BAT.shape.nodesb(:,3))/range(BUL.shape.nodesb(:,3));

% Scaled bullethead mesh
bul_sc_path = 'F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input';
bul_sc_file = 'bullethead_sc_colony_0.8mm.msh';

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


figure
plot3(BAT.shape.nodesb(:,1),BAT.shape.nodesb(:,2),BAT.shape.nodesb(:,3),...
    '.','markersize',5);
hold on
plot3(BUL.shape.nodesb(:,1),BUL.shape.nodesb(:,2),BUL.shape.nodesb(:,3),...
    '.','markersize',5);
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.','markersize',5)
legend('CT scan bat','bullethead original','bullethead scaled');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal
grid on

