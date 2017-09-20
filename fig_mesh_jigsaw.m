% 2016 08 10  Ise JIGSAW to plot bat head mesh

clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\MATLAB\jigsaw');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
else
    save_base_path = ['C:\Users\',usrn,'\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\jigsaw']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\brewermap']);
end

mesh_path = 'F:\Dropbox\0_CODE\OpenBEM_02-2015\3DBEM\input';

% mesh_file = 'bullethead_new_0.1.stl';
% view_angle = [85 45];
% flips = [3,1,2];

% mesh_file = 'Ra224_015k_0.8mm.stl';
% view_angle = [5,25];
% flips = [1,2,3];

% mesh_file = 'Ra-colony-0.5mm.stl';
% view_angle = [-30,25];
% flips = [1,2,3];

% mesh_file = 'Ra-colony-rotear-0.8mm.stl';
% view_angle = [-30,25];
% flips = [1,2,3];

mesh_file = 'Ra-colony-noear-0.8mm.stl';
view_angle = [-30,25];
flips = [1,2,3];

% mesh_file = 'bullethead_sc_colony_0.8mm.stl';
% view_angle = [-30,25];
% flips = [3,2,1];

ec = [0 0 0];
fc = [114 213 223]/255;

ss = strsplit(mesh_file,'.');
mesh_name = strjoin({ss{1:end-1}},'.');
mesh_file_type = ss{end};

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = [script_name,'_',mesh_name];
save_path = fullfile(save_base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Make MSH file for JIGSAW
geom = readstl(fullfile(mesh_path,mesh_file));
makemsh(fullfile(save_path,[mesh_name,'.msh']),geom);

%-- setup files for JIGSAW
opts.geom_file = ...            % geom file
    fullfile(save_path,[mesh_name,'.msh']);
opts.jcfg_file = ...            % config file
    fullfile(save_path,[mesh_name,'.jig']);
opts.mesh_file = ...            % output file
    fullfile(save_path,[mesh_name,'.msh']);

%-- read GEOM file for display
geom = readmsh(fullfile(save_path,[mesh_name,'.msh']));

% %-- meshing options for JIGSAW
% opts.mesh_kern = 'delaunay';
% opts.mesh_dims = 2 ;
% opts.hfun_hmax = 0.03 ;
% mesh_delaunay = jigsaw(opts) ;
% 
% %-- meshing options for JIGSAW
% opts.mesh_kern = 'delfront';
% opts.mesh_dims = 2 ;
% opts.hfun_hmax = 0.03 ;
% mesh_delfront = jigsaw(opts) ;

%-- draw the output
draw_jigsaw_mesh(geom,flips,fc,ec,view_angle);
title(regexprep(['Input geometry: ',mesh_name],'_','\\_'),...
    'fontsize',18);

% draw_jigsaw_mesh(mesh_delaunay,flips,fc,ec,view_angle);
% title('Delaunay mesh');
% 
% draw_jigsaw_mesh(mesh_delfront,flips,fc,ec,view_angle);
% title('Delfront mesh');

saveas(gcf,fullfile(save_path,[mesh_name,'.fig']),'fig');
epswrite(fullfile(save_path,[mesh_name,'.eps']));
saveSameSize(gcf,'file',fullfile(save_path,[mesh_name,'.png']),...
    'format','png','renderer','painters');


