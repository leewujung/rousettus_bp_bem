% 2015 03 14  3D BEM code directly from example
% 2016 08 03  Modify code for dealing with CT mesh
%             including mesh rotation and parfar save files

% clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

model_shape = 'Ra-colony-0.5mm';
model_shape_msh_file = 'Ra-colony-0.5mm.msh';
date_stamp = datestr(now,'yyyymmdd');  % date
time_stamp = datestr(now,'yyyy-mm-dd, HH:MM:SS');  % timestamp

freq_all = 35*1e3;

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s_%s',script_name,date_stamp,model_shape);

save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

shape.model_shape = model_shape;
shape.msh_file = model_shape_msh_file;
param.code_start_time = time_stamp;


%% Set parameters
R = 1e-3;       % Model shape scale factor
u0 = 1;         % Maximum velocity amplitude
Rfp = 2;        % Radius of the arc of field points (half a circle in the z-y plane)
nsingON = 1;	% Deal with near-singular integrals
% CONDICIONES AMBIENTALES
pa = 101325;	% Presion Atmosferica (Pa)
t = 20;         % Temperatura (ºC)
Hr = 50;        % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000);

param.R = R;
param.u0 = u0;
param.Rfp = Rfp;
param.nsingON = nsingON;
param.pa = pa;
param.t = t;
param.Hr = Hr;
param.rho = rho;
param.c = c;
param.cf = cf;
param.CpCv = CpCv;
param.nu = nu;
param.alfa = alfa;


%% Read nodes and topology
[nodes,elements,elementsQUAD]=readgeomGMSH(model_shape_msh_file); % Load mesh from GMSH file
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

% Rotate around Z-axis
z_rot_rad = -2/180*pi;
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
y_rot_rad = -5/180*pi;
Ry = [cos(y_rot_rad), 0, sin(y_rot_rad);...
      0, 1, 0;...
      -sin(y_rot_rad), 0, cos(y_rot_rad)];
nodes(:,1:3) = (Ry*nodes(:,1:3)')';
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);  % re-center

% Save rotation results and shape
shape.rot.R_prox = R_prox;
shape.rot.R_dist = R_dist;
shape.rot.T = T;
shape.rot.Rx = Rx;
shape.rot.Ry = Ry;
shape.rot.Rz = Rz;
shape.rot.x_rot_rad = x_rot_rad;
shape.rot.y_rot_rad = y_rot_rad;
shape.rot.z_rot_rad = z_rot_rad;
shape.rot.seq = 'Z-X-Y, centering after each';

shape.nodes = nodes;
shape.elements = elements;
shape.elementsQUAD = elementsQUAD;

% Check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements);
close all

M=size(nodesb,1);
N=size(topologyb,1);

shape.nodesb = nodesb;
shape.topologyb = topologyb;
shape.segments = segments;
shape.M = M;
shape.N = N;


%% Calculating BEM coeffs
for iF=1:length(freq_all)
    k=2*pi*freq_all(iF)/c;  % wavenumber
    %param.k = k;  % do this in the save function

    save_param_fname = sprintf('%s_freq%dkHz_param.mat',script_name,freq_all(iF)/1e3);
    save_param(param,shape,k,save_path,save_param_fname);
    %save(fullfile(save_path,save_param_fname),'param','shape');

    % Calculate the BEM matrices and solve the pressures on the surface
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    B=1i*k*rho*c*B;

    save_bem_fname = sprintf('%s_freq%dkHz_bem_mtx.mat',script_name,freq_all(iF)/1e3);
    save_bem(A,B,CConst,save_path,save_bem_fname);
    %save(fullfile(save_path,save_bem_fname),'bem','-v7.3');
end

