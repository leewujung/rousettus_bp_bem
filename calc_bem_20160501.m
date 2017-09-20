% 2015 03 14  3D BEM code directly from example
clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

model_shape = 'bullethead';
freq_all = [25:5:55,20,60]*1e3;
[theta,phi] = meshgrid(-90:5:90,-180:2:180);  % field point locations

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = [script_name,'_',model_shape];

save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Set parameters
R = 0.01;       % Model shape scale factor
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

% Read nodes and topology
[nodes,elements,elementsQUAD]=readgeomGMSH('bullethead_new_0.1.msh'); % Load mesh from GMSH file
nodes(:,1:3)=nodes(:,1:3)*R;    % change head size

% Rotate model shape to a different coordinate system
R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,0,1;...
          0,1,0;...
          -1,0,0]';
T = R_dist'*R_prox;
nodes(:,1:3) = (T*nodes(:,1:3)')';

% Center the mesh location and tilt rostrum downward toward -z direction
nodes_norm = nodes-repmat(mean(nodes,1),size(nodes,1),1);
theta_y = 10/180*pi;
R_y = [cos(theta_y), 0, sin(theta_y);...
       0, 1, 0;...
       -sin(theta_y), 0, cos(theta_y)];
nodes_rot = (R_y*nodes_norm')';
zshift = -0.005;
nodes_rot_zshift = nodes_rot-...
    [zeros(size(nodes_rot,1),2),...
    zshift*ones(size(nodes_rot,1),1)];
nodes = nodes_rot_zshift;

shape.nodes = nodes;
shape.elements = elements;
shape.elementsQUAD = elementsQUAD;

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements);
M=size(nodesb,1); N=size(topologyb,1);

shape.nodesb = nodesb;
shape.topologyb = topologyb;
shape.segments = segments;
shape.M = M;
shape.N = N;


for iF=1:length(freq_all)
    k=2*pi*freq_all(iF)/c;  % wavenumber
    param.k = k;

    save_param_fname = sprintf('%s_freq%dkHz_param.mat',script_name,freq_all(iF)/1e3);
    save(fullfile(save_path,save_param_fname),'param','shape');

    % Calculate the BEM matrices and solve the pressures on the surface
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    B=1i*k*rho*c*B;
    
    bem.A = A;
    bem.B = B;
    bem.CConst = CConst;

    save_bem_fname = sprintf('%s_freq%dkHz_bem_mtx.mat',script_name,freq_all(iF)/1e3);
    save(fullfile(save_path,save_bem_fname),'bem');
end

% for iF=1:length(freq_all)
%     % Define field points
%     [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
%     xyzFP = [x,y,z]*Rfp;
%     
%     bem_fp.theta = theta;
%     bem_fp.phi = phi;
%     bem_fp.x = x;
%     bem_fp.y = y;
%     bem_fp.z = z;
%     bem_fp.xyzFP = xyzFP;
%     
%     % Calculate field points
%     [Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
%     
%     bem_fp.Afp = Afp;
%     bem_fp.Bfp = Bfp;
%     bem_fp.Cfp = Cfp;
%     
%     save_fp_fname = sprintf('%s_freq%dkHz_fp_mtx.mat',script_name,freq_all(iF)/1e3);
%     save(fullfile(save_path,save_fp_fname),'bem_fp');
% end

