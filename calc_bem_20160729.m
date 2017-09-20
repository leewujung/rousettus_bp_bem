% 2015 03 14  3D BEM code directly from example
clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

model_shape = 'piston';
model_shape_msh_file = 'piston_circ_0.05.msh';
model_manipulation = 'centering mesh';
time_stamp = datestr(now,'yyyy-mm-dd, HH:MM:SS');  % timestamp

freq_all = 35*1e3;

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s',script_name,model_shape);

save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

shape.model_shape = model_shape;
shape.msh_file = model_shape_msh_file;
shape.model_manipulation = model_manipulation;
param.code_start_time = time_stamp;


%% Set parameters
R = 0.001;       % Model shape scale factor
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

% Rotate model shape to a different coordinate system
R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,0,1;...
          0,1,0;...
          -1,0,0]';
T = R_dist'*R_prox;
nodes(:,1:3) = (T*nodes(:,1:3)')';

% Center the mesh location
nodes = nodes-repmat(mean(nodes,1),size(nodes,1),1);

shape.nodes = nodes;
shape.elements = elements;
shape.elementsQUAD = elementsQUAD;

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements);
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

