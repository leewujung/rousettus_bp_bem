% Calculates three examples:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Radiation of a a piston on a sphere

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160309_bem3D_sphere'];
save_fpre = 'sphere';

if ~exist(save_path,'dir')
    mkdir(save_path);
end

freq_all = 35*1e3;

for iF = 1:length(freq_all)
    % Set parameters
    R = 0.01;   % Shape scaling factor
    Rfp = 3;  % Radius of the arc of field points [m]
    u0 = 1;     % Maximum velocity amplitude
    nsingON=1;  % Deal with near-singular integrals
    
    param.R = R;
    param.Rfp = Rfp;
    param.u0 = u0;
    param.nsingON = nsingON;
    
    % ambient conditions
    pa = 101325;  % Presion Atmosferica (Pa)
    t = 20;       % Temperatura (ºC)
    Hr = 50;      % Humedad relativa (%)
    [rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000);
    
    param.pa = pa;
    param.t = t;
    param.Hr = Hr;
    param.rho = rho;
    param.c = c;
    param.cf = cf;
    param.CpCv = CpCv;
    param.nu = nu;
    param.alfa = alfa;
    
    % More parameters
    freq = freq_all(iF);       % acoustic frequency [Hz]
    k = 2*pi*freq/c;               % Wavenumber, m-1
    
    param.freq = freq;
    param.k = k;
    
    % Read nodes and topology. Sphere of radius 1 m.
    %nodes=readnodes('sphere.nod');
    %elements=readelements('sphere.ele');
    %nodes=readnodes('sphere4tri.nod');
    %elements=readelements('sphere4tri.ele');
    % for more nodes use sphere8 or sphere16
    % for quadratic elements use sphqua4
    % [nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file
    shape_file = 'test1_sphere.msh';
    [nodes,elements,elementsQUAD]=readgeomGMSH(shape_file); % Load mesh from GMSH file
    nodes = nodes(:,[3 2 1]);
    
    % Change radius
    nodes(:,1:3)=nodes(:,1:3)*R;
    
    % check geometry and add body numbers
    %[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
    [nodesb,topologyb,segments]=meshcheck(nodes,elements);  % WJL: nodesb is the actually node points used
    
    M=size(nodesb,1); N=size(topologyb,1);
    
    shape.shape_file = shape_file;
    shape.nodes = nodes;
    shape.elements = elements;
    shape.elementsQUAD = elementsQUAD;
    shape.nodesb = nodesb;
    shape.topologyb = topologyb;
    shape.segments = segments;
    shape.M = M;
    shape.N = N;
    
    % Calculate the BEM matrices
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    
    bem.A = A;
    bem.B = B;
    bem.CConst = CConst;
    
    S.param = param;
    S.shape = shape;
    S.bem = bem;
    
    save_fname = sprintf('%s_%dkHz_bem_mtx.mat',save_fpre,freq/1e3);
    save(fullfile(save_path,save_fname),'-struct','S');
    
end



