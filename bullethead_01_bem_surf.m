% 2015 03 14  3D BEM code directly from example
clear

% Set save path
usrn = getenv('username');
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160317_bem3D_bullethead0.1_surf'];
% save_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160317_bem3D_bullethead0.1_surf';
save_fpre = 'bullethead';

if ~exist(save_path,'dir')
    mkdir(save_path);
end


R = 0.01;       % Radius of the sphere
u0 = 1;         % Maximum velocity amplitude
Rfp = 2;        % Radius of the arc of field points (half a circle in the z-y plane)
nsingON = 1;	% Deal with near-singular integrals
% CONDICIONES AMBIENTALES
pa = 101325;	% Presion Atmosferica (Pa)
t = 20;         % Temperatura (ºC)
Hr = 50;        % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000);

% freq_all = 35e3;
freq_all = (20:5:55)*1e3;

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

% Read nodes and topology. Sphere of radius 1 m.
[nodes,elements,elementsQUAD]=readgeomGMSH(['bullethead_new_0.1.msh']); % Load mesh from GMSH file
nodes(:,1:3)=nodes(:,1:3)*R;    % change head size

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
    
    % Calculate the BEM matrices and solve the pressures on the surface
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    B=1i*k*rho*c*B;
    
    bem.A = A;
    bem.B = B;
    bem.CConst = CConst;
    
%     save_fname = sprintf('%s_freq%dkHz_bem_mtx.mat',save_fpre,freq_all(iF)/1e3);
%     save(fullfile(save_path,save_fname),'param','shape');

    % Define field points
%     [theta,phi] = meshgrid(-90:5:90,0:10:180);
%     [theta,phi] = meshgrid(-180:5:180,80:10:100);
%     [theta,phi] = meshgrid(-180:2:180,0:5:180);
    [theta,phi] = meshgrid(-90:2:90,0:5:90);
    [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
    xyzFP = [x,y,z]*Rfp;
    
    bem_fp.theta = theta;
    bem_fp.phi = phi;
    bem_fp.x = x;
    bem_fp.y = y;
    bem_fp.z = z;
    bem_fp.xyzFP = xyzFP;
    
    % Calculate field points
    [Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
    
    bem_fp.Afp = Afp;
    bem_fp.Bfp = Bfp;
    bem_fp.Cfp = Cfp;
    
    save_fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',save_fpre,freq_all(iF)/1e3);
    save(fullfile(save_path,save_fname),'param','shape','bem','bem_fp');
    
    % Velocity
    nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644,...
               1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];

    for iN=1:length(nn_all)
        
        nn = nn_all(iN);
        vp = zeros(M,1); vp(nn)=u0;
        
        bem_results.src_index = nn;
        bem_results.vp = vp;
        
        % Pressure on surface
        pp = A\(-B*vp);
        
        bem_results.psurface = pp;
        
        % Pressure at field points
        ppFP = (Afp*pp+1j*k*rho*c*Bfp*vp)./Cfp;
        
        bem_results.pfield = ppFP;
        bem_results.pfield_dB = 20*log10(abs(ppFP));
        bem_results.pfield_dB_norm = bem_results.pfield_dB - max(bem_results.pfield_dB);
        
        bem_results.theta = bem_fp.theta;
        bem_results.phi = bem_fp.phi;
        bem_results.x = bem_fp.x;
        bem_results.y = bem_fp.y;
        bem_results.z = bem_fp.z;
        bem_results.xyzFP = bem_fp.xyzFP;

        save_fname = sprintf('%s_nn%04d_freq%dkHz.mat%',save_fpre,nn,freq_all(iF)/1e3);
        save(fullfile(save_path,save_fname),'bem_results');

%         figure(fig);
%         polar((-180:5:180)/180*pi,dd(10,:));
%         hold on
%         
%         figure;
%         axesm eckert4;
%         gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%         surfm(theta,phi-90,dd);
%         caxis([0 30])
    end

end
