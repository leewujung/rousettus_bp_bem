% Calculates three examples:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Radiation of a a piston on a sphere

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
% save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160313_bem3D_bullethead0.2'];
save_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160313_bem3D_bullethead0.2';
save_fpre = 'bullethead';

if ~exist(save_path,'dir')
    mkdir(save_path);
end

% freq_all = [20:5:30,40:5:55]*1e3;


R=0.01;               % Radius of the sphere
u0=1;              % Maximum velocity amplitude
anglepiston=45;    % Angle of the piston, degrees
Rfp=2;           % Radius of the arc of field points (half a circle in the z-y plane)
nsingON=1;         % Deal with near-singular integrals
% CONDICIONES AMBIENTALES
pa = 101325;         % Presion Atmosferica (Pa)
t = 20;              % Temperatura (ºC)
Hr = 50;              % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000);

% freq = 35e3;
freq_all = (20:5:55)*1e3;


% Read nodes and topology. Sphere of radius 1 m.
%nodes=readnodes('sphere.nod');
%elements=readelements('sphere.ele');
%nodes=readnodes('sphere4tri.nod');
%elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4
% [nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file
[nodes,elements,elementsQUAD]=readgeomGMSH(['bullethead_0.2.msh']); % Load mesh from GMSH file
nodes = nodes(:,[3,2,1]);

% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements);

M=size(nodesb,1); N=size(topologyb,1);

for iF=1:length(freq_all)
    
    k=2*pi*freq_all(iF)/c;  % wavenumber
    
    % Calculate the BEM matrices and solve the pressures on the surface
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    B=i*k*rho*c*B;
    
    save_fname = sprintf('bullethead_freq%dkHz_bem_mtx.mat',freq_all(iF)/1e3);
    save(fullfile(save_path,save_fname));

    % Define field points
    % [theta,phi] = meshgrid(-90:10:90,0:10:180);
    [theta,phi] = meshgrid(-90:5:90,90);
    [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
    xyzFP = [x,y,z]*Rfp;
    
    % Calculate field points
    [Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
    
    save_fname = sprintf('bullethead_freq%dkHz_bem_fp_mtx.mat',freq_all(iF)/1e3);
    save(fullfile(save_path,save_fname));
    
    % Velocity
    nn_all = [37    40    43    46    49];
    
    fig = figure;
    for iN=1:length(nn_all)
        
        nn = nn_all(iN);
        vp = zeros(M,1); vp(nn)=u0;
        
        % Pressure on surface
        pp = A\(-B*vp);
        
        % Pressure at field points
        ppFP(iF,iN,:) = (Afp*pp+j*k*rho*c*Bfp*vp)./Cfp;
        
        ppFP_curr_dB(iF,iN,:) = 20*log10(abs(ppFP(iF,iN,:)));
        ppFP_curr_dB_norm(iF,iN,:) = ppFP_curr_dB(iF,iN,:) - min(ppFP_curr_dB(iF,iN,:));
%         ppFP_curr_dB_norm(iF,iN,:) = ppFP_curr_dB_norm(iF,iN,:) - min(ppFP_curr_dB_norm(iF,iN,:));
        
        figure(fig);
        polar(theta'/180*pi,squeeze(ppFP_curr_dB_norm(iF,iN,:)));
        hold on
    end
    title(['freq = ',num2str(freq_all(iF)/1e3),'kHz']);
end
