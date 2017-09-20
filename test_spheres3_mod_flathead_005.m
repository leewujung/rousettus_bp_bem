% Calculates three examples:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Radiation of a a piston on a sphere

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
% save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160313_bem3D_bullethead0.2'];
save_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160314_bem3D_flathead0.2';
save_fpre = 'flathead';

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

freq_all = 35e3;
% freq_all = (20:5:55)*1e3;


% Read nodes and topology. Sphere of radius 1 m.
%nodes=readnodes('sphere.nod');
%elements=readelements('sphere.ele');
%nodes=readnodes('sphere4tri.nod');
%elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4
% [nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file
[nodes,elements,elementsQUAD]=readgeomGMSH(['flathead_0.2.msh']); % Load mesh from GMSH file
nodes = nodes(:,[3,1,2]);

% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements);

M=size(nodesb,1); N=size(topologyb,1);

% Select nodes
[az,el,r] = cart2sph(nodesb(:,1),nodesb(:,2),nodesb(:,3));

% % Below are for msh 0.05
nn = find(az>-92/180*pi & az<-88/180*pi);
nn = find(az<92/180*pi & az>88/180*pi);

nn = find(el<65/180*pi & el>64/180*pi);
nn = find(el<50/180*pi & el>49/180*pi);
nn = find(el<30/180*pi & el>29/180*pi);
nn = find(el<5/180*pi & el>3/180*pi);
nn = find(el<-15/180*pi & el>-16/180*pi);

% % Below are for msh 0.1
% nn = find(az>-95/180*pi & az<-85/180*pi);
% nn = find(az<95/180*pi & az>85/180*pi);
% 
% nn = find(el<65/180*pi & el>62/180*pi);
% nn = find(el<50/180*pi & el>47/180*pi);
% nn = find(el<30/180*pi & el>27/180*pi);
% nn = find(el<5/180*pi & el>2/180*pi);
% nn = find(el<-15/180*pi & el>-17/180*pi);


% nn_all  = [158, 341, 162, 210, 203, 875, 1046, 1003, 906, 990];

for iF=1:length(freq_all)
    
    k=2*pi*freq_all(iF)/c;  % wavenumber
    
    % Calculate the BEM matrices and solve the pressures on the surface
    [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
    B=i*k*rho*c*B;
    
    save_fname = sprintf('%s_freq%dkHz_bem_mtx.mat',save_fpre,freq_all(iF)/1e3);
    save(fullfile(save_path,save_fname));

    % Define field points
%     [theta,phi] = meshgrid(-90:5:90,0:10:180);
%     [theta,phi] = meshgrid(-180:5:180,80:10:100);
    [theta,phi] = meshgrid(-180:5:180,0:10:180);
    [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
    xyzFP = [x,y,z]*Rfp;
    
    % Calculate field points
    [Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
    
    save_fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',save_fpre,freq_all(iF)/1e3);
    save(fullfile(save_path,save_fname));
    
    % Velocity
    nn_all = [158, 341, 162, 210, 203];
    
    fig = figure;
    for iN=1:length(nn_all)
        
        nn = nn_all(iN);
        vp = zeros(M,1); vp(nn)=u0;
        
        % Pressure on surface
        pp = A\(-B*vp);
        
        % Pressure at field points
        ppFP = (Afp*pp+j*k*rho*c*Bfp*vp)./Cfp;
        
        ppFP_curr_dB = 20*log10(abs(ppFP));
        ppFP_curr_dB_norm = ppFP_curr_dB - max(ppFP_curr_dB);
        ppFP_curr_dB_norm(ppFP_curr_dB_norm<-30) = NaN;
        ppFP_curr_dB_norm = ppFP_curr_dB_norm+30;
%         ppFP_curr_dB_norm(iF,iN,:) = ppFP_curr_dB_norm(iF,iN,:) - min(ppFP_curr_dB_norm(iF,iN,:));
        dd=reshape(ppFP_curr_dB_norm,size(phi));        

        figure(fig);
        polar((-180:5:180)/180*pi,dd(10,:));
        hold on
        
        figure;
        axesm eckert4;
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        surfm(theta,phi-90,dd);
        caxis([0 30])
    end
    title(['freq = ',num2str(freq_all(iF)/1e3),'kHz']);
end
