% Calculates three examples:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Radiation of a a piston on a sphere

% The pressure at an arc of field points is also calculated
clear

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

freq = 35e3;
k=2*pi*freq/c;               % Wavenumber, m-1


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


% Velocity (zeroth-order pulsating sphere)
vz=ones(M,1)*u0;

% Velocity (first-order oscillating sphere)
vf=u0*nodesb(:,3)/R;

% Velocity (piston on sphere)
% nn=find(acos(nodesb(:,3)/R)<=anglepiston*pi/180);

% v0 = [0,0,1];
% vv = nodesb(:,1:3)./repmat(sqrt(diag(nodesb(:,1:3)*nodesb(:,1:3)')),1,3);
% tt = acos(vv*v0');
% nn = find(tt/pi*180<anglepiston);

nn = 42;  % for mesh 0.2
% nn = 181;  % for mesh 0.1

vp=zeros(M,1); vp(nn)=u0;

% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
B=i*k*rho*c*B;

% Define field points
[theta,phi] = meshgrid(-90:10:90,0:10:180);  % for mesh 0.2
% [theta,phi] = meshgrid(-90:10:90,90:10:270);  % for mesh 0.1
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
xyzFP = [x,y,z]*Rfp;

% % Field points
% Mfp=50;  % Number of field points
% theta=linspace(0,2*pi,Mfp)';
% % phi = linspace(0,pi,10);
% phi = pi/6;
% [TT,PP] = meshgrid(theta,phi);
% % xyzFP=Rfp*[zeros(Mfp,1) sin(theta) cos(theta)];
% xyzFP=Rfp*[cos(PP(:)).*sin(TT(:)) sin(PP(:)).*sin(TT(:)) cos(TT(:))];

% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
%[Afp,Bfp,Cfp]=point(nodesb,topologyb,topologyb,ones(N,1),k,xyzFP,nsingON);  % triangular

% Pressure on surface
pz=A\(-B*vz); % the three test cases
pf=A\(-B*vf);
pp=A\(-B*vp); 

% Pressure at field points
pzFP=(Afp*pz+j*k*rho*c*Bfp*vz)./Cfp;
pfFP=(Afp*pf+j*k*rho*c*Bfp*vf)./Cfp;
ppFP=(Afp*pp+j*k*rho*c*Bfp*vp)./Cfp;


figure;
axesm eckert4;
% framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%     contourfm(theta_fp,phi_fp,pfield_tot_dB,0:-3:-30);
% surfm(flipud(TT(:)/pi*180-90),PP(:)/pi*180,20*log10(abs(ppFP)));
surfm(theta,phi-90,reshape(20*log10(abs(ppFP))-min(20*log10(abs(ppFP))),size(phi)));
colorbar



