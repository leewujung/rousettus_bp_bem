

Rfp = 2;
[theta,phi] = meshgrid(-180:20:180,0:20:180);
[x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
xyzFP = [x,y,z]*Rfp;

[tt,pp,rr] = cart2sph(x,y,z);


R_prox = [1,0,0;...
          0,1,0;...
          0,0,1]';
R_dist = [0,1,0;...
          0,0,1;...
          1,0,0]';
T = R_dist'*R_prox;


nodesb = shape.nodesb(:,1:3);
nodesb_rot = (T*nodesb')';


