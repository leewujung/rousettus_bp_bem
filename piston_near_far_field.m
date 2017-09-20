% 2016 03 09 BEM3D debug

usrn = getenv('username');
% bem_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\\20160309_bem3D_sphere'];
bem_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160312_bem3D_piston';
% save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\\20160309_bem3D_sphere_2src'];
save_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160312_bem3D_piston';
save_fpre = 'piston';

if ~exist(save_path,'dir')
    mkdir(save_path);
end

bem = load(fullfile(bem_path,'piston_35kHz_bem_mtx.mat'));

Rfp_all = [0.05,0.1,0.2,0.5,1,2,5];

for iR = 1:length(Rfp_all)
    Rfp = Rfp_all(iR);
    
    % Define field points
    [theta,phi] = meshgrid(0,-180:10:180);
    [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,Rfp);
    xyzFP = [x,y,z];
    
    % Calculate matrices for field points
    [Afp,Bfp,Cfp]=point(bem.shape.nodesb,bem.shape.topologyb,bem.param.k,xyzFP,bem.param.nsingON); % quadrilateral
    
    bem_fp.theta = theta;
    bem_fp.phi = phi;
    bem_fp.xyzFP = xyzFP;
    bem_fp.Afp = Afp;
    bem_fp.Bfp = Bfp;
    bem_fp.Cfp = Cfp;
    
    S.bem_fp = bem_fp;  % save BEM results again
    
    save_fname = sprintf('Rfp_%.2f.mat',Rfp);
    save(fullfile(save_path,save_fname),'-struct','S');
    
    % find source location with a radius
    v0 = [1,0,0];
    vv = bem.shape.nodesb(:,1:3)./repmat(sqrt(diag(bem.shape.nodesb(:,1:3)*bem.shape.nodesb(:,1:3)')),1,3);
    tt = acos(vv*v0');
    nn = find(tt/pi*180<5);
    
    % set source
%     vp_idx = [12,334];
    vp = zeros(bem.shape.M,1);
    vp(nn) = bem.param.u0;
%     vp_theta = atan2(bem.shape.nodesb(vp_idx,3),bem.shape.nodesb(vp_idx,2))/pi*180;
%     vp_phi = atan2(bem.shape.nodesb(vp_idx,2),bem.shape.nodesb(vp_idx,1))/pi*180;
    
    % calculate fieldpoints
    S.vp = vp;
%     S.vp_theta = theta;
%     S.vp_phi = phi;
    S.pp = bem.bem.A\(-bem.bem.B*vp);
    S.ppFP = (bem_fp.Afp*S.pp+1j*bem.param.k*bem.param.rho*bem.param.c*bem_fp.Bfp*vp)./bem_fp.Cfp;

    save_fname = sprintf('Rfp_%.2f_fp.mat',Rfp);
    save(fullfile(save_path,save_fname),'-struct','S');

    ppFP_all(:,iR) = S.ppFP;
end