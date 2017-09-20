% 2016 03 09 BEM3D debug

usrn = getenv('username');
bem_path1 = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160308_bem3D_bullethead_0.1'];
bem_path2 = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160308_bem3D_bullethead_0.1_Rfp2'];

bem1 = load(fullfile(bem_path1,'head_src_55kHz_bem_mtx.mat'));
bem2 = load(fullfile(bem_path2,'head_src_55kHz_bem_mtx.mat'));

% set source
vp_idx = 954;
vp = zeros(bem1.shape.M,1);
vp(vp_idx) = bem1.param.u0;
vp_theta = atan2(bem1.shape.nodesb(vp_idx,3),bem1.shape.nodesb(vp_idx,1))/pi*180;
vp_phi = atan2(bem1.shape.nodesb(vp_idx,2),bem1.shape.nodesb(vp_idx,1))/pi*180;

% calculate fieldpoints
bem1.pp = bem1.bem.A\(-bem1.bem.B*vp);
bem2.pp = bem2.bem.A\(-bem2.bem.B*vp);
bem1.ppFP = (bem1.bem_fp.Afp*bem1.pp+1j*bem1.param.k*bem1.param.rho*bem1.param.c*bem1.bem_fp.Bfp*vp)./bem1.bem_fp.Cfp;
bem2.ppFP = (bem2.bem_fp.Afp*bem2.pp+1j*bem2.param.k*bem2.param.rho*bem2.param.c*bem2.bem_fp.Bfp*vp)./bem2.bem_fp.Cfp;

% plot to compare
figure;
axesm eckert4;
surfm(bem1.bem_fp.theta,bem1.bem_fp.phi,20*log10(abs(reshape(bem1.ppFP,size(bem1.bem_fp.theta)))));
hold on
plotm(vp_theta,vp_phi,'r*');
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
ttext = sprintf('Rfp %.2f, vp %d, %.2f deg',bem1.param.Rfp,vp_idx,vp_theta);
title(ttext);

figure;
axesm eckert4;
surfm(bem2.bem_fp.theta,bem2.bem_fp.phi,20*log10(abs(reshape(bem2.ppFP,size(bem1.bem_fp.theta)))));
hold on
plotm(vp_theta,vp_phi,'r*');
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
ttext = sprintf('Rfp %.2f, vp %d, %.2f deg',bem2.param.Rfp,vp_idx,vp_theta);
title(ttext);