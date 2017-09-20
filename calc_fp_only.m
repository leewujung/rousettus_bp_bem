% 2015 03 14  3D BEM code directly from example
% 2016 08 03  Parallelized, also modified to work with new calc_bem.m

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    mtx_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling\calc_bem_20160801b_Ra_224_14k';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    mtx_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling\calc_bem_20160801b_Ra_224_14k'];
end

freq_all = 35*1e3;
[~,script_name,~] = fileparts(mtx_path);

[theta,phi] = meshgrid(-90:5:90,0:5:180);  % field point locations
theta = theta(:);
phi = phi(:);
save_path = mtx_path;

time_stamp = datestr(now,'yyyy-mm-dd, HH:MM:SS');  % timestamp

parfor iF=1:length(freq_all)
    % Load param
    param_fname = sprintf('%s_freq%dkHz_param.mat',script_name,freq_all(iF)/1e3);
    [shape,param] = load_shape_param(save_path,param_fname);
    %load(fullfile(save_path,param_fname));
    
    % Define field points
    [x,y,z] = sph2cart(phi/180*pi,theta/180*pi,param.Rfp);
    xyzFP = [x,y,z];
    
    % Calculate field points
    [Afp,Bfp,Cfp]=point(shape.nodesb,shape.topologyb,param.k,xyzFP,param.nsingON); % quadrilateral
    
    save_fp_fname = sprintf('%s_freq%dkHz_fp_mtx.mat',script_name,freq_all(iF)/1e3);
    save_fp(save_path,save_fp_fname,...
        time_stamp,param_fname,theta,phi,x,y,z,xyzFP,...
        Afp,Bfp,Cfp);
    %save(fullfile(save_path,save_fp_fname),'bem_fp','-v7.3');
end

