% 2015 03 14  3D BEM code directly from example
clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    mtx_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling\calc_bem_20160729_piston';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    mtx_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling\calc_bem_20160729_piston'];
end

freq_all = 35*1e3;
[~,script_name,~] = fileparts(mtx_path);

[theta,phi] = meshgrid(-90:2:90,-90:2:90);  % field point locations
save_path = mtx_path;

time_stamp = datestr(now,'yyyy-mm-dd, HH:MM:SS');  % timestamp

for iF=1:length(freq_all)
    % Load param
    param_fname = sprintf('%s_freq%dkHz_param.mat',script_name,freq_all(iF)/1e3);
    load(fullfile(save_path,param_fname));
    
    bem_fp.code_start_time = time_stamp;
    bem_fp.param_fname = param_fname;  % param file used
    
    % Define field points
    [x,y,z] = sph2cart(phi(:)/180*pi,theta(:)/180*pi,1);
    xyzFP = [x,y,z]*param.Rfp;
    
    bem_fp.theta = theta;
    bem_fp.phi = phi;
    bem_fp.x = x;
    bem_fp.y = y;
    bem_fp.z = z;
    bem_fp.xyzFP = xyzFP;
    
    % Calculate field points
    [Afp,Bfp,Cfp]=point(shape.nodesb,shape.topologyb,param.k,xyzFP,param.nsingON); % quadrilateral
    
    bem_fp.Afp = Afp;
    bem_fp.Bfp = Bfp;
    bem_fp.Cfp = Cfp;
    
    save_fp_fname = sprintf('%s_freq%dkHz_fp_mtx.mat',script_name,freq_all(iF)/1e3);
    save(fullfile(save_path,save_fp_fname),'bem_fp');
end

