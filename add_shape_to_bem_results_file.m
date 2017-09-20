% 2016 03 16    Add shape info in to the bem_results field for easy
%               assembling of pfield after calculation

% Set save path
usrn = getenv('username');
% data_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160317_bem3D_bullethead0.1_surf'];
data_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160317_bem3D_bullethead0.1_surf';
data_fpre = 'bullethead';

freq_all = (20:5:55)*1e3;

for iFreq=1:length(freq_all)
    
    data_fname = sprintf('%s_nn*_freq%dkHz.mat',data_fpre,freq_all(iFreq)/1e3);
    files = dir(fullfile(data_path,data_fname));
    
    bem_fp_fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',data_fpre,freq_all(iFreq)/1e3);
    bem_fp = load(fullfile(data_path,bem_fp_fname),'bem_fp');
    bem_fp = bem_fp.bem_fp;
    
    for iFile=1:length(files)
        
        A = load(fullfile(data_path,files(iFile).name));
        
        if ~isfield(A.bem_results,'phi')
            disp(['Updating file: ',files(iFile).name]);
            A.bem_results.theta = bem_fp.theta;
            A.bem_results.phi = bem_fp.phi;
            A.bem_results.x = bem_fp.x;
            A.bem_results.y = bem_fp.y;
            A.bem_results.z = bem_fp.z;
            A.bem_results.xyzFP = bem_fp.xyzFP;
        end
        
        save(fullfile(data_path,files(iFile).name),'-struct','A');
        
    end
end