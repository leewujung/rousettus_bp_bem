% 2016 03 23  Save the bem matrices and fp matrices to different files, so
%             tha the file size is more manageable

clear
usrn = getenv('username');
% data_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling\20160319_bem3D_ellipsehead0.1_surf'];
data_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling\20160319_bem3D_ellipsehead0.1_surf';
data_fpre = 'ellipsehead';

freq_all = (20:5:55)*1e3;

for iFR=1:length(freq_all)
    disp(['Separating files at ',num2str(freq_all(iFR)/1e3),'kHz']);
    
    fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',data_fpre,freq_all(iFR)/1e3);
    load(fullfile(data_path,fname));
    
    fname_param = sprintf('%s_freq%dkHz_param.mat',data_fpre,freq_all(iFR)/1e3);
    fname_bem_mtx = sprintf('%s_freq%dkHz_bem_mtx.mat',data_fpre,freq_all(iFR)/1e3);
    fname_fp_mtx = sprintf('%s_freq%dkHz_fp_mtx.mat',data_fpre,freq_all(iFR)/1e3);
    
    save(fullfile(data_path,fname_param),'param','shape');
    save(fullfile(data_path,fname_bem_mtx),'bem');
    save(fullfile(data_path,fname_fp_mtx),'bem_fp');
end
