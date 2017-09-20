% 2016 03 08  Steer beam in 3DBEM results

clear

usrn = getenv('username');
bem_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160307_bem3D_5src'];
fp_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160307_bem3D_5src'];
fname_pre = 'head_src';

freq_all = 35*1e3;
% freq_all = 35*1e3;
nn_all = [1206, 1223, 1079, 1260, 992];

load(fullfile(bem_path,'head_src_bem_output2.mat'));

src_loc = [0.01,0,0];

for iF=1:length(freq_all)
    for iN=1:length(nn_all)
       
        % field points
        fp_file = sprintf('%s_nn%04d_freq_%dkHz.mat',fname_pre,nn_all(iN),freq_all(iF)/1e3);
        load(fullfile(fp_path,fp_file));
        
        % phase delay
        gap_loc = shape.nodesb(nn_all(iN),1:3);
        src_gap_dist = sqrt((gap_loc(:,1)-src_loc(1)).^2 + (gap_loc(:,2)-src_loc(2)).^2);
        phase_delay = exp(-1i*param.k*src_gap_dist);

        % steered pfield
        pfield(:,iN) = ppFP*phase_delay;
    end
    
    pfield_tot = sum(pfield,2);
    
    pfield_tot = reshape(pfield_tot,size(theta));
    pfield_tot_dB = 20*log10(abs(pfield_tot));
    pfield_tot_dB = pfield_tot_dB-max(max(pfield_tot_dB));
    pfield_tot_dB(pfield_tot_dB<-20) = NaN;
    
    % Geoshow
    figure;
    axesm eckert4;
    % framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%     contourfm(theta_fp,phi_fp,pfield_tot_dB,0:-3:-30);
    surfm(theta,phi,pfield_tot_dB);
    colorbar
    title([num2str(freq_all(iF)/1e3),'kHz']);
end
