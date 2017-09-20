% 2016 03 13  Assemble bullethead msh 0.2 results

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
% data_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160314_bem3D_flathead0.1'];
data_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160314_bem3D_flathead0.1';
data_fpre = 'flathead';

% get shape and param
load(fullfile(data_path,[data_fpre,'_freq35kHz_bem_fp_mtx.mat']),'param','shape','bem_fp');

freq_all = 35e3;

% gap index
% nn_all  = [389, 308, 310, 395, 825];
nn_all = [3138, 3141, 3576, 3224, 3148];

src_all = -[0:0.002:0.016];

fig = figure;
for iS=1:length(src_all)
    
%     src_loc = [0,0.005,-0.02];
    src_loc = [0,src_all(iS),-0.02];
    
    ppFP_all = zeros(length(freq_all),length(nn_all),size(bem_fp.Afp,1));
    for iF=1:length(freq_all)
        k = 2*pi*freq_all(iF)/param.c;
        
        for iN=1:length(nn_all)
            % load bem results
            fname = sprintf('%s_nn%04d_freq%dkHz.mat',data_fpre,nn_all(iN),freq_all(iF)/1e3);
            load(fullfile(data_path,fname));
            
            % phase delay
            gap_loc = shape.nodesb(nn_all(iN),1:3);
            src_gap_dist = sqrt((src_loc-gap_loc)*(src_loc-gap_loc)');
            phase_delay = exp(-1i*k*src_gap_dist);
            
            pfield(:,iN) = bem_results.pfield*phase_delay;
        end
        
        pfield_tot = sum(pfield,2);
        
        %     pfield_tot = reshape(pfield_tot,size(theta));
        pfield_tot_dB = 20*log10(abs(pfield_tot));
        pfield_tot_dB = pfield_tot_dB-(max(max(pfield_tot_dB))-3);
        pfield_tot_dB(pfield_tot_dB<-40) = NaN;
        pfield_tot_dB = pfield_tot_dB+40;
        
        % Plot
        figure(fig)
        polar(bem_fp.theta'/180*pi,pfield_tot_dB);
        title(['steer y = ',num2str(src_loc(2))])
        pause(0.5)
    end
end

