% 2016 03 13  Assemble bullethead msh 0.2 results

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
data_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160313_bem3D_bullethead0.2'];
% data_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160313_bem3D_bullethead0.2';
data_fpre = 'bullethead';

% get shape and param
load(fullfile(data_path,[data_fpre,'_freq20kHz_bem_mtx.mat']));

% gap index
nn_all = [37    40    43    46    49];

for iF=1:length(freq_all)
    
    % load bem results
    fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',data_fpre,freq_all(iF)/1e3);
    load(fullfile(data_path,fname));
    
    
    for iN=1:length(nn_all)
        
        nn = nn_all(iN);
        vp = zeros(M,1); vp(nn)=u0;
        
        % Pressure on surface
        pp = A\(-B*vp);
        
        % Pressure at field points
        ppFP(iF,iN,:) = (Afp*pp+j*k*rho*c*Bfp*vp)./Cfp;
        
        ppFP_curr_dB(iF,iN,:) = 20*log10(abs(ppFP(iF,iN,:)));
        ppFP_curr_dB_norm(iF,iN,:) = ppFP_curr_dB(iF,iN,:) - min(ppFP_curr_dB(iF,iN,:));
%         ppFP_curr_dB_norm(iF,iN,:) = ppFP_curr_dB_norm(iF,iN,:) - min(ppFP_curr_dB_norm(iF,iN,:));
        
    end
end    
   

src_loc = [0,0.02,0];

fig = figure;
for iF=1:length(freq_all)
    k = 2*pi*freq_all(iF)/c;
    
    for iN=1:length(nn_all)
       
%         % field points
%         fp_file = sprintf('%s_nn%04d_freq_%dkHz.mat',fname_pre,nn_all(iN),freq_all(iF)/1e3);
%         load(fullfile(fp_path,fp_file));
        
        % phase delay
        gap_loc = nodesb(nn_all(iN),1:3);
        src_gap_dist = sqrt((gap_loc(:,1)-src_loc(1)).^2 + (gap_loc(:,2)-src_loc(2)).^2);
        phase_delay = exp(-1i*k*src_gap_dist);

        % steered pfield
        pfield(:,iN) = squeeze(ppFP(iF,iN,:))*phase_delay;
    end
    
    pfield_tot = sum(pfield,2);
    
%     pfield_tot = reshape(pfield_tot,size(theta));
    pfield_tot_dB = 20*log10(abs(pfield_tot));
    pfield_tot_dB = pfield_tot_dB-max(max(pfield_tot_dB));
    pfield_tot_dB(pfield_tot_dB<-20) = NaN;
    pfield_tot_dB = pfield_tot_dB+20;
    
    % Geoshow
    figure(fig);
    polar(theta'/180*pi,pfield_tot_dB);
    hold on
    
%     axesm eckert4;
%     % framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-120 120]);
%     gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
% %     contourfm(theta_fp,phi_fp,pfield_tot_dB,0:-3:-30);
%     surfm(theta,phi,pfield_tot_dB);
%     colorbar
%     title([num2str(freq_all(iF)/1e3),'kHz']);
end
