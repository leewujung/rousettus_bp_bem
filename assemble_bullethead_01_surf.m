% 2016 03 13  Assemble bullethead msh 0.2 results

% The pressure at an arc of field points is also calculated
clear

% Set save path
usrn = getenv('username');
% data_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling\20160317_bem3D_bullethead0.1_surf'];
data_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling\20160317_bem3D_bullethead0.1_surf';
data_fpre = 'bullethead';

freq_all = (25:10:55)*1e3;
% freq_all = 35*1e3;

% gap index
nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644];  % right side 2 rows
% nn_all  = [1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];  % left side 2 rows
% nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644, ...
%            1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];  % all sources

% src_all = -[0:0.002:0.016];
src_all = 0.005;

fig = figure;
for iF=1:length(freq_all)
    % get shape and param
    load(fullfile(data_path,[data_fpre,'_freq',num2str(freq_all(iF)/1e3),'kHz_param.mat']));
    
    for iS=1:length(src_all)
        %     src_loc = [0,0.005,-0.02];
        src_loc = [0,src_all(iS),-0.02];
        
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
        
        pp = reshape(pfield_tot_dB,size(bem_results.phi));
        
        % Define rotation matrix
        R_prox = [1,0,0;...
            0,1,0;...
            0,0,1]';
        R_dist = [0,1,0;...
            0,0,1;...
            1,0,0]';
        T = R_dist'*R_prox;
        
        % Rotate fieldpoints to new coordinate system
        xyzFP_rot = (T*bem_results.xyzFP')';
        [lon,lat,~] = cart2sph(xyzFP_rot(:,1),xyzFP_rot(:,2),xyzFP_rot(:,3));
        lon = reshape(lon,size(bem_results.theta));
        lat = reshape(lat,size(bem_results.phi));
        
%         figure(fig)
%         subplot(4,1,iF);
% %         axesm eckert4;
%         axesm ortho;
%         gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%         surfm(lat/pi*180,lon/pi*180-90,pp);
%         tightmap
        
        figure(fig)
        subplot(length(freq_all),1,iF);
        contourf(lon/pi*180-90,lat/pi*180,pp,[43:-3:0]);
        if iF==1
            title(['steer x = ',num2str(src_loc(2)),' m']);
        end
        xlim([-90 90])
        
        %         figure(fig2)
        %         polar(bem_fp.theta(19,:)/180*pi,pp(19,:));
    end
end

