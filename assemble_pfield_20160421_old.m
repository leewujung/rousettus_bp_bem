% 2016 04 21  Assemble pfield from individual src

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

indiv_src_path = 'calc_indiv_src_pfield_20160421_bullethead';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};

% Individual src locations
nn_all = [344, 337, 331, 325, 12];    % midline left
% nn_all = [251, 244, 238, 232, 13];    % midline right
% nn_all = [344, 337, 331, 325, 12,...  % midline left
%           251, 244, 238, 232, 13];    % midline right

freq_all = 35e3;
cvec = 0:-3:-30;

% Tongue clicking location
tongue_loc_all = 0.005;
tongue_loc_all = [-0.02*ones(length(tongue_loc_all),1),...
                  tongue_loc_all',...
                  0*ones(length(tongue_loc_all),1)];
% tongue_loc_all = [0*ones(length(tongue_loc_all),1),...
%                   tongue_loc_all',...
%                   -0.02*ones(length(tongue_loc_all),1)];

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq_all(1)/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results

fig_all = figure;
for iT=1:size(tongue_loc_all,1)
    load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param
        
    tongue_loc = tongue_loc_all(iT,:);

    % Plot tongue and src location on model
    figure(fig_all)
    subplot(length(freq_all),3,1:3:length(freq_all)*3)
    plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
    hold on
    hsrc = plot3(shape.nodesb(nn_all,1),shape.nodesb(nn_all,2),shape.nodesb(nn_all,3),...
                'ko','linewidth',2);
    hton = plot3(tongue_loc(:,1),tongue_loc(:,2),tongue_loc(:,3),...
                'r*','linewidth',2,'markersize',10);
%     ht = text(tongue_loc(:,1),tongue_loc(:,2)+0.002,tongue_loc(:,3),...
%              sprintf('%5.3f,%5.3f,%5.3f',tongue_loc(:,1),tongue_loc(:,2),tongue_loc(:,3)));
% 	set(ht,'backgroundcolor','w','fontsize',12);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    grid on
    title(model_shape)
    tongue_legend = sprintf('tongue loc (%5.3f,%5.3f,%5.3f)',...
                            tongue_loc(:,1),tongue_loc(:,2),tongue_loc(:,3));
    ll = legend([hsrc,hton],{'src loc',tongue_legend},'location','southoutside');
    set(ll,'fontsize',12);
    view([-90 90]);
    
    for iF=1:length(freq_all)
        
        k = 2*pi*freq_all(iF)/param.c;
        
        pfield = nan(size(bem_results.pfield));
        for iN=1:length(nn_all)
            % load individual src pfield
            indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
                indiv_src_path,nn_all(iN),freq_all(iF)/1e3);
            load(fullfile(base_path,indiv_src_path,indiv_src_file));
            
            % phase delay
            gap_loc = shape.nodesb(nn_all(iN),1:3);
            src_gap_dist = sqrt((tongue_loc-gap_loc)*(tongue_loc-gap_loc)');
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
        pp = pp-max(pp(:));
        
%         % Define rotation matrix
%         R_prox = [1,0,0;...
%                   0,1,0;...
%                   0,0,1]';
% %         R_dist = [0,1,0;...
% %                   0,0,1;...
% %                   1,0,0]';
%         R_dist = [0,0,1;...
%                   0,-1,0;...
%                   1,0,0]';
%         T = R_dist'*R_prox;
        
        figure
        axesm eckert4;
        contourfm(bem_results.theta,-bem_results.phi,pp);
        surfm(bem_results.theta,-bem_results.phi,pp);

        % Rotate fieldpoints to new coordinate system
        xyzFP_rot = (T*bem_results.xyzFP')';
        [lon,lat,~] = cart2sph(xyzFP_rot(:,1),xyzFP_rot(:,2),xyzFP_rot(:,3));
        lon = reshape(lon,size(bem_results.theta))/pi*180;
        lat = reshape(lat,size(bem_results.phi))/pi*180;
        
        [lat_new,lon_new] = rotatem(lat,lon,[0,90],'forward','degree');
        axesm eckert4
        [x,y] = mfwdtran(lat,lon);
        [x,y] = mfwdtran(lat_new,lon_new);
        
        figure
        axesm eckert4;
        contourfm(lat_new,lon_new,pp);
        
        lon = reshape(lon,size(bem_results.theta))/pi*180-90;
        lon_correct = lon;
        lon_correct(lon_correct<=-180) = lon_correct(lon_correct<=-180)+360;
        lat = reshape(lat,size(bem_results.phi))/pi*180;
        [x,y] = mfwdtran(defaultm('eckert4'),lat,lon);
        
        figure(fig_all)
        subplot(length(freq_all),3,(iF-1)*3+2:(iF-1)*3+3)
        figure
        axesm eckert4;
%         axesm ortho;
%         gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        surfm(lat,lon,pp);
        surfm(lat,lon_correct,pp);
        contourfm(lat,lon_correct,pp);  % don't plot 0 dB contour
        colormap(parula(length(cvec)-1))
        colorbar
        caxis(cvec([end 1]))
        tightmap
        
        figure(fig_all)
        subplot(length(freq_all),3,(iF-1)*3+2:(iF-1)*3+3)
%         contourf(lon/pi*180-90,lat/pi*180,pp,[43:-3:0]);
        contourf(lon/pi*180-90,lat/pi*180,pp,cvec(2:end));
        grid on
        xlim([-180 180])
        set(gca,'xtick',-180:30:180,'ytick',-90:30:90);
        colormap(parula(length(cvec)-1))
        colorbar
        caxis(cvec([end 1]))
        %         figure(fig2)
        %         polar(bem_fp.theta(19,:)/180*pi,pp(19,:));
    end
end

