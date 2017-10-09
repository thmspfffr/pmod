%% nc_src_powcorr_degree

% compute degree across frequencies of interest.

clear all

% --------------------------------------------------------
% VERSION 1 - same thresholding for both HC and MS
% --------------------------------------------------------
v_in          = 1; % amplitude corr with log-scaling
v             = 1;
foi_range     = unique(round(2.^[1:.25:7]));
NSUBJ         = 18;
gridsize      = 'coarse';
% --------------------------------------------------------
restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

ft_defaults

indir   = '/home/tpfeffer/pconn/proc/src/';
outdir = '/home/tpfeffer/pconn/proc/conn/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];

%%
% load sa_meg_template;
% grid = sa_meg_template.grid_medium;

for ifoi = 1 : length(foi_range)
  

    allpow = zeros(2113,2113,2,18);
    
    for im = 1 : 3% first compute threshold on icond = 2
      
      disp(sprintf('Computing f%d c%d ...',ifoi,icond))
      
      cnt = 0;
      for isubj = SUBJLIST

        % average the two blocks
        for iblock = 1 : 2
          if exist([indir sprintf('pconn_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v_in)])
            load([indir sprintf('pconn_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v_in)])
            br = 0;
            % correction factor for correlations and fisher z-transform
            res(:,:,iblock) = (resout+resout')./2; clear resout
          else
            br = 1;
            break
          end
        end
        % break & continue with next subject
        if br
          continue
        end
        disp(sprintf('Processig s%df%dm%d',isubj,ifoi,im));
        % count only if subject data exists
        cnt = cnt + 1;
        allpow(:,:,im,cnt)	= nanmean(res,3); clear res
      end
      

  end
end

save([outdir sprintf('nc_powcorr_allpow_thresh%d_v%d.mat',fixed_thresh,v)],'allpow','-v7.3')


%% plot degree (with statistical threshold)
% %
% plt = 1;
% 
% if plt == 1
%   NFREQ = 23;
% end
% 
% if plt
%   if v == 3
%     
%     for ifoi = 1 : NFREQ
%       ifoi
%       for icond = 1 : 2
%         
%         load([outdir sprintf('nc_powcorr_threshold_c%d_f%d_thresh%d_v%d.mat',icond,ifoi,fixed_thresh,v)])
%         c(ifoi,icond) = sum(th(:))/(size(th,1)*size(th,1));
%         
%       end
%     end
%     
%     figure; hold on
%     set(gcf,'color','w')
%     
%     plot(log2(1:NFREQ),c(1:end,1),'color',[1 .55 0.2],'LineWidth',4)
%     plot(log2(1:NFREQ),c(1:end,2),'color',[.11 .56 1],'LineWidth',4)
%     
%     ylabel('Number of connections (in %)');
%     xlabel('Carrier frequency (in Hz)');
%     title('Degree Distribution');
%     
%     set(gca,'TickDir','out','XTick',log2([2,4,8,16,32,64]),'XTickLabel',[2,4,8,16,32,64])
%     
%     b = area(log2(1:NFREQ),c(1:end,2),'FaceColor',[.11 .56 1]);
%     a = area(log2(1:NFREQ),c(1:end,1),'FaceColor',[1 .55 0.2]);
%     
%   elseif v == 4
%     
%     load([outdir sprintf('nc_powcorr_allpow_thresh%d_v%d.mat',fixed_thresh,v)])
%     
%     m = nanmean(allpow,3);
%     s = nanstd(allpow,[],3)./sqrt(size(allpow,3));
%     
%     figure; hold on
%     set(gcf,'color','w')
%     
%     plot(1:NFREQ,m(1,:)','color',[1 .55 0.2],'LineWidth',4)
%     plot(1:NFREQ,m(2,:)','color',[.11 .56 1],'LineWidth',4)
%     
%     ylabel('Correlation');
%     xlabel('Carrier frequency (Hz)');
%     title('Average Correlation');
%     
%     set(gca,'TickDir','out','XTick',[1,3,7,11,15,19,23],'XTickLabel',[2,4,8,16,32,64 128])
%     
%     X = [(1:NFREQ),fliplr(1:NFREQ)];
%     Y = [m(1,:)-s(1,:),fliplr(m(1,:)+s(1,:))];
%     b = fill(X',Y,[1 .55 0.2],'EdgeColor','none'); alpha(b,0.2)
%     
%     X = [(1:NFREQ),fliplr(1:NFREQ)];
%     Y = [m(2,:)-s(2,:),fliplr(m(2,:)+s(2,:))];
%     c = fill(X',Y,[.11 .56 1],'EdgeColor','none'); alpha(c,0.2)
%     
%     set(gca,'xlim',[0 25.101])
%     set(gca,'ylim',[-0.0101 0.0502])
%     
%   end
%   
%   for ifoi = 1 : NFREQ
%     [~, p(ifoi)] = ttest(allpow(1,ifoi,:),allpow(2,ifoi,:),'tail');
%   end
%   
%   
%   
%   
%   %
%   %     sem(1,:) = nanstd(c(5:30,1,:),[],3)./sqrt(isubj);
%   %     sem(2,:) = nanstd(c(5:30,2,:),[],3)./sqrt(isubj);
%   %
%   % X1 = [5:30,fliplr(5:30)];
%   % Y1 = [nanmean(c(5:30,1,:),3)'-sem(1,:),fliplr(nanmean(c(5:30,1,:),3)'+sem(1,:))];
%   %
%   % h = figure; hold on;
%   % set(h,'color','w')
%   %
%   % a = plot(5:30,nanmean(c(5:30,1,:),3),'LineWidth',3); hold on
%   % b = patch(X1,Y1,'b','EdgeColor','none'); alpha(b,0.2);
%   %
%   % X1 = [5:30,fliplr(5:30)];
%   % Y1 = [nanmean(c(5:30,2,:),3)'-sem(2,:),fliplr(nanmean(c(5:30,2,:),3)'+sem(2,:))];
%   %
%   % d = plot(5:30,nanmean(c(5:30,2,:),3),'color','r','LineWidth',3); hold on
%   % e = patch(X1,Y1,'r','EdgeColor','none'); alpha(e,0.2);
%   %
%   % legend([a,d],'MS','HC')
%   %
%   % title('Mean correlation across frequency bands')
%   % xlabel('Frequency (in Hz)')
%   % ylabel('Mean correlation')
%   
%   
%   
%   %   pow = nanmean(allpow(1,:,:),3) - nanmean(allpow(2,:,:),3);
%   
%   
%   %   para = [];
%   %   para.mydotmarkersize=20;
%   %   para.orientation='sagittal';
%   %     para.colorlimits = [-0.075 0.075];
%   
%   %   M = [nanmean(nanmean(allpow(1,:,:),3)) nanmean(nanmean(allpow(2,:,:),3))];
%   %   Mall = nanmean(M);
%   
%   
%   %   for j = 1 : 2
%   %     for is = 1 : cnt
%   %       z(is,:,j) = (allpow(j,:,is) - nanmean(allpow(j,:,is),2))./nanstd(allpow(j,:,is),[],2);
%   %     end
%   %
%   %     t = (nanmean(z(:,:,j))-0)./(nanstd(z(:,:,j))/sqrt(is));
%   %     p = tpdf(t,15);
%   %     th = fdr(p,0.05);
%   %
%   %     allpow_plt(j,:) = nanmean(allpow(j,:,:),3);
%   %     allpow_plt(j,~th | allpow_plt(j,:)<M(j)) = 0;
%   %
%   %   end
%   
%   
%   %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid pow'],grid(1,:));
%   %     saveas(h,sprintf([plotdir 'nc_powcorr_diff_f%d_loc%d_v%d.fig'],ifoi,idx,1));
%   %     close
%   
%   %   lim(1) = min([min(allpow_plt(1,:)) min(allpow_plt(2,:))]);
%   %   lim(2) = max([max(allpow_plt(1,:)) max(allpow_plt(2,:))]);
%   
%   %   para.colorlimits = [lim(1) lim(2)];
%   
%   %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(1,:)'],grid(i,:));
%   %     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],1,ifoi,idx,1));
%   %     close
%   
%   %   h=figure;showmri_transp_v3(sa_meg_template.mri,para,[grid allpow_plt(2,:)'],grid(i,:));
%   %     saveas(h,sprintf([plotdir 'nc_powcorr_c%d_f%d_loc%d_v%d.fig'],2,ifoi,idx,1),'fig');
%   %     close
%   
%   % end
%   
%   
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
