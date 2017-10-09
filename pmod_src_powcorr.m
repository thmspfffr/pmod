
v_in = 1;

% --------------------------------------------------------
% VERSION 1 - same thresholding for both HC and MS
% --------------------------------------------------------
v_in          = 1; % amplitude corr with log-scaling
v             = 1;
foi_range     = unique(round(2.^[1:.25:7]));
NSUBJ         = 28;
gridsize      = 'cortex';
% --------------------------------------------------------
restoredefaultpath

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath /home/tpfeffer/Documents/MATLAB/fieldtrip-20130925/

tp_addpaths

indir   = '/home/tpfeffer/pmod/proc/src/';
outdir = '/home/tpfeffer/pconn/proc/conn/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

%%
% 
% 

ord       = pconn_randomization;

for ifoi = 1 : 23;
  
%   if ~exist([outdir sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)])

  if strcmp(gridsize,'coarse')
    allpow = zeros(2113,2113,3,18);
  elseif strcmp(gridsize,'cortex')
    allpow = zeros(3000,3000,3,18);
  end

  
  for isubj = SUBJLIST
    for m = 1 : 3 % first compute threshold on icond = 2
      
      im = find(ord(isubj,:)==m);

      disp(sprintf('Computing f%d c%d ...',ifoi,m))
      
      % average the two blocks
      for iblock = 1 : 2
        load([indir sprintf('pconn_powcorr_s%d_m%d_b%d_f%d_v%d.mat',isubj,im,iblock,ifoi,v_in)])
        % correction factor for correlations and fisher z-transform
        tmp(:,:,iblock) = (resout+resout')./2; clear resout
        
      end
      
      allpow(:,:,m,isubj) = nanmean(tmp,3);
      
    end
  end
  save([outdir sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)],'allpow','-v7.3')
  clear allpow
%   end
  
end


%%

v = 1;

clear pos neg t d

for ifoi = 1 : 23;

  
  ifoi 
  
	load([outdir sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)])

  allpow = reshape(single(allpow(:,:,:,SUBJLIST)),[2113*2113 3 18]);
  
  m(ifoi,1) = nanmean(nanmean(allpow(:,1,:),3),1);
  m(ifoi,2) = nanmean(nanmean(allpow(:,2,:),3),1);
  m(ifoi,3) = nanmean(nanmean(allpow(:,3,:),3),1);
  
  s(ifoi,1) = nanstd(nanmean(allpow(:,1,:),1),[],3)/sqrt(18);
  s(ifoi,2) = nanstd(nanmean(allpow(:,2,:),1),[],3)/sqrt(18);
  s(ifoi,3) = nanstd(nanmean(allpow(:,3,:),1),[],3)/sqrt(18);


end
% 
% figure; set(gcf,'color','white'); hold on
% 
% plot(pos,'linewidth',5); 
% plot(neg,'linewidth',5); 
% 
% ylabel('Number of sign. altered connections'); xlabel('Frequency [Hz]');
% title('ATOMOXETINE VS PLACEBO')
% 
set(gca,'XTick',[3 7 11 15 19 23],'XTickLabel',[4 8 16 32 64 128])
% 
% print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_src_powcorr_total_a-p_v%d.eps',v))


%% 
v = 1;

clear pos neg t d

for ifoi = 1 : 23;

  
  ifoi 
  
	load([outdir sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)])

  allpow = allpow(:,:,:,SUBJLIST);
  
  [t,~,~,s]=ttest(allpow(:,:,2,:),allpow(:,:,1,:),'dim',4);
  
  d = nanmean(allpow(:,:,2,:)-allpow(:,:,1,:),4);
  
  pos(ifoi) = sum(nansum(t.*(d>0)));
 	neg(ifoi) = sum(nansum(t.*(d<0)));

end

figure; set(gcf,'color','white'); hold on

plot(pos,'linewidth',5); 
plot(neg,'linewidth',5); 

ylabel('Number of sign. altered connections'); xlabel('Frequency [Hz]');
title('ATOMOXETINE VS PLACEBO')

set(gca,'XTick',[3 7 11 15 19 23],'XTickLabel',[4 8 16 32 64 128])

print(gcf,'-depsc2',sprintf('~/pconn/proc/plots/pconn_src_powcorr_total_a-p_v%d.eps',v))


%%

for ifoi = 1 : 23
  
  rng('default')
  %   pause(ceil(rand*100));
  
  if ~exist([outdir sprintf('pmod_src_powcorr_allpow_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('pmod_src_powcorr_allpow_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
	load([outdir sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  allpow        = single(allpow(:,:,:,1:NSUBJ));
  
  a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:))));
  b     = a(logical(triu(ones(size(a)),1))); clear a
  pn    = [b>0 b<0];
  p     = 2*normcdf(-abs(b));
  cnt 	= [sum(p<thresh) sum(p<thresh)/length(p) sum(p(pn(:,1))<thresh) sum(p(pn(:,2))<thresh)]; clear a
  
  if perm
    
    cntperm = nan(NPERM,4);
    
    switch test
      
      % if statistical test is unpaired ttest
      case 'ttest'
        allpow	= squeeze(cat(3,allpow(:,1,:),allpow(:,2,:)));
        
        for iperm = 1 : NPERM
          
          clear a_perm
          disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
          
          permidx           = randperm(NSUBJ*2);
          [a_perm,~,~,t]    = ttest2(allpow(:,permidx(1:NSUBJ)),allpow(:,permidx(NSUBJ+1:end)),'dim',2);
          t                 = t.tstat;
          cnt_perm(iperm,:)	= [nansum(a_perm(:)) sum(t(logical(a_perm))>0) sum(t(logical(a_perm))<0)];
          
        end
        
        % if statistical test is ranksum test
      case 'ranksum'
        
        allpow	= cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)));
        
        for iperm = 1 : NPERM
          
          clear a_perm
          disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
          
          permidx           = randperm(NSUBJ*2);
          a_perm            = jh_ranksum(cat(3,squeeze(allpow(:,:,permidx(1:NSUBJ))),squeeze(allpow(:,:,permidx(NSUBJ+1:end)))));
          b_perm            = a_perm(logical(triu(ones(size(a_perm)),1))); clear a
          pn_perm           = [b_perm>0 b_perm<0];
          p_perm            = 2*normcdf(-abs(b_perm));
          cnt_perm(iperm,:)	= [sum(p_perm<thresh) sum(p_perm<thresh)/length(p_perm) sum(p_perm(pn_perm(:,1))<thresh) sum(p_perm(pn_perm(:,2))<thresh)];
          
        end
    end
    
    save([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',ifoi,perm_v)],'cnt','cnt_perm','-v7.3');
    
  end
end

error('!!!')
%% MULTIPLE COMPARISONS CORRECTION AND PLOTTING

perm_v  = 3;
v = perm_v;
clear all_cnt acnt

plt = 1;
NFREQ = 23;

  
  clear corp thresh acnt all_cnt srt idx_D idx_R idx_D_max
  
  it = 0;
  for ifoi = 1 : 2 : NFREQ*2
    
    it = it + 1;
    try
    load([outdir sprintf('nc_src_ana_all2all_stat_f%d_v%d.mat',it,perm_v)]);
    all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
    acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
    catch me
    end
  end
  
%   ranks, where 1 = min, NPERM = max
%   for ifreq = 1 : NFREQ*2
%     [~,~,idx_D(:,ifreq)] = unique(all_cnt(:,ifreq));
%   end
%   
 for ifreq = 1 : NFREQ*2
    [idx_D(:,ifreq)] = floor(tiedrank(all_cnt(:,ifreq)));
 end
%   
  % get maximum rank across frequencies and directions
  idx_R     = max(idx_D,[],2);
  idx_D_max = sort(all_cnt,'ascend');
  idx_D_max = idx_D_max(idx_R,:);
  
  % not sure about this part
  % ----------------------------------
  % thresh	= prctile(idx_R,95);
  % all_srt = sort(all_cnt,'ascend');
  % threshs = all_srt(thresh,:);
  % ----------------------------------
  figure; set(gcf,'color','white'); hold on; box on
  % ----------------------------------
  % PLOT CORRECTED P-VALUES
  % ----------------------------------
  % even freqs: HC > MS
  cnt = 0;
  for ifoi = 2 : 2 : 46
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(foi_range),-log10(corp),'b','LineWidth',3); hold on; 
  
  % odd freqs: MS > HC
  cnt = 0;
  for ifoi = 1 : 2 : 45
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/NPERM;
    
  end
  
  plot(log10(foi_range),-log10(corp),'r','LineWidth',3); hold on
  line(log10([foi_range(1) foi_range(end)]),-log10([.05 .05]),'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[0 -log10(0.5) -log10(0.25) 1 2 3],'YTickLabel',[1 0.5 0.25 0.1 0.01 0.001])
  set(gca,'TickDir','out')

  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')
set(gca,'YLim',[0 2.5]);
saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_corr_v%d.fig'],v),'fig')





