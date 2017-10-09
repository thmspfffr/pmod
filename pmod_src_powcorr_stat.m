%% COMPUTES STATISTICS IN THE STYLE OF HAWELLEK ET AL. (2013)

% Counts number of altered connections (both directons) and compares
% those against a permutation distribution. Subsequently, the result
% is corrected for multiple comparisons using a variant of the omnibus
% test.

% This version: 24/03/2016

% pmod_src_powcorr_stat

v_in = 1;

% --------------------------------------------------------
% VERSION 1 - same thresholding for both HC and MS
% --------------------------------------------------------
v_in          = 1; % amplitude corr with log-scaling
v             = 1;
foi_range     = unique(round(2.^[1:.25:7]));
NSUBJ         = 28;
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
outdir = '/home/tpfeffer/pmod/proc/';
plotdir = '/home/tpfeffer/pconn/proc/plots/';

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24];



%%
v = 1;
perm_v = 1;
para.perm = 1;
para.nperm = 500;
para.nsubj = 28;
para.thresh = 0.05;

for ifoi = 1 : 23
  
  rng('default')
  %   pause(ceil(rand*100));
  
  if ~exist([outdir sprintf('pmod_src_powcorr_allpow_f%d_v%d_processing.txt',ifoi,perm_v)])
    system(['touch ' outdir sprintf('pmod_src_powcorr_allpow_f%d_v%d_processing.txt',ifoi,perm_v)]);
  else
    continue
  end
  
  disp(sprintf('Processing freq %d ... Loading data ...',ifoi));
  
  load(['/home/tpfeffer/pconn/proc/conn/' sprintf('pconn_src_powcorr_allpow_f%d_v%d.mat',ifoi,v)])
  
  disp(sprintf('Processing freq %d ... Data loaded ...',ifoi));
  
  allpow        = single(allpow(:,:,:,1:para.nsubj));
  
  a     = jh_ranksum(cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:))));
  b     = a(logical(triu(ones(size(a)),1))); clear a
  pn    = [b>0 b<0];
  p     = 2*normcdf(-abs(b));
  cnt 	= [sum(p<para.thresh) sum(p<para.thresh)/length(p) sum(p(pn(:,1))<para.thresh) sum(p(pn(:,2))<para.thresh)]; clear a
  
  if para.perm
    
    cntperm = nan(para.nperm,4);
 
      allpow	= cat(3,squeeze(allpow(:,:,1,:)),squeeze(allpow(:,:,2,:)));
      
      for iperm = 1 : para.nperm
        
        clear a_perm
        disp(sprintf('processing f%d / p%d ...',ifoi,iperm))
        
        permidx           = randperm(para.nsubj*2);
        a_perm            = jh_ranksum(cat(3,squeeze(allpow(:,:,permidx(1:para.nsubj))),squeeze(allpow(:,:,permidx(para.nsubj+1:end)))));
        b_perm            = a_perm(logical(triu(ones(size(a_perm)),1))); clear a
        pn_perm           = [b_perm>0 b_perm<0];
        p_perm            = 2*normcdf(-abs(b_perm));
        cnt_perm(iperm,:)	= [sum(p_perm<para.thresh) sum(p_perm<para.thresh)/length(p_perm) sum(p_perm(pn_perm(:,1))<para.thresh) sum(p_perm(pn_perm(:,2))<para.thresh)];
        
      end
  end
  
  
  save([outdir sprintf('pmod_src_powcorr_allpow_f%d_v%d.mat',ifoi,1)],'cnt','cnt_perm','-v7.3');
     
end

error('!!!')
%% MULTIPLE COMPARISONS CORRECTION AND PLOTTING
% 
perm_v  = 1;
v = perm_v;
clear all_cnt acnt

plt = 1;
NFREQ = 21;

  
  clear corp thresh acnt all_cnt srt idx_D idx_R idx_D_max
  
  it = 0;
  for ifoi = 1 : 2 : NFREQ*2
    
    it = it + 1;
    try
    load([['~/pmod/proc/'] sprintf('pmod_src_powcorr_allpow_f%d_v%d.mat',it,perm_v)]);
    all_cnt(:,ifoi:ifoi+1)  = cnt_perm(:,3:4);
    acnt(:,ifoi:ifoi+1)     = cnt(:,3:4);
    catch me
    end
  end
  
%   ranks, where 1 = min, para.nperm = max
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
%   
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
  for ifoi = 2 : 2 : NFREQ*2
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/para.nperm;
    
  end
%   
  plot(log10(foi_range(1:NFREQ)),-log10(corp),'b','LineWidth',3); hold on; 
  
%   odd freqs: MS > HC
  cnt = 0;
  for ifoi = 1 : 2 : NFREQ*2-1
    cnt = cnt + 1;
    corp(cnt)  = sum(acnt(ifoi)<idx_D_max(:,ifoi))/para.nperm;
    
  end
%   
  plot(log10(foi_range(1:NFREQ)),-log10(corp),'r','LineWidth',3); hold on
  line(log10([foi_range(1) foi_range(NFREQ)]),-log10([.05 .05]),'LineStyle','--','color','k','LineWidth',2)
  
  set(gca,'XTick',[log10([2 4 8 16 32 64 128])],'XTickLabel',[2 4 8 16 32 64 128])
  set(gca,'YTick',[0 -log10(0.5) -log10(0.25) 1 2 3],'YTickLabel',[1 0.5 0.25 0.1 0.01 0.001])
  set(gca,'TickDir','out')

  title('Connectivity differences')
  xlabel('CARRIER FREQUENCY (HZ)')
  ylabel('P-VALUE')
set(gca,'YLim',[0 2.5]);
% saveas(gcf,sprintf([plotdir 'nc_src_all2all_stats_corr_v%d.fig'],v),'fig')

% 
% 
% 

