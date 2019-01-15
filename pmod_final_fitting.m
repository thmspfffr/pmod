%% FITTING
% pmod_final_fitting

clear
% %-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 0:0.05:1;
% Gains       = 0; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% %-------------------------------------------------------------------------
% VERSION 11: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %-------------------------------------------------------------------------
v           = 11;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0:0.05:1;
Gains       = 0; 
nTrials     = 1;
tmax        = 6500; % in units of tauE
EC = 0;
%-------------------------------------------------------------------------
% VERSION 1: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v           = 2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 0.87; % this is where correlation peaks 
% Gains       = [0 0.025:0.025:0.4 -0.025:-0.025:-0.4]; 
% nTrials     = 1;
% tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

v_sim = v;
% connectivity, AAL
v_conn =  1;

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
Gs = 1;

% ---------
% LOAD EMPIRICAL DFA (LCMV)
% ---------
% load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],v_dfa));
% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
load ~/M.mat

clear outp

addpath ~/pupmod/matlab
cleandat = pupmod_loadpowcorr(v_conn);

load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat

peakfreq_rest = m_res(1);
peakfreq_task = m_tsk(1);

load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

mask = logical(tril(ones(76,76),-1));

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

ifoi = 6;
fc_rest = squeeze(nanmean(nanmean(cleandat(:,:,:,1,1,ifoi),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(:,:,:,1,2,ifoi),3),6));

para = [];
para.transfer = 'to_bcn';
para.N = 90;

fc_rest = tp_match_aal(para,fc_rest);
fc_task = tp_match_aal(para,fc_task);

fc_rest = fc_rest(include_bcn,include_bcn);
fc_task = fc_task(include_bcn,include_bcn);

% transform indiv. subj. matrices to AAL BCN 
for isubj = 1 : 28
  fc_rest_indiv(:,isubj) = reshape(squeeze(nanmean(cleandat(include_bcn,include_bcn,isubj,1,1,ifoi),6)),[size(fc_rest,1)*size(fc_rest,1) 1]);
  fc_task_indiv(:,isubj) = reshape(squeeze(nanmean(cleandat(include_bcn,include_bcn,isubj,1,2,ifoi),6)),[size(fc_rest,1)*size(fc_rest,1) 1]);
end

fc_rest_indiv = fc_rest_indiv(mask,:);
fc_task_indiv = fc_task_indiv(mask,:);

k = [0 0];

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))

  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          %         igain
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))

          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;

          pars.dim = 2;
          
%           outp.fc_sim_tmp = out.FC;
          outp.fc_sim_env_tmp = out.FC_env;
                    
%           [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
          [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
          [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
          
%           [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
          [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv);
          %
          pars.dim = 1;
%           [~,~,outp.kds(iies,iiis,iG,igain)] = kstest2(fc_rest(mask),outp.fc_sim_env_tmp(mask));
   
%           outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
          outp.dist(iies,iiis,iG,igain) = 1-(outp.r_env_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
          outp.dist_indiv(iies,iiis,iG,igain) = 1-(squeeze(mean(outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain)))-(mean(fc_rest(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
          
%           outp.fc_sim_mean(iies,iiis,iG,igain)    = mean(outp.fc_sim_tmp(mask));
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          
%           if Gains(igain)>0     
%             k(1) = k(1) + 1;
%             all_FC_env_pos(:,:,iies,iiis,iG,k(1)) = single(out.FC_env); 
%           elseif Gains(igain)<0 
%             k(2) = k(2) + 1;
%             all_FC_env_neg(:,:,iies,iiis,iG,k(2)) = single(out.FC_env);
%           else
%             all_FC_env_zer(:,:,iies,iiis,iG) = single(out.FC_env);
%             save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_FC_zer_v%d.mat',v_sim),'all_FC_env_zer')
%             clear all_FC_env_zer
%           end
          fclose all;
          %
        end
      end
    end
  end

  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim),'outp')
%   save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_FC_pos_v%d.mat',v_sim),'all_FC_env_pos')
%   save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_FC_neg_v%d.mat',v_sim),'all_FC_env_neg')

%   % DELETE OUTOUT FILES
%   for iies = 1 : length(Ies)
%     iies
%     for iiis = 1 : length(Iis)
%       %     iiis
%       for iG = 1:length(Gg)
%         for igain = 1:length(Gains)
%           
%           fn =  sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v_sim);
%           delete(sprintf('%s_processing.txt',fn))
%           delete(sprintf('%s.mat',fn))
%           
%         end
%       end
%     end
%   end
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim))
% 
end

clear cleandat


error('!')

%% DETERMINE GLOBAL COUPLING PARAMETER FIRST

osc = osc1;
oscthresh = 0;

figure; set (gcf,'color','w');
subplot(2,2,1)
par = outp.dist_indiv;
% par(osc>oscthresh) = nan;
par1 = squeeze(nanmean(nanmean(par)));
plot(par1,'linewidth',2);
tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',1:0.05:1.2,'YTickLabel',1:0.05:1.2);
xlabel('Global coupling'); ylabel('Mean FC')
axis([1 21 0.90 1.2])
axis(gca,'square')

subplot(2,2,2)
par = outp.r_env_rest_corr;
% par(osc>oscthresh) = nan;
par1 = squeeze(nanmean(nanmean(par)))
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-0.02 0.15],'linestyle',':','color',[0.3 0.3 0.3])
tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0:0.1:0.4,'YTickLabel',0:0.1:0.4);
xlabel('Global coupling'); ylabel('Correlation (WC, MEG)')
axis([1 21 -0.02 0.15])
axis(gca,'square')

subplot(2,2,3)
par1 = squeeze(mean(mean(outp.fc_sim_env_mean)))
plot(par1,'linewidth',2);
tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0:0.1:0.4,'YTickLabel',0:0.1:0.4);
xlabel('Global coupling'); ylabel('Correlation (WC, MEG)')
axis([1 21 0 0.4])
[i,j]=max(par1);
line([j j],[-0.02 0.4],'linestyle',':','color',[0.3 0.3 0.3])
axis(gca,'square')


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_fitG_v%d.pdf',v_sim))


%%
%   load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',4))
osc = osc1;
oscthres = 0;
nancol = [0.97 0.97 0.97];

clear par ax
figure; set(gcf,'color','w')
% clear 

igain = 1;
iG = 1;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.peakfreq(:,:,iG,igain));
par(osc>oscthres)=nan;
imagescnan(par,[2 20],'NanColor',nancol)
title('Peak freq: Difference');

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = peakfreq_rest-squeeze(outp.peakfreq(:,:,iG,igain));
par(abs(par)<3)=1; par(abs(par)>3)=0; 
par(osc>oscthres)=0;
m1 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('Peak freq: Masked');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = -log10(squeeze(outp.p_env_rest_corr(:,:,iG,igain)));
par(osc>oscthres)=nan;
imagescnan(par,[0 3],'NanColor',nancol)
title('r(FC_{sim},FC_{MEG})');

% plot peak freq model
fdr_p = fdr1(reshape(squeeze(outp.p_env_rest_corr(:,:,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.01);
ax{4} = subplot(2,2,4); hold on
par = -log10(squeeze(outp.p_env_rest_corr(:,:,iG,igain)));
par(par<-log10(fdr_p))=0; par(par>=-log10(fdr_p))=1;
% par=par.*squeeze(outp.r_env_rest_corr(:,:,iG,igain))>0.15;
par(osc>oscthres)=0;
clust = bwlabel(par,8);
for iclust = 1 : max(clust(:))
  n_clust(iclust)=sum(clust(:)==iclust);
end
[~,k] = max(n_clust);
par = clust == k
m2 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('Peak freq: Masked');

[i,j]= find(m2);
idx = round([mean(i) mean(j)]);

for iax = 1 : 4
 scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
%   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  elseif iax ==4
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  end
  
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  if iax == 3
    c.TickLabels = {'0'; sprintf('%.3f',10^-max(c.Limits))};
  end
   axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
    c.Ticks = [min(c.Ticks) max(c.Ticks)];
   axis(ax{iax},'square')
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_fc_v%d.pdf',v_sim))

%% FIT TASK PARAMETERS (MEAN - see below for indiv subjects)
load redblue.mat
gain = 1; 
mask = logical(tril(ones(76,76),-1));
fc_tvr_sim = 100*(outp.fc_sim_mean(:,:,:,gain)-outp.fc_sim_mean(idx(1),idx(2),:,gain))./outp.fc_sim_mean(idx(1),idx(2),:,gain)
fc_tvr_exp = 100*(mean(fc_task(mask))-mean(fc_rest(mask)))./mean(fc_rest(mask));

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = fc_tvr_sim;
par(par<(fc_tvr_exp-10))=0; par(par>0)=0; par(par<(fc_tvr_exp+9))=1; 
par(1:idx(1),1:idx(2))=nan
% par=par.*squeeze(outp.r_env_rest_corr(:,:,iG,igain))>0.15;
par(osc>oscthres)=0;
m2 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('FC decr. Task (+/- 10%)');
colormap(redblue)

scatter(gca,idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')

[idx_task(:,1),idx_task(:,2)]=find(par>0)

scatter(gca,idx_task(:,2),idx_task(:,1),20,'markerfacecolor','r','markeredgecolor','k')

iax = 4
colormap(ax{iax},redblue)
ylabel(ax{iax},'Excitatory input');
set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
xlabel(ax{iax},'Inhibitory input'); tp_editplots;
axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
%   c.Ticks = [min(c.Ticks) max(c.Ticks)];
axis(ax{iax},'square')

clear par
par.rest = [Ies(idx(1)) Iis(idx(2))];
par.task = [Ies(idx_task(1,1)) Iis(idx_task(1,2))];
par.task_alt = [Ies(idx_task(2:end,1)); Iis(idx_task(2:end,2))]';
par.descr = 'First entry: Excitation (Ies), Second entry: Inhibition (Iis)';

save(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',v_sim),'par')


%% INDIV SUBJ
iG = 1;

% osc=zeros(301,401)
oscthresh = 0
h=figure; set(gcf,'color','w')

width = 10;
bif_mask = zeros(size(osc));
for i = 1 : size(bif_mask,1)
  if isempty(find(osc(i,:)>oscthresh,1,'last'))
    idx=0;
  else
    idx = find(osc(i,:)>oscthresh,1,'last')
    if idx-width < 1

      bif_mask(i,idx-(width+(idx-(width+1))):idx+width)=1;
    else
      bif_mask(i,idx-width:idx+width)=1;
    end
  end
end
bif_mask = bif_mask(1:size(osc,1),1:size(osc,1))
clear idx
ax{1}=subplot(6,5,1); hold on

imagescnan(bif_mask,[-1 1])

colormap(redblue)
set(gca,'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(gca,'YTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots

for isubj = 1 : 28
  % plot lambda
  ax{1} = subplot(6,5,isubj+1); hold on
  p_corr = fdr1(reshape(squeeze(outp.p_env_rest_indiv_corr(isubj,:,:,iG,1)),[size(osc,1)*size(osc,2) 1]),0.05);
  par = squeeze(outp.p_env_rest_indiv_corr(isubj,:,:,iG,1))<p_corr; par=double(par);

  par(osc>0)=0;
  par(~bif_mask)=0;
  
  bw = bwlabel(par,8);
  bw(osc>0)=0;
  cnt = [];
  for iclust = 1 : max(bw(:))
    cnt(iclust) = sum(bw(:)==iclust);
  end
  
  if isempty(cnt)
    continue
  end
  
  
  [~,k(isubj)]=max(cnt);
%   
%   if sum(cnt==max(cnt)>1)
%     for iclust = 1 : sum(cnt==max(cnt))
%       find(bw==
%     end
%   end
%   

  par(bw==k(isubj))=1;
  par(bw~=k(isubj))=nan;
  par(osc>0)=nan;

  imagescnan(par,[-1 1])
  colormap(redblue)

  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))

  [i,kk]=find(par>0);
  idx(isubj,:) = [round(mean(i)) round(mean(kk))];
  scatter(gca,idx(isubj,2),idx(isubj,1),20,'markerfacecolor','w','markeredgecolor','k');
  tp_editplots
end

% set(h,'Position',[50 50 1200 900])
% set(h,'Renderer','Painters')
% print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv.pdf'))

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_v%d.mat',v_sim),'idx')

% close all
%%
figure; set(gcf,'color','white')

ax{1} = subplot(2,2,1); hold on

imagescnan(bif_mask,[-1 1])

colormap(redblue)
set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
set(gca,'YTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
tp_editplots
title('Bifurcation Mask')
set(ax{1},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{1},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
tp_editplots
xlabel(ax{1},'Inhibitory input')
ylabel(ax{1},'Excitatory input');
axis(ax{1},'square'); axis tight
colorbar
axis(gca,[1 size(osc1,1) 1 size(osc1,2)])


ax{2} = subplot(2,2,4); hold on
scatter(gca,idx(:,2),idx(:,1),10,'markerfacecolor','w','markeredgecolor','k');
axis(gca,[1 size(osc1,1) 1 size(osc1,2)])

set(ax{2},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{2},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
tp_editplots
xlabel(ax{2},'Inhibitory input')
ylabel(ax{2},'Excitatory input');
axis(ax{2},'square')
title('Indiv subjects')
colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_sum.pdf'))


%% GET INDIV TASK PARAMETERS
% iG = 1;
clear idx_task
indiv_change_prc = 100*(mean(fc_task_indiv,1)-mean(fc_rest_indiv,1))./mean(fc_rest_indiv,1);

for isubj = 1:28
  
  if isnan(idx(isubj,1)) || idx(isubj,1)==0
    continue
  end
  fc_tvr_sim = 100*(outp.fc_sim_mean(:,:,1,1)-outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),1,1))./outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),1,1)

  d = indiv_change_prc(isubj)-fc_tvr_sim;
  d(osc>0) = nan;

  mask = zeros(size(bif_mask,1),size(bif_mask,1));
  i = 0;
  while 1
    if idx(isubj,1)+i <= 31
      mask(idx(isubj,1)+i,idx(isubj,2)+i) = 1;
    else
      break
    end
    mask(idx(isubj,1)+i,1:idx(isubj,2)+i)=1;
    i = i + 1;
  end
  mask(1:idx(isubj,1),1:idx(isubj,2)) = 1;
  mask(1:idx(isubj,1),:) = 1; mask = ~mask;
  
  mask = mask&bif_mask;

thresh = 1;
% d(~logical(bif_mask)) = Inf;
d(~logical(mask)) = Inf;
% d(1:idx(isubj,1),1:idx(isubj,2)) = Inf;

siz = [121 121];
warning('Change size!!!!!\n!!!!\n!!!!')

while 1
  thresh
  clust = reshape(bwlabel(abs(d)<thresh,8),[siz(1)*siz(2) 1]);
  if max(clust)<1
    thresh = thresh+1;
    continue
  else
    clear cnt
    for iclust = 1 : max(clust)
      cnt(iclust) = sum(clust==iclust);
    end
    if any(cnt>7)
      break
    else
      thresh = thresh+1;
      continue
    end
  end
end

[i,j] = find(bwlabel(abs(d)<thresh)==find(cnt>1,1,'first'))
idx_task(isubj,:) = round([mean(i) mean(j)]);

end


figure; set(gcf,'color','w')

ax{1}= subplot(2,2,1); hold on

imagesc(zeros(109,121))
scatter(ax{1},idx(:,2),idx(:,1),20,'markerfacecolor','w','markeredgecolor','k');
scatter(ax{1},idx_task(:,2),idx_task(:,1),20,'markerfacecolor','y','markeredgecolor','k');
set(ax{1},'YTick',1:20:length(Ies),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{1},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
axis([1 109 1 121])
tp_editplots; colormap(redblue); colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_task_sum.pdf'))

clear par

indiv_idx.rest = idx;
indiv_idx.task = idx_task;

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim),'indiv_idx')

%% Number of altered correlations






%% FIT GAIN FOR EACH SUBJECT
% do for number of altered correaltions
%    

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v_conn));
cleandat = cleandat(include_bcn,include_bcn,:,:,:,:);

%     
% %
k = 0;
GAINS = [18:33 1:17];
EI =  [-16:1:16];
for gain = 1:length(GAINS)
  igain = GAINS(gain);
  for ei = 1 : length(EI)
    iei = EI(ei);
    for isubj = 1:28
      
      fc_frac(1) = sum(sum(triu(cleandat(:,:,isubj,2,1)-cleandat(:,:,isubj,1,1),1)>0))/((76*76-76)/2);
      fc_frac(2) = sum(sum(triu(cleandat(:,:,isubj,2,2)-cleandat(:,:,isubj,1,2),1)>0))/((76*76-76)/2);
       
      ie = idx(isubj,1);
      ii = idx(isubj,2);
      load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',ie,ii,1,igain,v_sim))
%       fc1 = 
%       fc_gain(isubj,1) = sum(sum(triu(outp.FC-outp.FC>0,1))/((76*76-76)/2);
%       fc_gain(isubj,2) = sum(sum(triu(outp.FC-outp.FC>0,1))/((76*76-76)/2);
            
      err(gain,ei,isubj) = (fc_gain(isubj,2)-fc_rest(isubj,2)).^2;
%       (fc_gain(isubj,1)-fc_rest(isubj,1)).^2 + 
    end
  end
end


%% LOAD PARAMETERS AND PLOT GAIN VS NO GAIN
gain = 3;

figure; set(gcf,'color','w')

load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim))

ax{1}= subplot(2,2,1); hold on

imagesc(zeros(size(osc1,1),size(osc1,2)))
scatter(ax{1},indiv_idx.rest(:,2),indiv_idx.rest(:,1),20,'markerfacecolor','w','markeredgecolor','k');
scatter(ax{1},indiv_idx.task(:,2),indiv_idx.task(:,1),20,'markerfacecolor','y','markeredgecolor','k');
set(ax{1},'YTick',1:20:length(Ies),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{1},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
axis([1 size(osc1,1) 1 size(osc1,2)])
tp_editplots; colormap(redblue); colorbar

ax{2}= subplot(2,2,2); hold on
par = outp.fc_sim_env_mean(:,:,:,gain)-outp.fc_sim_env_mean(:,:,:,1);
par(osc1>oscthresh)=nan;
imagescnan(par,[-0.02 0.02])
scatter(ax{2},indiv_idx.rest(:,2),indiv_idx.rest(:,1),20,'markerfacecolor','w','markeredgecolor','k');
scatter(ax{2},indiv_idx.task(:,2),indiv_idx.task(:,1),20,'markerfacecolor','y','markeredgecolor','k');
set(ax{2},'YTick',1:20:length(Ies),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{2},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
xlabel(ax{2},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
axis(ax{2},'square')
axis([1 size(osc1,1) 1 size(osc1,2)])
tp_editplots; colormap(redblue); colorbar


% for isubj = 1 : 28
%   fc_sim(:,:,isubj,1) = squeeze(outp.fc_env(:,:,indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,gain))-squeeze(outp.fc_env(:,:,indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,1));
%   fc_sim(:,:,isubj,2) = squeeze(outp.fc_env(:,:,indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,gain))-squeeze(outp.fc_env(:,:,indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,1));
%   fc_sim_mean(isubj,1) = squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,gain));
% end
% 
% for gain = 1 : 21
% sim1 = squeeze(outp.fc_sim_env_mean(indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,1));
% emp1 = nanmean(nanmean(cleandat(include_bcn,include_bcn,isubj,2,1,6)))-nanmean(nanmean(cleandat(include_bcn,include_bcn,isubj,1,1,6)));
% err(gain,1) = (sim1-emp1).^2
% 
% sim2 = squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,1));
% emp2 = nanmean(nanmean(cleandat(include_bcn,include_bcn,isubj,2,2,6)))-nanmean(nanmean(cleandat(include_bcn,include_bcn,isubj,2,1,6)));
% err(gain,2) = (sim2-emp2).^2
% 
% % err(gain,2) = (squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,1)).^2)
% end
% 
% [h,~,~,s]=ttest(fc_sim(:,:,:,1),zeros(76,76,28),'dim',3)
% altered_corr_p(:,1) = sum(sum(triu(h>0&s.tstat>0,1)))/((76*76-76)/2);
% altered_corr_n(:,1) = sum(sum(triu(h>0&s.tstat<0,1)))/((76*76-76)/2);
% 
% [h,~,~,s]=ttest(fc_sim(:,:,:,2),zeros(76,76,28),'dim',3)
% altered_corr_p(:,2) = sum(sum(triu(h>0&s.tstat>0,1)))/((76*76-76)/2);
% altered_corr_n(:,2) = sum(sum(triu(h>0&s.tstat<0,1)))/((76*76-76)/2);
% close
% COMPUTE NUMBER OF ALTERED CORRELATIONS!!!
 for isubj = 1 : 28
  fc_sim_mean(isubj,1) = squeeze(outp.fc_sim_env_mean(indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.rest(isubj,1),indiv_idx.rest(isubj,2),:,1));
  fc_sim_mean(isubj,2) = squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.task(isubj,1),indiv_idx.task(isubj,2),:,1));
 end

 [~,p]=ttest(fc_sim_mean,zeros(28,2),'dim',1)
 
 
 [r(1),p(1)]=corr(fc_sim_mean(:,1),fc_rest(:,1))
 [r(2),p(2)]=corr(fc_sim_mean(:,2),fc_rest(:,2))

% print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_task_sum.pdf'))

% clear par




%% BOLD

figure; set(gcf,'color','w')
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.r_rest_BOLD_corr(:,:,iG,igain));
par(osc>oscthres)=nan;
imagescnan(par,[0 0.1],'NanColor',nancol)
title('r(FC_{sim},FC_{MEG})');

% plot peak freq model
fdr_p = fdr1(reshape(squeeze(outp.p_rest_BOLD_corr(:,:,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.05);
ax{2} = subplot(2,2,2); hold on
par = -log10(squeeze(outp.p_rest_BOLD_corr(:,:,iG,igain)));
par(par<-log(fdr_p))=0; par(par>=-log(fdr_p))=1;
% par=par.*squeeze(outp.r_env_rest_corr(:,:,iG,igain))>0.15;
par(osc>oscthres)=0;
m2 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('Peak freq: Masked');

[i,j]= find(m2);
idx = round([mean(i) mean(j)]);

for iax = 1 : 2 
%  scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
%   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  elseif iax == 2
    colormap(ax{iax},redblue)
    set(ax{iax},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
    set(ax{iax},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
%   elseif iax == 3
%     colormap(ax{iax},plasma)
%     xlabel(ax{iax},'Inhibitory input')
%     ylabel(ax{iax},'Excitatory input');
%     set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
%     set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
%   elseif iax ==4
%     colormap(ax{iax},redblue)
%     xlabel(ax{iax},'Inhibitory input')
%     set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
%     set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  if iax == 3
    c.TickLabels = {'0'; sprintf('%.3f',10^-max(c.Limits))};
  end
   axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
   axis(ax{iax},'square')
  
end
