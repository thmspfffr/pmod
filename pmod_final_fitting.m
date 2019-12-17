%% FITTING
% pmod_final_fitting

clear
% %-------------------------------------------------------------------------
% VERSION 1:
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.05:-1;
% Iis         = -5:0.05:-2;
% Gg          = 0:0.05:2;
% Gains       = [0 0.45]; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
% dt          = 0.01;
%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 1.4;
Gains       = [0:0.05:0.6 -0.2:0.05:-0.05 0.65:0.05:1]; 
nTrials     = 1;
tmax        = 6500;  % in units of tauE
EC          = 0;
dt          = 0.01;
%-------------------------------------------------------------------------

% load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_parameters_v%d.mat',v,v))

% Gains = [0 0.025:0.025:0.4 -0.025:-0.025:-0.3 ]
v_sim = v;
% connectivity, AAL
v_conn =  25;
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

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
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pupmod/matlab
cleandat = pupmod_loadpowcorr(v_conn,SUBJLIST,1);

mask = logical(tril(ones(76,76),-1));

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

% LOAD SC MATRIX
load ~/sc90.mat
SC = SC(include_bcn,include_bcn);

% ifoi = 1;
fc_rest = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,1,2),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,2,2),3),6));

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

fc_rest = tp_match_aal(para,fc_rest);
fc_task = tp_match_aal(para,fc_task);

fc_rest = fc_rest(include_bcn,include_bcn);
fc_task = fc_task(include_bcn,include_bcn);

% transform indiv. subj. matrices to AAL BCN
for isubj = 1 : 28
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,1,:),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  [corrwithfc_rest(isubj), p_corrwithfc_rest(isubj)]  = corr(tmp(mask),SC(mask));
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,2,:),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  fc_task_indiv(:,isubj) = tmp(mask);
  [corrwithfc_task(isubj), p_corrwithfc_task(isubj)] = corr(tmp(mask),SC(mask));
end

k = [0 0];
%%
if ~exist(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  
  for iies = 1 : length(Ies)
    iies
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        for iiis = 1 : length(Iis)
          
          load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,iies,iiis,iG,igain,v_sim))
%           out.peakfreq
          
          outp.fc_sim_env_tmp = out.rc_wl;
          %             outp.fc_sim_env_ad = out(iiis).FC_env_ad;
          
          [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
          [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
          
          %             [outp.r_env_rest_corr_ad(iies,iiis,iG,igain), outp.p_env_rest_corr_ad(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_ad(mask),fc_rest(mask));
          
          
          [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv);
          [outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_task_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task_indiv);
          
          outp.dist_rest(iies,iiis,iG,igain) = 1-(outp.r_env_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
          outp.dist_task(iies,iiis,iG,igain) = 1-(outp.r_env_task_corr(iies,iiis,iG,igain)-(mean(fc_task(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
          
          outp.dist_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
          outp.dist_task_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_task_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
          
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          %             outp.fc_sim_env_mean_ad(iies,iiis,iG,igain) = mean( outp.fc_sim_env_ad(mask));
          
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          if isempty( out.peakfreq)
            
             out.peakfreq = nan;
          end
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          outp.alphapow(iies,iiis,iG,igain) = mean(out.alphapow);
        end
        fclose all;
        
      end
    end
  end
  
  save(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim),'outp')

else
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  load(sprintf('~/pmod/proc/detosc/v%d/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim,v_sim))
end
% osc1=osc1(1:121,1:121,:,:);
clear cleandat

error('!')

%%
lim = 121;

oscthres = 0;
nancol = [0.97 0.97 0.97];


clear par ax
figure; set(gcf,'color','w')

igain = 1;
iG = 1;
osc = osc1(1:lim,1:lim,iG);

% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = outp.fc_sim_env_mean(1:lim,1:lim,iG,igain);
par(osc>oscthres)=nan;
b=imagesc(par,[0 0.02]); set(b,'AlphaData',~isnan(par))
title('Mean FC (Envelopes)');

% % plot corr fc_model with fc_meg
% ax{2} = subplot(2,2,2); hold on
% par = outp.r_env_rest_corr(1:lim,1:lim,iG,igain);
% par(osc>oscthres)=nan;
% b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isnan(par))
% title('r(FC_{sim},FC_{MEG})');

% plot corr fc_model with fc_meg
ax{3} = subplot(2,2,2); hold on
par = outp.r_env_rest_corr(1:lim,1:lim,iG,igain);
par(osc>oscthres)=nan;
b=imagesc(par,[0 0.4]); set(b,'AlphaData',~isnan(par))
title('r(FC_{sim},FC_{MEG})');

k=(par>prctile(par(:),99));
clust = bwlabel(k,8);
for iclust = 1 : max(clust(:))
  n_clust(iclust)=sum(clust(:)==iclust);
end

[~,k] = max(n_clust);
par = clust == k
m2 = par>0;
[i,j]= find(m2);
idx = round([mean(i) mean(j)]);

scatter(idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')

par = outp.dist_rest(:,:,iG,igain); par(osc>oscthres)=nan;
prc = prctile(reshape(par,[lim*lim 1 1]),1);

ax{4} = subplot(2,2,4); hold on
par = par<prc;
par(osc>oscthres)=0;
clust = bwlabel(par,8);
for iclust = 1 : max(clust(:))
  n_clust(iclust)=sum(clust(:)==iclust);
end

[~,k] = max(n_clust);
par = clust == k
m2 = par>0;
par = double(par); par(par<1)=NaN;
b=imagesc(par,[0 1]); set(b,'AlphaData',~isnan(par))
title('Correlation: p-values (log10)');

[i,j]= find(m2);
idx_avg = round([mean(i) mean(j)]);

for iax = [1 3 4]
  scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax ==4
    colormap(ax{iax},[nancol; 1 0 0])
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  end
  
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  if iax == 3
    c.TickLabels = {'0';'1'};
  end
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
  axis(ax{iax},'square')
  
end

set(gcf,'Renderer','Painters')
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_fc_v%d.eps',v_sim))

%% FIT INDIVIDUAL SUBJECTS
% --------------------------------------
% --------------------------------------

clear idx par ax cnt idx_rest par

iG         = 1; 
igain      = 1;
osc        = osc1(1:lim,1:lim,iG,igain);
oscthresh  = 0;
prc_thresh = 99;


h=figure; set(gcf,'color','w')

% Loop through all subjects
for isubj = 1 : 28
  clear par bw
  
    clear k
    par                 = squeeze(outp.r_env_rest_indiv_corr(isubj,1:lim,1:lim,iG,igain));    
    par(osc>oscthresh)  = nan;
    binmap              = par>prctile(par(:),prc_thresh);
    bw                  = bwlabel(binmap,8);
    cnt                 = [];
    
    for iclust = 1 : max(bw(:))
      cnt(iclust) = sum(bw(:)==iclust);
    end
    
    % if no cluster(s), omit subject from further analysis
    if isempty(cnt)
      idx_rest.inh(isubj) = nan;
      idx_rest.exc(isubj) = nan;
      continue
    end
    
    % find single largest cluster
    [~,k(isubj)]       = max(cnt);
    par(bw==k(isubj))  = 1;
    par(bw~=k(isubj))  = 0;
    r = squeeze(outp.r_env_rest_indiv_corr(isubj,1:lim,1:lim,iG,igain));
    [i,kk]=find(par==1);
    idx(isubj,:) = round([mean(i) mean(kk)]);
%     idx(isubj,:) = [round(sum(r(logical(par)).*i)/sum(r(logical(par)))) round(sum(r(par).*kk)/sum(r(par)))];
   % PLOT RESULT
  ax{1}   = subplot(6,5,isubj); hold on
  par = squeeze(outp.dist_rest_indiv(isubj,1:lim,1:lim,iG,igain));
  par(osc>oscthresh)  = nan;
  b=imagesc(par,[0.6 1.4]); colormap(plasma)
  set(b,'AlphaData',~isnan(par))
  scatter(gca,idx(isubj,2),idx(isubj,1),20,'markerfacecolor','w','markeredgecolor','k');
  tp_editplots; axis square
  
  set(h,'Renderer','Painters')
  set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
  set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
  idx_rest.inh(isubj) = idx(isubj,2);
  idx_rest.exc(isubj) = idx(isubj,1);
end

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',v_sim),'idx_rest')

clear idx

%% GET INDIV TASK PARAMETERS

clear idx_task idx
task_effect = 6;

for isubj = 1:28
idx_task.inh(isubj)= idx_rest.inh(isubj)+task_effect;
  idx_task.exc(isubj)= idx_rest.exc(isubj)+task_effect;

  
end

indiv_idx.rest = idx_rest;
indiv_idx.task = idx_task;
%%
clear d
[sorted_gains,ii]=sort(Gains);
for igain = ii
  
  for isubj = 1 : 28
    isubj
  d(1,igain,isubj) = squeeze(outp.fc_sim_env_mean(round(idx_rest.exc(isubj)),round(idx_rest.inh(isubj)),:,igain)-outp.fc_sim_env_mean(round(idx_rest.exc(isubj)),round(idx_rest.inh(isubj)),:,1));
  d(2,igain,isubj) = squeeze(outp.fc_sim_env_mean(round(idx_task.exc(isubj)),round(idx_task.inh(isubj)),:,igain)-outp.fc_sim_env_mean(round(idx_task.exc(isubj)),round(idx_task.inh(isubj)),:,1));
  end
end

figure; set(gcf,'color','w');hold on
plot(sorted_gains,nanmean(d(1,:,:),3))
plot(sorted_gains,nanmean(d(2,:,:),3))
line([-0.3 1.1],[0 0],'linestyle','--','color',[.5 .5 .5])
axis([-0.3 1.1 -0.8 0.8]); tp_editplots
%% PLOT REST AND TASK

figure; set(gcf,'color','white')
load redblue.mat

ax{1} = subplot(2,2,1); hold on
imagesc(zeros(lim,lim))
set(gca,'ydir','normal')
scatter(gca,idx_rest.inh,idx_rest.exc,25,'markerfacecolor','r','markeredgecolor','k');
axis(ax{1},[1 lim 1 lim])
set(ax{1},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{1},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
title('Indiv subjects')
axis(ax{1},'square')
tp_editplots; colormap([0.9 0.9 0.9])

ax{2}= subplot(2,2,2); hold on

imagesc(zeros(lim,lim))
set(gca,'ydir','normal')
scatter(gca,idx_task.inh,idx_task.exc,25,'markerfacecolor','y','markeredgecolor','k');
axis(ax{2},[1 lim 1 lim])
set(ax{2},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{2},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{2},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
title('Indiv subjects')
axis(ax{2},'square')
tp_editplots

ax{3}= subplot(2,2,3); hold on
imagesc(zeros(lim,lim))
scatter(gca,round(nanmean(indiv_idx.rest.inh)),round(nanmean(indiv_idx.rest.exc)),35,'markerfacecolor','r','markeredgecolor','k');
scatter(gca,round(nanmean(indiv_idx.task.inh)),round(nanmean(indiv_idx.task.exc)),35,'markerfacecolor','y','markeredgecolor','k');
axis(ax{3},[1 lim 1 lim])
set(ax{3},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{3},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{3},'Input to inhibitory population'); ylabel(ax{1},'Input to excitatory population');
title('Indiv subjects')
axis(ax{3},'square')
tp_editplots

% set(gca,'ydir','normal')
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_taskandrest_v%d.eps',v_sim))

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim),'indiv_idx')

%% PLOT HIGH GAIN VS BASELINE GAIN
Ies_ad_rest = -2.85;% -1.85]; 
Iis_ad_rest = -3.50;% -2.2];

Ies_ad_task = -1.85; 
Iis_ad_task = -2.2;

for igain = 1:25
  
  
  osc=osc1(:,:,:,1);
  par = outp.fc_sim_env_mean(:,:,:,igain)-outp.fc_sim_env_mean(:,:,:,1);
  par(osc>oscthres)=0;

  figure; set(gcf,'color','w'); hold on
  subplot(2,2,1); hold on
  b=imagesc(Iis,Ies,par,[-0.03 0.03]); set(gca,'ydir','normal');  set(b,'AlphaData',~isnan(par))
  hold on
  colormap(redblue); axis square
  scatter(gca,Iis(round(nanmean(indiv_idx.rest.inh))),Ies(round(nanmean(indiv_idx.rest.exc))),35,'markerfacecolor','w','markeredgecolor','k');
  scatter(gca,Iis(round(nanmean(indiv_idx.task.inh))),Ies(round(nanmean(indiv_idx.task.exc))),35,'markerfacecolor','y','markeredgecolor','k');
%   scatter(gca,Iis(idx_rest.inh),Ies(idx_rest.exc),10,'markerfacecolor','r','markeredgecolor','k');

  set(gca,'YTick',Ies(1:40:end),'YTickLabels',num2cell(Ies(1:40:end)))
  set(gca,'XTick',Iis(1:40:end),'XTickLabels',num2cell(Iis(1:40:end)))
  axis([-5 -2 -4 -1]); xlabel('Input to I');  ylabel('Input to E'); 
  tp_editplots
  
  subplot(2,2,2); hold on
  b=imagesc(Iis,Ies,par,[-0.03 0.03]); set(gca,'ydir','normal');  set(b,'AlphaData',~isnan(par))
  hold on
  colormap(redblue); axis square
  scatter(gca,Iis(round(nanmean(indiv_idx.rest.inh))),Ies(round(nanmean(indiv_idx.rest.exc))),35,'markerfacecolor','w','markeredgecolor','k');
  scatter(gca,Iis(round(nanmean(indiv_idx.task.inh))),Ies(round(nanmean(indiv_idx.task.exc))),35,'markerfacecolor','y','markeredgecolor','k');
%   scatter(gca,Iis(idx_task.inh),Ies(idx_task.exc),10,'markerfacecolor','y','markeredgecolor','k');
  tp_editplots
  set(gca,'YTick',Ies(1:40:end),'YTickLabels',num2cell(Ies(1:40:end)))
  set(gca,'XTick',Iis(1:40:end),'XTickLabels',num2cell(Iis(1:40:end)))
  axis([-5 -2 -4 -1]); xlabel('Input to I');  ylabel('Input to E'); 
  
  print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_highvslowgain_gain%d_v%d.eps',igain,v_sim))
end
%% LOAD FRACTIION OF ALTERED CORRELATIONS FROM DATA
v = 25;
para.nfreq = 1:3;
para.alpha = 0.05;
SUBJLIST1  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));
para.transfer = 'to_bcn'
if ~exist('emp','var')
  if ~exist('cleandat','var')
    cleandat = pupmod_loadpowcorr(v,SUBJLIST1,1);
  end
  % transform indiv. subj. matrices to AAL BCN
  for isubj = 1 : 28
    isubj
    for im = 1 : 3
      for icont = 1 : 2
        for ifreq = 1 : max(para.nfreq)
          tmp = squeeze(cleandat(:,:,isubj,im,icont,ifreq));
          tmp = tp_match_aal(para,tmp);
          fc_indiv(:,:,isubj,im,icont,ifreq) = tmp(include_bcn,include_bcn);
      
        end
      end
    end
  end
  emp = pupmod_compute_altered_correlations(fc_indiv,para);
end


%% PLOT ALTERED CORRELATIONS AS FUNCTION OF GAIN
clear FC_simrest FC_simtask pp_pos_atx_sim pp_neg_atx_sim

SUBJLIST = find(idx_rest.inh>15);

mask = logical(tril(ones(76,76),-1));
[~,idx_gain]=sort(Gains)

idx_gain1 = idx_gain(1:1:end);

pp_pos_atx_sim = nan(length(idx_gain1),2);
pp_neg_atx_sim = nan(length(idx_gain1),2);

for igain = 1:length(idx_gain(1:1:end))
  
  if Gains(idx_gain1(igain)) == 0
    pp_pos_atx_sim(igain,:) = [nan nan];
    pp_neg_atx_sim(igain,:) = [nan nan];
    continue
  end
  
  gain = idx_gain1(igain)
  
  for isubj = 1:length(SUBJLIST)
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simrest(:,:,isubj,1) = out.rc_wl; 
    pf_rest(isubj,1) = out.peakfreq; clear out
    load(sprintf('~/pmod/proc/numerical/v%d//pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simrest(:,:,isubj,2) = out.rc_wl; 
    pf_rest(isubj,2) = out.peakfreq; clear out
    load(sprintf('~/pmod/proc/numerical/v%d//pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simtask(:,:,isubj,1) = out.rc_wl; 
    pf_task(isubj,1) = out.peakfreq; clear out
    load(sprintf('~/pmod/proc/numerical/v%d//pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simtask(:,:,isubj,2) = out.rc_wl; 
    pf_task(isubj,2) = out.peakfreq; clear out
  end
  
  fc_sim_rest(igain,:) = mean(squeeze(mean(mean(FC_simrest),2)));
  fc_sim_task(igain,:) = mean(squeeze(mean(mean(FC_simtask),2)));
  
  [h,~,~,s]=ttest(atanh(FC_simrest(:,:,:,2)),atanh(FC_simrest(:,:,:,1)),'dim',3);
  pp_pos_atx_sim(igain,1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
  pp_neg_atx_sim(igain,1) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
  
  [h,~,~,s]=ttest(atanh(FC_simtask(:,:,:,2)),atanh(FC_simtask(:,:,:,1)),'dim',3);
  pp_pos_atx_sim(igain,2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
  pp_neg_atx_sim(igain,2) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
  %
  %   [h,~,~,s]=ttest(atanh(FC_simtask(:,:,:,2)),atanh(FC_simrest(:,:,:,2)),'dim',3);
  %   p_pos_tvr_sim(igain) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
  %   p_neg_tvr_sim(igain) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
  %
  %   diff_tvr(igain) = sum(sum(triu(nanmean(FC_simtask(:,:,:,2),3)-nanmean(FC_simrest(:,:,:,2),3),1)))./sum(mask(:));
  
end
%
% [h,~,~,s]=ttest(atanh(FC_simtask(:,:,:,1)),atanh(FC_simrest(:,:,:,1)),'dim',3);
% p_pos_tvr_sim(Gains(idx_gain)==0) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
% p_neg_tvr_sim(Gains(idx_gain)==0) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
%%
% PLOT
% -----------------------
figure; set(gcf,'color','w');
subplot(4,3,1); hold on
plot(pp_pos_atx_sim(:,1),'k');
plot(pp_neg_atx_sim(:,1),'color',[.5 .5 .5],'linestyle','-');
plot(length(idx_gain1)+2,nanmean(emp.n_p_atx(:,1)),'o','markerfacecolor','k')
plot(length(idx_gain1)+2,nanmean(emp.n_n_atx(:,1)),'o','markerfacecolor',[.5 .5 .5])
line([0 length(idx_gain1)+2],[nanmean(emp.n_p_atx(:,1)) nanmean(emp.n_p_atx(:,1))],'color','k','linestyle','--')
line([0 length(idx_gain1)+2],[nanmean(emp.n_n_atx(:,1)) nanmean(emp.n_n_atx(:,2))],'color',[.5 .5 .5],'linestyle','--')

axis([0 length(idx_gain1)+2  -0.05 1]);
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[0 0.50  1],'yTickLabels',num2cell([0 50 100]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

subplot(4,3,3); hold on
plot(pp_pos_atx_sim(:,2),'k');
plot(pp_neg_atx_sim(:,2),'color',[.5 .5 .5],'linestyle','-');
plot(length(idx_gain1)+2,nanmean(emp.n_p_atx(:,2)),'o','markerfacecolor','k')
plot(length(idx_gain1)+2,nanmean(emp.n_n_atx(:,2)),'o','markerfacecolor',[.5 .5 .5])
axis([0 length(idx_gain1)+2 -0.05 1]);
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain1)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[0  0.50  1],'yTickLabels',num2cell([0  50  100]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots
line([0 length(idx_gain1)+2],[nanmean(emp.n_p_atx(:,2)) nanmean(emp.n_p_atx(:,2))],'color','k','linestyle','--')
line([0 length(idx_gain1)+2],[nanmean(emp.n_n_atx(:,2)) nanmean(emp.n_n_atx(:,2))],'color',[.5 .5 .5],'linestyle','--')

% PLOT CONTEXT DEPENDENCE
subplot(4,3,5); hold on
plot(pp_pos_atx_sim(:,1)-pp_pos_atx_sim(:,2),'k');
plot(pp_neg_atx_sim(:,1)-pp_neg_atx_sim(:,2),'color',[.5 .5 .5],'linestyle','-');
plot(length(idx_gain1)+2,nanmean(emp.n_n_atx(:,1))-nanmean(emp.n_n_atx(:,2)),'o','markerfacecolor','b')
plot(length(idx_gain1)+2,nanmean(emp.n_p_atx(:,1))-nanmean(emp.n_p_atx(:,2)),'o','markerfacecolor','r')
axis([0 length(idx_gain1)+2  -1 1]);
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain1)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[-1 -0.5 0 0.5 1],'yTickLabels',num2cell([-100 -50 0 50 100]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots
line([0 length(idx_gain1)+2],[nanmean(emp.n_n_atx(:,1))-nanmean(emp.n_n_atx(:,2)) nanmean(emp.n_n_atx(:,1))-nanmean(emp.n_n_atx(:,2))],'color','b','linestyle',':')
line([0 length(idx_gain1)+2],[nanmean(emp.n_p_atx(:,1))-nanmean(emp.n_p_atx(:,2)) nanmean(emp.n_p_atx(:,1))-nanmean(emp.n_p_atx(:,2))],'color','r','linestyle',':')
line([find(round(Gains(idx_gain)*100)==55) find(round(Gains(idx_gain)*100)==55)],[-1 0],'color',[0.5 0.5 0.5],'linestyle','-')


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_fraction_as_gain_v%d.pdf',v_sim))
%% CORRELATE DIFFERENCE MATRICES
fc_task = squeeze(nanmean(cleandat(:,:,:,:,2,:),6));%-nanmean(cleandat(:,:,:,1,2,:),6));
v = 2;
para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

clear tmp
for isubj = 1:length(SUBJLIST)
  for icond = 1 : 2
  tmp(:,:,isubj,icond) = tp_match_aal(para,fc_task(:,:,isubj,icond));
  end
  for igain = 1 : length(Gains)
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simtask(:,:,isubj,1,igain) = out.rc_wl; 
    pf(isubj,1,igain) = out.peakfreq; clear out
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,igain,v_sim))
    FC_simtask(:,:,isubj,2,igain) = out.rc_wl; 
    pf(isubj,2,igain) = out.peakfreq; clear out
  end
  
end

fc_task = tmp(include_bcn,include_bcn,:,:);

%%
igain = 8;

for isubj = 1:length(SUBJLIST)
  tmp_sim = FC_simtask(:,:,isubj,2,igain)-FC_simtask(:,:,isubj,1,igain);
  tmp_emp = fc_task(:,:,isubj,2)-fc_task(:,:,isubj,1);
  
  r_diff(isubj) = corr(tmp_sim(mask),tmp_emp(mask));
  
end

[h,p]=ttest(r_diff)
% Correlate the mean difference matrices
mean_fcd_sim = nanmean(FC_simtask(:,:,:,2,igain)-FC_simtask(:,:,:,1,igain),3);
mean_fcd_emp = nanmean(fc_task(:,:,:,2)-fc_task(:,:,:,1),3);

[r_mean,p_mean] = corr(mean_fcd_sim(mask),mean_fcd_emp(mask));

figure; set(gcf,'color','w'); 
subplot(2,2,1); 
imagesc(mean_fcd_sim,[-0.2 0.2]); axis square off
title('High gain vs baseline (Model)');
colorbar; colormap(redblue)

subplot(2,2,2); 
imagesc(mean_fcd_emp,[-0.05 0.05]); axis square off
title('Atx vs Pbo (MEG)')
colorbar; colormap(redblue)

%% CHANGE VS BASELINE

figure; set(gcf,'color','w');
change_sim = nanmean(FC_simtask(:,:,:,2,igain)-FC_simtask(:,:,:,1,1),3);
baseline_sim = nanmean(FC_simtask(:,:,:,2,igain)+FC_simtask(:,:,:,1,1),3);



subplot(2,2,1); hold on
plot(baseline_sim(mask),change_sim(mask),'.')
axis square; tp_editplots

change_emp = nanmean(fc_task(:,:,:,2)-fc_task(:,:,:,1),3);
baseline_emp = nanmean(fc_task(:,:,:,2)+fc_task(:,:,:,1),3);

subplot(2,2,2); hold on
plot(baseline_emp(mask),change_emp(mask),'.');
axis square; tp_editplots

subplot(2,2,3); hold on
plot(baseline_sim(mask),change_emp(mask),'.');
axis square; tp_editplots


%% PLOT ALTERED CORRELATIONS AS FUNCTION OF E-I RATIO
if exist(sprintf('~/pmod/proc/pmod_alteredcorr_asgain_and_ei_v%d.mat',v_sim))
  load(sprintf('~/pmod/proc/pmod_alteredcorr_asgain_and_ei_v%d.mat',v_sim))
  
else

clear p_pos_atx_sim p_neg_atx_sim
SUBJLIST = find(idx_rest.inh>15);
% gain = 18;
mask = logical(tril(ones(76,76),-1));

es = -18:1:18;
is = -18:1:18;


p_pos_atx_sim = nan(length(es),length(is),length(idx_gain1),2);
p_neg_atx_sim = nan(length(es),length(is),length(idx_gain1),2);


for ie = 1:length(es)
  ie
  for ii = 1:length(is)
    ii
    for igain = 1:length(idx_gain1)
%       igain
      if es(ie) == 0 && is(ii) == 0 && Gains(idx_gain1(igain)) == 0
        p_pos_atx_sim(ie,ii,igain,:) = [nan nan];
        p_neg_atx_sim(ie,ii,igain,:) = [nan nan];
        continue
      end
      
      gain = idx_gain1(igain);
      
      exc = es(ie);
      inh = is(ii);
      
      for isubj = SUBJLIST
%         isubj
        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),1,1,v_sim))
        FC_simrest(:,:,isubj,1) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(isubj)+exc,indiv_idx.rest.inh(isubj)+inh,1,gain,v_sim))
        FC_simrest(:,:,isubj,2) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),1,1,v_sim))
        FC_simtask(:,:,isubj,1) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(isubj)+exc,indiv_idx.task.inh(isubj)+inh,1,gain,v_sim))
        FC_simtask(:,:,isubj,2) = out.FC_env; clear out
      end
      
      [h,~,~,s]=ttest(atanh(FC_simrest(:,:,:,2)),atanh(FC_simrest(:,:,:,1)),'dim',3);
      p_pos_atx_sim(ie,ii,igain,1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
      p_neg_atx_sim(ie,ii,igain,1) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
      
      [h,~,~,s]=ttest(atanh(FC_simtask(:,:,:,2)),atanh(FC_simtask(:,:,:,1)),'dim',3);
      p_pos_atx_sim(ie,ii,igain,2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
      p_neg_atx_sim(ie,ii,igain,2) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
      
      %       [h,~,~,s]=ttest(FC_simtask(:,:,:,2),FC_simrest(:,:,:,2),'dim',3);
      %       p_pos_tvr_sim(ie,ii,igain) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
      %       p_neg_tvr_sim(ie,ii,igain) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
      %
      %       diff_tvr(ie,ii,igain) = sum(sum(triu(nanmean(FC_simtask(:,:,:,2),3)-nanmean(FC_simrest(:,:,:,2),3),1)))./sum(mask(:));
      %       fc_sim_rest(ie,ii,igain,:) = mean(squeeze(mean(mean(FC_simrest),2)));
      %       fc_sim_task(ie,ii,igain,:) = mean(squeeze(mean(mean(FC_simtask),2)));
      
    end
  end
end

% save "settings"
paar = [];
para.oscthresh = oscthres;
para.task_effect = task_effect;
para.prc_thresh = prc_thresh;
para.indiv_idx = indiv_idx;

save(sprintf('~/pmod/proc/pmod_alteredcorr_asgain_and_ei_v%d.mat',v_sim),'para','p_pos_atx_sim','p_neg_atx_sim');
end

%%
gg = Gains(idx_gain1);
gain = 6;
% PLOT
% -----------------------
figure; set(gcf,'color','w');

subplot(2,3,1); hold on
imagesc(p_pos_atx_sim(:,:,gain,1),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Rest - Pos (g=%.2f)',gg(gain)'))

subplot(2,3,4); hold on
imagesc(p_neg_atx_sim(:,:,gain,1),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Rest - Neg (g=%.2f)',gg(gain)'))

subplot(2,3,2); hold on
imagesc(p_pos_atx_sim(:,:,gain,2),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Task - Pos (g=%.2f)',gg(gain)'))

subplot(2,3,5); hold on
imagesc(p_neg_atx_sim(:,:,gain,2),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Task - Neg (g=%.2f)',gg(gain)'))

colormap(plasma)


subplot(2,3,3); hold on
imagesc(p_pos_atx_sim(:,:,gain,1)-p_pos_atx_sim(:,:,gain,2),[-0.2 0.2]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Ctx - Pos (g=%.2f)',gg(gain)')); colormap(gca,redblue)

subplot(2,3,6); hold on
imagesc(p_neg_atx_sim(:,:,gain,1)-p_neg_atx_sim(:,:,gain,2),[-0.2 0.2]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Ctx - Neg (g=%.2f)',gg(gain)'));colormap(gca,redblue)

for i = 1 : length(es)
  neg(i) = p_neg_atx_sim(i,i,gain,1)-p_neg_atx_sim(i,i,gain,2);
  pos(i) = p_pos_atx_sim(i,i,gain,1)-p_pos_atx_sim(i,i,gain,2);
end

figure; set(gcf,'color','w');
subplot(2,3,3); hold on
plot(pos,'r')
plot(neg,'b')
axis square; tp_editplots
axis([-1 38 -0.2 0.2]);
xlabel('\Delta(Background input)'); ylabel('Fraction of altered correlations [%]')
% axis([])



%%


%% LOAD PARAMETERS AND PLOT GAIN VS NO GAIN
% For specific level of gain
% ----------
clear p_* FC_simrest FC_simtask p_pos_atx_sim p_neg_atx_sim permdat_atx_rest permdat_atx_task perm
isperm = 1;
SUBJLIST
SUBJLIST = find(indiv_idx.rest.inh>15);

mask = logical(tril(ones(76,76),-1));
gain = 13;
%
if isperm
  load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim))
  
  for isubj = 1:length(SUBJLIST)
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simrest(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simrest(:,:,isubj,2) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simtask(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simtask(:,:,isubj,2) = out.FC_env; clear out
  end
  
  [h,~,~,s]=ttest(FC_simrest(:,:,:,2),FC_simrest(:,:,:,1),'dim',3);
  p_pos_atx_sim(1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:))
  p_neg_atx_sim(1) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:))
  
  [h,~,~,s]=ttest(FC_simtask(:,:,:,2),FC_simtask(:,:,:,1),'dim',3);
  p_pos_atx_sim(2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:))
  p_neg_atx_sim(2) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:))
  
  % --------------
  % PERMUTATION TEST: Test sign. altered correlations
  % ----------------------------
  
  nperm = 500;
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  
  for iperm = 1  : nperm
    fprintf('Permutation%d\n',iperm)
    % Atomoxetine
    idx1 = all_idx1(:,iperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat_atx_rest(:,:,i,1) = FC_simrest(:,:,i,idx1(i));
      permdat_atx_rest(:,:,i,2) = FC_simrest(:,:,i,idx2(i));
      permdat_atx_task(:,:,i,1) = FC_simtask(:,:,i,idx1(i));
      permdat_atx_task(:,:,i,2) = FC_simtask(:,:,i,idx2(i));
    end
    
    [h,~,~,s]=ttest(permdat_atx_rest(:,:,:,2),permdat_atx_rest(:,:,:,1),'dim',3);
    perm.p_pos_atx_sim(iperm,1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
    perm.p_neg_atx_sim(iperm,1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
    [h,~,~,s]=ttest(permdat_atx_task(:,:,:,2),permdat_atx_task(:,:,:,1),'dim',3);
    perm.p_pos_atx_sim(iperm,2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
    perm.p_neg_atx_sim(iperm,2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
    
  end
  
  for isubj = 1 : 28
    fc_sim_mean(isubj,1) = squeeze(outp.fc_sim_env_mean(indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),:,1));
    fc_sim_mean(isubj,2) = squeeze(outp.fc_sim_env_mean(indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),:,gain))-squeeze(outp.fc_sim_env_mean(indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),:,1));
  end
end

% --------------
% PLOT EVERYTHING
% --------------
figure; set(gcf,'color','w')
ax{1}= subplot(2,2,1); hold on

imagesc(zeros(size(osc,1),size(osc,2)))
par = outp.fc_sim_env_mean(:,:,:,gain)-outp.fc_sim_env_mean(:,:,:,1);
par(osc>oscthresh)=0;
imagesc(par,[-0.01 0.01])
set(ax{1},'YTick',1:40:length(Ies),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{1},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{1},'External input (to I)'); ylabel(ax{1},'External input (to E)');
axis(ax{1},'square')
axis([1 size(osc,1) 1 size(osc,2)])
tp_editplots; colormap(redblue); colorbar
% ax{2}= subplot(2,2,2); hold on
% m = mean(fc_sim_mean);
% s = std(fc_sim_mean);
% bar([1 2],p_pos_atx_sim,0.4,'facecolor','r','edgecolor','w');
% bar([1.5 2.5],p_neg_atx_sim,0.4,'facecolor','b','edgecolor','w');
% tp_editplots; axis([0.5 3 0 0.8])
% ylabel('Fraction of altered correlations')
% --------------

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_fits_indiv_gainmodulation_frac_gain%d_v%d.eps',gain,v_sim))


%%

% PLOT CONTEXT DEPENDENCE
subplot(4,2,5); hold on
plot(p_pos_atx_sim(:,1)-p_pos_atx_sim(:,2),'r');
plot(p_neg_atx_sim(:,1)-p_neg_atx_sim(:,2),'b');
plot(length(idx_gain)+2,emp.n_n_atx(6,1)-emp.n_n_atx(6,2),'o','markerfacecolor','b')
plot(length(idx_gain)+2,emp.n_p_atx(6,1)-emp.n_p_atx(6,2),'o','markerfacecolor','r')
axis([0 length(idx_gain)+2  -1 1]);
line([17 17],[-1 1],'color',[0.8 0.8 0.8])
set(gca,'XTick',[1:4:length(idx_gain) length(idx_gain)+2],'XTickLabels',[num2cell(Gains(idx_gain(1:4:length(idx_gain)))) 'Emp.'])
set(gca,'yTick',[-1 -0.5 0 0.5 1],'yTickLabels',num2cell([-1 0 1]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

subplot(4,2,2); hold on
plot(p_pos_tvr_sim,'r');
plot(p_neg_tvr_sim,'b');
axis([0 41 -0.05 0.65]);
line([17 17],[0 0.5],'color','k')
set(gca,'XTick',1:4:length(idx_gain),'XTickLabels',num2cell(Gains(idx_gain(1:4:length(idx_gain)))))
set(gca,'yTick',[-0.5 -0.25 0 0.25 0.50],'yTickLabels',num2cell([-0.5 -0.25 0 0.25 0.50]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

subplot(4,2,4); hold on
plot(diff_tvr,'k');
axis([0 31 -0.01 0.1]);
line([17 17],[0 0.5],'color','k')
set(gca,'XTick',1:4:length(idx_gain),'XTickLabels',num2cell(Gains(idx_gain(1:4:length(idx_gain)))))
set(gca,'yTick',[0 0.05 0.10],'yTickLabels',num2cell([0 0.25 0.50]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_percchange_as_gain_v%d.pdf',v_sim))

%% PLOT PEAK FREQ AND ALPHA FOR SUBJECT PARAMETERS

for isubj = 1 : 28
  if ~isnan(idx_rest.exc(isubj))
    pf(isubj,1) = squeeze(outp.peakfreq(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1));
    pf(isubj,2) = squeeze(outp.peakfreq(idx_task.exc(isubj),idx_task.inh(isubj),1,1));
    
    alp(isubj,1) = squeeze(outp.alphapow(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1));
    alp(isubj,2) = squeeze(outp.alphapow(idx_task.exc(isubj),idx_task.inh(isubj),1,1));
    
  else
    pf(isubj,:) = nan;
    alp(isubj,:) = nan;
  end
end

figure; set(gcf,'color','w'); hold on
subplot(3,2,1); hold on
m1 = nanmean(pf);
s1 = nanstd(pf)/sqrt(sum(sum(~isnan(pf)))/2);

bar([1,2],m1); hold on
line([1 1],[m1(1)-s1(1) m1(1)+s1(1)])
line([2 2],[m1(2)-s1(2) m1(2)+s1(2)])

axis([0.5 2.5 10 15]); tp_editplots; ylabel('Peak frequency [Hz]')

subplot(3,2,2); hold on
m1 = nanmean(alp);
s1 = nanstd(alp)/sqrt(sum(sum(~isnan(alp)))/2);

bar([1,2],m1); hold on
line([1 1],[m1(1)-s1(1) m1(1)+s1(1)])
line([2 2],[m1(2)-s1(2) m1(2)+s1(2)])

axis([0.5 2.5 0.1e-04 1.1e-04]); tp_editplots; ylabel('Peak frequency [Hz]')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_tvr_alpha_and_peakfreq_v%d.pdf',v_sim))

%%
figure; set(gcf,'color','white');
subplot(2,2,1)
par = outp.fc_sim_env_mean(:,:,:,11)-outp.fc_sim_env_mean(:,:,:,1);
par(osc>oscthresh)=0;
imagescnan(par,[-0.02 0.02]);
set(gca,'ydir','normal'); colormap(redblue);
axis([1 size(outp.fc_sim_env_mean,1) 1 size(outp.fc_sim_env_mean,1)])
%  set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel('Background input (to I)'); ylabel('Background input (to E)');
axis tight square;tp_editplots;
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_wc_final_fitting_gain_phase_v%d.eps',v));


%% INDIVIDUAL *GAIN* FITS BASED ON CORRELATION
mask = logical(tril(ones(76,76),-1));

v = 1;
para.nfreq = 1:13;
para.alpha = 0.05;

% if ~exist('emp','var')
if ~exist('cleandat','var')
  cleandat = pupmod_loadpowcorr(v,1);
end
% transform indiv. subj. matrices to AAL BCN
for isubj = 1 : 28
  isubj
  for im = 1 : 3
    for icont = 1 : 2
      for ifoi = 6:7
        tmp = squeeze(cleandat(:,:,isubj,im,icont,ifoi));
        tmp = tp_match_aal(para,tmp);
        fc_indiv(:,:,isubj,im,icont,ifoi) = tmp(include_bcn,include_bcn);
      end
    end
  end
end
% end
%%
clear cleandat
[~,idx_gain]=sort(Gains)

SUBJLIST = find(idx_rest.inh>15);

for icond = 1 : 2
  icond
  for isubj = SUBJLIST
    for igain = 1 : length(Gains)
      
      gain = idx_gain(igain)
      
      load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),1,gain,v_sim))
      tmp = fc_indiv(:,:,isubj,icond,1,6);
      if isubj == 1; figure;imagesc(out.FC_env); end;
      
      r(1,icond,igain,isubj) = 1-(corr(tmp(mask),out.FC_env(mask))-(mean(tmp(mask))-mean(out.FC_env(mask)))^2); clear out
      
      load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),1,gain,v_sim))
      tmp = fc_indiv(:,:,isubj,icond,2,6);
      r(2,icond,igain,isubj) = 1-(corr(tmp(mask),out.FC_env(mask))-(mean(tmp(mask))-mean(out.FC_env(mask)))^2); clear out
      
    end
  end
end

%% PLOT
figure; set(gcf,'color','w');

subplot(2,2,1); hold on
imagesc(squeeze(r(1,1,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))
colorbar; xlabel('Gains'); ylabel('Subjects'); xtickangle(90); title('Pbo & Rest')

subplot(2,2,3); hold on
imagesc(squeeze(r(1,2,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))
colorbar; xlabel('Gains'); ylabel('Subjects'); xtickangle(90); title('Atx & Rest')

subplot(2,2,2); hold on
imagesc(squeeze(r(2,1,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))
colorbar; xlabel('Gains'); ylabel('Subjects'); xtickangle(90); title('Pbo & Task')

subplot(2,2,4); hold on
imagesc(squeeze(r(2,2,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))
colorbar; xlabel('Gains'); ylabel('Subjects'); xtickangle(90); title('Atx & Task')

for i = 1 : 42
  for isubj = 1 : 28
    kk = i:i+2;
    smoothed_task(i,isubj,1) = mean((r(2,1,kk,isubj)));
    smoothed_task(i,isubj,2) = mean((r(2,2,kk,isubj)));
  end
end
%%
[i,j]=min(squeeze(smoothed_task(:,:,2)));
% [i,jj]=min(r(2,1,:,:)); jj = squeeze(jj);

SUBJLIST1  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
para.str_behav = 'count';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav1 = nanmean(behav,3)';


para.str_behav = 'numb_switches';
behav = pconn_read_behavioral_data(SUBJLIST1,para);
behav = nanmean(behav,3)';

behav = (behav1 + behav)./2;

ss = find(~isnan(behav(:,2)));

ss = intersect(ss,SUBJLIST);
scatter(behav(ss,2)-behav(ss,1),j(ss))


%%
isubj = 1;
figure;
subplot(1,2,1);
imagesc(abs(mean(fc_indiv(:,:,:,3,1,7),3)-mean(fc_indiv(:,:,:,1,1,7),3)),[0 0.0075]); axis square off
subplot(1,2,2)
imagesc(abs(mean(fc_indiv(:,:,:,2,2,6),3)-mean(fc_indiv(:,:,:,1,2,6),3)),[0 0.015]); axis square off

atx = mean(fc_indiv(:,:,:,2,2,7),3)-mean(fc_indiv(:,:,:,1,2,7),3);
atx = atx(mask);

atx_bl = mean(fc_indiv(:,:,:,1,2,7),3);
atx_bl = atx_bl(mask);


r_emp_atx = corr(atx_bl,atx);

dpz = mean(fc_indiv(:,:,:,3,1,7),3)-mean(fc_indiv(:,:,:,1,1,7),3);
dpz = dpz(mask);

dpz_bl = mean(fc_indiv(:,:,:,1,1,7),3);
dpz_bl = dpz_bl(mask);

r_emp_dpz = corr(dpz_bl,dpz);

nperm = 3000;

all_idx1 = randi(2,[28,nperm]);

for iiperm = 1 : nperm
  iiperm
  idx1 = all_idx1(:,iiperm);
  idx2 = 3-idx1;
  
  for i = 1 : length(idx1)
    permdat_cnt1(:,:,i,1) = fc_indiv(:,:,i,idx1(i),2,7);
    permdat_cnt1(:,:,i,2) = fc_indiv(:,:,i,idx2(i),2,7);
  end
  
  idx1 = all_idx1(:,iiperm);
  idx2 = 3-idx1;
  idx1(idx1 == 2) = 3;
  idx2(idx2 == 2) = 3;
  
  for i = 1 : length(idx1)
    permdat_cnt2(:,:,i,1) = fc_indiv(:,:,i,idx1(i),1,7);
    permdat_cnt2(:,:,i,2) = fc_indiv(:,:,i,idx2(i),1,7);
  end
  
  
  atx_perm = mean(permdat_cnt1(:,:,:,2),3)-mean(permdat_cnt1(:,:,:,1),3);
  atx_perm = atx_perm(mask);
  
  atx_bl_perm = mean(permdat_cnt1(:,:,:,1),3);
  atx_bl_perm = atx_bl_perm(mask);
  
  r_perm_atx(iiperm) = corr(atx_bl_perm,atx_perm);
  
  dpz_perm = mean(permdat_cnt2(:,:,:,2),3)-mean(permdat_cnt2(:,:,:,1),3);
  dpz_perm = dpz_perm(mask);
  
  dpz_bl_perm = mean(permdat_cnt2(:,:,:,1),3);
  dpz_bl_perm = dpz_bl_perm(mask);
  
  r_perm_dpz(iiperm) = corr(dpz_bl_perm,dpz_perm);
  
end

%% LOAD CLEAD PUPIL FROM ALENA

ord   = pconn_randomization;
for isubj = 4 : 34
  for iblock = [1 2]
    for m = 1 :3
      im = find(ord(isubj,:)==m);
      try
        load(sprintf('/home/arussmann/pupil_data/pconn_pup_samples_s%d_b%d_m%d_no_blinks.mat',isubj,iblock,im))
        isubj
        pp(isubj,m,iblock) = mean(pupil);
      catch me
        pp(isubj,m,iblock) = nan'
      end
    end
  end
end

pp = pp(SUBJLIST,:,:);

for ia = 1 : 28
  
  ppp(ia,:,:) = reshape((reshape(pp(ia,:,:),[6 1])-nanmean(reshape(pp(ia,:,:),[6 1])))./nanstd(reshape(pp(ia,:,:),[6 1])),[3 2])
end

%% PLOT JUST FC IN MODEL FOR VARIOUS GLOBAL COUPLING AND 2 GAINS
% v1
% osc=osc

% oscthres = 0.5;
for ig = 31
%   osc = osc1(:,:,ig,1);

  figure; set(gcf,'color','w');
  
  subplot(2,3,1); 
  
  imagesc(outp.fc_sim_env_mean(:,:,ig,1),[0 0.5]);
  axis square; colormap(gca, 'plasma')
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal');
  xlabel('Input to E'); ylabel('Input to I')
  title('FC: Gain = 1')
 
  subplot(2,3,2);
  
  imagesc(outp.fc_sim_env_mean(:,:,ig,2),[0 0.5]);
  axis square; colormap(gca, 'plasma')
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('FC: Gain = 2')
  
  subplot(2,3,3);
  
  par = outp.fc_sim_env_mean(:,:,ig,2)-(outp.fc_sim_env_mean(:,:,ig,1));
  par(osc>oscthres)=nan;
  b=imagesc(par,[-0.02 0.02]); set(b,'AlphaData',~isnan(par))
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('\DeltaFC (2-1)')
  
  subplot(2,3,4);
  par = outp.r_env_rest_corr(:,:,ig,1);
  par(osc>oscthres)=nan;
  b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isnan(par)) 
  axis square; colormap(gca, plasma)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('Corr w FC_{emp} (Gain: 1)')
  
  subplot(2,3,5);
  par = outp.r_env_rest_corr(:,:,ig,2);
  par(osc>oscthres)=nan;
  b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isnan(par)) 
  axis square; colormap(gca, plasma)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('Corr w FC_{emp} (Gain: 2)')
  
  subplot(2,3,6);
  
  imagesc(r_task(:,:,ig),[-0.5 0.5]);
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('Corr w diff')
  
  
  print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_effect_of_alpha_and_gain_on_FC_G%d_v%d.pdf',ig,v))
  
  figure; set(gcf,'color','w');
  subplot(2,3,2);
  par = r_rest(:,:,ig);
  par(osc1(:,:,ig)>oscthres)=nan
  b=imagesc(par,[-0.5 0.5]); set(b,'AlphaData',~isnan(par))
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('Corr w diff (Rest/Corr)')
  
  subplot(2,3,3);
  par = r_task(:,:,ig);
  par(osc1(:,:,ig)>0.5)=nan
  b=imagesc(par,[-0.5 0.5]); set(b,'AlphaData',~isnan(par))
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E')
  title('Corr w diff (Task/Corr)')
  
  k=subplot(2,3,5);
  par = dist_rest(:,:,ig);
  par(osc1(:,:,ig)>oscthres)=nan
  b=imagesc(par,[0.7 1.3]); set(b,'AlphaData',~isnan(par))
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E'); colormap(k,plasma)
  title('Corr w diff (Rest/Dist)')
  
  k=subplot(2,3,6);
  par = dist_task(:,:,ig);
  par(osc1(:,:,ig)>oscthres)=nan
  b=imagesc(par,[0.7 1.3]); set(b,'AlphaData',~isnan(par))
  axis square; colormap(gca, cmap)
  set(gca,'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  tp_editplots; set(gca,'ydir','normal')
  xlabel('Input to E'); colormap(k,plasma)
  title('Corr w diff (Task/Dist)')
  
    print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_effect_of_alpha_and_gain_on_correlationswithFC_G%d_v%d.pdf',ig,v))

end

%%
fcd_task = squeeze(nanmean(nanmean(cleandat(:,:,:,2,2,:),3),6)-nanmean(nanmean(cleandat(:,:,:,1,2,:),3),6));
fcd_rest = squeeze(nanmean(nanmean(cleandat(:,:,:,2,1,:),3),6)-nanmean(nanmean(cleandat(:,:,:,1,1,:),3),6));
v = 2;
para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

fcd_task = tp_match_aal(para,fcd_task);
fcd_task = fcd_task(include_bcn,include_bcn);
fcd_rest = tp_match_aal(para,fcd_rest);
fcd_rest = fcd_rest(include_bcn,include_bcn);

for ig = 1: length(Gg)
  for igain = 1 : length(Gains)
    for i = 1 : length(Ies)
      i
      for j = 1 :length(Iis)

        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v,i,j,ig,1,v))
        fc_sim1 = out.rc_wl;

        load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v,i,j,ig,igain,v))
        fc_sim2 = out.rc_wl;

        fcd_sim = fc_sim2-fc_sim1;

        [r_task(i,j,ig,igain),p_task(i,j,ig,igain)]=corr(fcd_sim(mask),fcd_task(mask));
        [r_rest(i,j,ig,igain),p_rest(i,j,ig,igain)]=corr(fcd_sim(mask),fcd_rest(mask));

        dist_rest(i,j,ig,igain) = 1-(r_rest(i,j,ig,igain)-(mean(fcd_sim(mask))-mean(fcd_rest(mask))).^2);
        dist_task(i,j,ig,igain) = 1-(r_task(i,j,ig,igain)-(mean(fcd_sim(mask))-mean(fcd_task(mask))).^2);
      end
    end
  end
end
    
save(sprintf('~/pmod/proc/pmod_correlationmodeldifferences_v%d.mat',v),'r_task','r_rest','dist_rest','dist_task','p_task','p_rest')

%% OLD TASK FITTING CODE
% ------------------------------------
% clear idx_task idx
% 
% indiv_change_prc = 100*(mean(fc_task_indiv,1)-mean(fc_rest_indiv,1))./mean(fc_rest_indiv,1);
% 
% for isubj = 1:28
%   
%   if isnan(idx_rest.exc(isubj))
%     idx_task.exc(isubj) = nan;
%     idx_task.inh(isubj) = nan;
%     continue
%   end
%   
%   if dist
%     
%     par                 = squeeze(outp.dist_task_indiv(isubj,1:lim,1:lim,iG,igain));
%     par(osc>oscthresh)  = nan;
%     dist_indiv(isubj)   = mean(par(par<prctile(par(:),prc_thresh)));
%     par                 = par<prctile(par(:),prc_thresh);
%     bw                  = bwlabel(double(par),8);
%     cnt                 = [];
%     
%     for iclust = 1 : max(bw(:))
%       cnt(iclust) = sum(bw(:)==iclust);
%     end
%     
%     % if no cluster(s), omit subject from further analysis
%     if isempty(cnt)
%       idx_rest.inh(isubj) = nan;
%       idx_rest.exc(isubj) = nan;
%       continue
%     end
%     
%     % find single largest cluster
%     [~,k(isubj)]       = max(cnt);
%     par(bw==k(isubj))  = 1;
%     par(bw~=k(isubj))  = 0;
%     
%     [i,kk]=find(par==1);
%     idx(isubj,:) = [round(mean(i)) round(mean(kk))];
%     
%     idx_task.inh(isubj) = idx(isubj,2);
%     idx_task.exc(isubj) = idx(isubj,1);
%     %     [i,j]            
%   else
%   
% %   p_corr  = fdr1(reshape(squeeze(outp.p_env_task_indiv_corr(isubj,1:lim,1:lim,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.05);
% %   par     = squeeze(outp.p_env_task_indiv_corr(isubj,1:lim,1:lim,iG,igain)) < 0.05;
% %   par(osc>oscthresh) = 0;
%   
% %   fc_tvr_sim = 100*(outp.fc_sim_env_mean(1:lim,1:lim,1,1)-outp.fc_sim_env_mean(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1))./outp.fc_sim_env_mean(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1);
% %   d = indiv_change_prc(isubj)-fc_tvr_sim;
% %   d(osc>oscthresh) = inf;
% %   d(par==0) = inf;
%   
%   % DISTANCE METRIC FROM DEMIRTAS??
% %   mask = zeros(lim,lim);
% %   mask(idx_rest.exc(isubj):end,idx_rest.inh(isubj):end) = 1;
% %   d(~logical(mask)) = Inf;
% %   
% %   %   d = abs(d.*(d<0));
% %   %   d(d==0)=Inf;
% %   m = min(reshape(abs(d),[lim*lim 1]));
% %   [i,j]=find(abs(d)==m);
% %   idx = round([mean(i) mean(j)]);
% %   idx_task.inh(isubj) = idx(2);
% %   idx_task.exc(isubj) = idx(1);
%   end
%   
% end
% 
% indiv_idx.rest = idx_rest;
% indiv_idx.task = idx_task;
