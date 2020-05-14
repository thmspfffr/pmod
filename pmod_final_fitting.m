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
% Gains       = [0]; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
% dt          = 0.01;
%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 1.4;
% Gains       = [0:0.05:0.6 -0.2:0.05:-0.05 0.65:0.05:1]; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
% dt          = 0.01;
%-------------------------------------------------------------------------
% VERSION 3: 27/01/2020
% %-------------------------------------------------------------------------
v           = 3;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = [1.2:-0.01:1.15];
Gains       = [-0.1:0.02:0.12]; 
nTrials     = 1;
tmax        = 6500;  % in units of tauE
dt          = 0.01;
%-------------------------------------------------------------------------

% load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_parameters_v%d.mat',v,v))

% Gains = [0 0.025:0.025:0.4 -0.025:-0.025:-0.3 ]
v_sim = v;
% connectivity, AAL
v_conn =  33;
addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

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
fc_rest = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,1,[3 4 5 6]),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,2,[3 4 5 6]),3),6));

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

fc_rest = tp_match_aal(para,fc_rest);
fc_task = tp_match_aal(para,fc_task);

fc_rest = fc_rest(include_bcn,include_bcn);
fc_task = fc_task(include_bcn,include_bcn);

% transform indiv. subj. matrices to AAL BCN
% for ifreq = 5
for isubj = 1 : 28
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,1,:),6));
  tmp = tp_match_aal(para,tmp); 
  fc(:,:,isubj)=tmp;
  tmp = tmp(include_bcn,include_bcn);
%   [corrwithfc_rest(isubj,ifreq), p_corrwithfc_rest(isubj)]  = corr(tmp(mask),SC(mask));
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,2,:),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  fc_task_indiv(:,isubj) = tmp(mask);
%   [corrwithfc_task(isubj,ifreq), p_corrwithfc_task(isubj)] = corr(tmp(mask),SC(mask));
end
% end
k = [0 0];

%%
% v_sim=3
% if ~exist(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  
  for iies = 1 : length(Ies)
    iies
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        for iiis = 1 :length(Iis)
          
          load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,iies,iiis,iG,igain,v_sim))
          
          outp.fc_sim_fr_tmp = out.fc_FR;
          
          [outp.r_fr_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest(mask));
%           [outp.r_fr_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_task(mask));          
          
          [outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest_indiv);
%           [outp.r_fr_task_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_task_indiv);
          
          outp.dist_fr_rest(iies,iiis,iG,igain) = 1-(outp.r_fr_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_fr_tmp(mask))).^2);
%           outp.dist_fr_task(iies,iiis,iG,igain) = 1-(outp.r_fr_task_corr(iies,iiis,iG,igain)-(mean(fc_task(mask))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          
          outp.dist_fr_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_fr_tmp(mask))).^2);
%           outp.dist_fr_task_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_fr_task_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_task_indiv))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          
          outp.fc_sim_fr_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_fr_tmp(mask));
          
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          if isempty( out.peakfreq)
             out.peakfreq = nan;
          end
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          outp.alphapow(iies,iiis,iG,igain) = mean(out.alphapow);   
          
          % DO THE SAME AS ABOVE WITH ENVELOPES
          % ------------------------
          
%           outp.fc_sim_env_tmp = out.rc_wl;
          
%           [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
%           [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));          
%           
%           [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv);
%           [outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_task_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task_indiv);
%           
%           outp.dist_env_rest(iies,iiis,iG,igain) = 1-(outp.r_env_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
%           outp.dist_env_task(iies,iiis,iG,igain) = 1-(outp.r_env_task_corr(iies,iiis,iG,igain)-(mean(fc_task(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
%           
%           outp.dist_env_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
%           outp.dist_env_task_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_task_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
%           
%           outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));      
%           
        end
      end
    end
  end
%   
  save(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim),'outp')
% % 
% else
%   load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  load(sprintf('~/pmod/proc/detosc/v%d/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim,v_sim))
% end
% osc1=osc1(1:121,1:121,:,:);
% clear cleandat

error('!')

%% PEAK FREQUENCY

igain = find(Gains==0);

figure_w
ax{1} = subplot(2,3,1);
par=outp.peakfreq(:,:,:,igain);
par(osc1(:,:,:,igain)>0)=nan;
b=imagesc(par,[5 20]); colormap(plasma);
set(b,'AlphaData',~isnan(par)); set(gca,'ydir','normal')

ax{2} = subplot(2,3,2);
b=imagesc(par>7.8814 & par<18.7452)
set(b,'AlphaData',~isnan(par)); set(gca,'ydir','normal')


for iax = [1 2 ]
%   scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'External input (to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax ==4
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E')
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  end
  
  tp_editplots(ax{iax})
  
%   c = colorbar(ax{iax}); axis(ax{iax},'tight')
%   c.Ticks = [min(c.Limits) max(c.Limits)];
%   if iax == 3
%     c.TickLabels = {'0';'1'};
%   end
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
%   c.Ticks = [min(c.Ticks) max(c.Ticks)];
  axis(ax{iax},'square')
  
end

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_peakfreq_igain%d_v%d.eps',igain,v_sim))

%%
lim = 121;

oscthres = 0;
nancol = [0.97 0.97 0.97];

clear par ax
figure; set(gcf,'color','w')

igain = 6;
iG = 1;
osc = osc1(1:lim,1:lim,iG,igain);

% plot peak freq model
ax{1} = subplot(3,2,1); hold on
par = outp.fc_sim_fr_mean(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.02]); set(b,'AlphaData',~isnan(par))

% plot corr fc_model with fc_meg
ax{3} = subplot(3,2,3); hold on
par = outp.r_fr_rest_corr(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isnan(par))
% title('r(FC_{im},FC_{MEG})');

% plot peak freq model
ax{2} = subplot(3,2,2); hold on
par = outp.fc_sim_env_mean(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.02]); set(b,'AlphaData',~isnan(par))

% plot corr fc_model with fc_meg
ax{4} = subplot(3,2,4); hold on
par = outp.r_env_rest_corr(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.4]); set(b,'AlphaData',~isnan(par))
% title('r(FC_{sim},FC_{MEG})');

% 
% k=(par>prctile(par(:),99));
% clust = bwlabel(k,8);
% for iclust = 1 : max(clust(:))
%   n_clust(iclust)=sum(clust(:)==iclust);
% end
% 
% [~,k] = max(n_clust);
% par = clust == k
% m2 = par>0;
% [i,j]= find(m2);
% idx = round([mean(i) mean(j)]);
% 
% scatter(idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
% 
% par = outp.dist_rest(:,:,iG,igain); par(osc>oscthres)=nan;
% prc = prctile(reshape(par,[lim*lim 1 1]),1);
% 
% ax{4} = subplot(2,3,3); hold on
% par = par<prc;
% par(osc>oscthres)=0;
% clust = bwlabel(par,8);
% for iclust = 1 : max(clust(:))
%   n_clust(iclust)=sum(clust(:)==iclust);
% end

% [~,k] = max(n_clust);
% par = clust == k
% m2 = par>0;
% par = double(par); par(par<1)=NaN;
% b=imagesc(par,[0 1]); set(b,'AlphaData',~isnan(par))
% title('Correlation: p-values (log10)');

% [i,j]= find(m2);
% idx_avg = round([mean(i) mean(j)]);

for iax = [1 2 3 4]
%   scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 2
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  elseif iax ==4
    colormap(ax{iax},plasma)
%     colormap(ax{iax},[nancol; 1 0 0])
    xlabel(ax{iax},'Background input to I')
    ylabel(ax{iax},'Background input to E');
    set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  end
  
  tp_editplots(ax{iax})
  
%   c = colorbar(ax{iax}); axis(ax{iax},'tight')
%   c.Ticks = [min(c.Limits) max(c.Limits)];
%   if iax == 3
%     c.TickLabels = {'0';'1'};
%   end
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
%   c.Ticks = [min(c.Ticks) max(c.Ticks)];
  axis(ax{iax},'square')
  
end

set(gcf,'Renderer','Painters')
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_fc_frandenv_iG%d_v%d.eps',iG,v_sim))

%% FIT RESTING STATE PLACEBO 
% (THROUGH CORRELATION OR DISTANCE)
% ----------------------------------
load ~/pmod/proc/numerical/ksdistance.mat

% osc=permute(repmat(osc,[1 1 1 16 28]),[5,1,2,3,4]);
% outp.r_fr_rest_indiv_corr(osc>0)=nan;
clear idx par ax cnt idx_rest par
lim = 121;
iG         = 1; 
igain      = find(Gains==0);
osc        = osc1(1:lim,1:lim,iG,igain);
oscthresh  = 0;
prc_thresh = 2.5;

h=figure; set(gcf,'color','w')

% Loop through all subjects
for isubj = 1 : 28
  clear par bw
  
    clear k
    par                 = squeeze(outp.dist_fr_rest_indiv(isubj,1:lim,1:lim,iG,igain));    
    par(nanmean(osc1(:,:,iG,6:12),4)>0)  = nan;
    binmap              = par<prctile(par(:),prc_thresh);
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
    
    r = squeeze(outp.dist_fr_rest_indiv(isubj,1:lim,1:lim,iG,igain));
    [i,kk]=find(par==1);
    idx(isubj,:) = round([mean(i) mean(kk)]);
    
    % PLOT INDIVIDUAL SUBJECT FITTING RESULT
    % ------------------------------------------------
    if isubj == 1
      k=figure_w; subplot(3,2,1)
      par1 = par;
      par1(par1==0)=nan;
      imagesc(par1)
      b=imagesc(par1,[0 0.4]); set(b,'AlphaData',~isnan(par1))
      colormap(jet)

      xlabel(gca,'Background input to I')
      ylabel(gca,'Background input to E');
      set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
      set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)),'ydir','normal')
      tp_editplots

      axis(gca,[1 length(Iis) 1 length(Ies) ])
      axis(gca,'square')
      
      hold on
      scatter(idx(isubj,2),idx(isubj,1),5,'o')
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_indivfitthresh_iG%d_v%d.pdf',igain,v))
close
    end
    % ------------------------------------------------
    
   % PLOT RESULT
  ax{1}   = subplot(6,5,isubj); hold on
  par = squeeze(outp.dist_fr_rest_indiv(isubj,1:lim,1:lim,iG,igain));
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



%% PLOT HIGH GAIN VS BASELINE GAIN

[Mu,ia,ic]=unique([idx_rest.inh' idx_rest.exc'],'rows','stable')
c=accumarray(ic,1);
cols = cbrewer('seq', 'YlOrRd', max(c)+3);
cols = cols(4:end,:);

SUBJ = 1:28; %SUBJ(17)=[];
clear fc_mod_rest fc_mod_task

igain = find(Gains==0);

inh_avg(1) = round(mean(idx_rest.inh(SUBJ)));
exc_avg(1) = round(mean(idx_rest.exc(SUBJ)));

inh_avg(2) = round(mean(idx_rest.inh(SUBJ)))+15;
exc_avg(2) = round(mean(idx_rest.exc(SUBJ)))+10;

h = figure_w; 

subplot(2,3,1); hold on
par = outp.fc_sim_fr_mean(:,:,:,igain);
par(osc1(:,:,iG,igain)==1)=nan;
% 
% b=imagesc(par,[0 0.02]); colormap(plasma)
% set(b,'AlphaData',~isnan(par))

for isubj = 1:size(Mu,1)
  
  scatter(Mu(isubj,1),Mu(isubj,2),10,'o','markerfacecolor',cols(c(isubj),:),'markeredgecolor','none');
  
end

line([40 40],[40 50],'color','k')
line([40 60],[50 50],'color','k')
line([60 60],[40 50],'color','k')
line([40 60],[40 40],'color','k')

axis(gca,[1 121 1 121 ])
tp_editplots; axis square
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel('Background input to I'); ylabel('Background input to E'); 

subplot(2,3,4); hold on
for isubj = 1:size(Mu,1)
  scatter(Mu(isubj,1),Mu(isubj,2),20,'o','markerfacecolor',cols(c(isubj),:),'markeredgecolor','none');
end

axis(gca,[40 61 40 51 ])
tp_editplots; %axis square
set(gca,'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
set(gca,'XTick',1:5:length(Iis),'XTickLabels',num2cell(Iis(1:5:end)))
xlabel('Background input to I'); ylabel('Background input to E'); 

baseline_gain=find(Gains==0);

load redblue.mat

for igain = 8
  
  for isubj = SUBJ
      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,baseline_gain,v_sim))
      fc_mod_rest(isubj,1) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,1,1) = out.fc_FR(mask);
      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,igain,v_sim))
      fc_mod_rest(isubj,2) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,2,1) = out.fc_FR(mask);
      
      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+10,idx_rest.inh(isubj)+15,iG,baseline_gain,v_sim))
      fc_mod_task(isubj,1) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,1,2) = out.fc_FR(mask);
      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+10,idx_rest.inh(isubj)+15,iG,igain,v_sim))
      fc_mod_task(isubj,2) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,2,2) = out.fc_FR(mask);
      
  end
  
  for isubj = SUBJ

      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)-14,idx_rest.inh(isubj)-19,iG,igain,v_sim))
      fc_mod_rest(isubj,3) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,3,1) = out.fc_FR(mask);
      
      load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)-4,idx_rest.inh(isubj)-4,iG,igain,v_sim))
      fc_mod_task(isubj,3) = mean(out.fc_FR(mask));
      fc_all_mod(:,isubj,3,2) = out.fc_FR(mask);
      
  end

  par = outp.fc_sim_fr_mean(:,:,:,igain)-outp.fc_sim_fr_mean(:,:,:,baseline_gain);
  par(osc1(:,:,1,baseline_gain)>0)=nan;
  par(osc1(:,:,1,igain)==1)=nan;

  subplot(2,3,2); hold on
  b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');  
%   for isubj = 1 : 28
% 
%     scatter(idx_rest.inh(isubj),idx_rest.exc(isubj),8,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','w');
% 
%   end

  scatter(round(mean(idx_rest.inh(SUBJ))),round(mean(idx_rest.exc(SUBJ))),20,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');
  scatter(round(mean(idx_rest.inh(SUBJ)))+15,round(mean(idx_rest.exc(SUBJ)))+10,20,'o','markerfacecolor','y','markeredgecolor','k');

  set(b,'AlphaData',~isnan(par));
  colormap(gca,redblue)
  ylabel(gca,'External input (to E)');
  set(gca,'YTick',1:40:121,'YTickLabels',num2cell(Ies(1:40:end)))
  set(gca,'XTick',1:40:121,'XTickLabels',num2cell(Iis(1:40:end)))

  tp_editplots; axis square; axis([1 121 1 121])
  
  par = outp.fc_sim_fr_mean(:,:,:,igain)-outp.fc_sim_fr_mean(:,:,:,baseline_gain);
  par(osc1(:,:,1,baseline_gain)>0)=nan;
  par(osc1(:,:,1,igain)==1)=nan;
  
  subplot(2,3,5); hold on
  b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');  

  scatter(round(mean(idx_rest.inh))-4,round(mean(idx_rest.exc))-4,20,'o','markerfacecolor','y','markeredgecolor','k');
  scatter(round(mean(idx_rest.inh))-19,round(mean(idx_rest.exc))-14,20,'o','markerfacecolor','w','markeredgecolor','k');
%   scatter(round(mean(idx_rest.inh))-12,round(mean(idx_rest.exc))-10,20,'o','markerfacecolor','y','markeredgecolor','k');

  set(b,'AlphaData',~isnan(par));
  colormap(gca,redblue)
  ylabel(gca,'External input (to E)');
  set(gca,'YTick',1:40:121,'YTickLabels',num2cell(Ies(1:40:end)))
  set(gca,'XTick',1:40:121,'XTickLabels',num2cell(Iis(1:40:end)))

  tp_editplots; axis square; axis([1 121 1 121])
  
end

  subplot(2,3,3); hold on
  
  bar(1,mean(fc_mod_rest(:,1)),'facecolor',[.8 .8 .8])
  bar(2,mean(fc_mod_rest(:,2)),'facecolor',[1 1 1])
  bar(3,mean(fc_mod_task(:,1)),'facecolor',[1 0 0])
  bar(4,mean(fc_mod_task(:,2)),'facecolor',[1 .5 .5])
  
  axis square; tp_editplots
  axis([0 5 0 0.03]);
  set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Rest (Pbo)';'Rest (Atx)';'Task (Pbo)';'Task (Atx)'})
  xtickangle(gca,45)
  
  subplot(2,3,6); hold on
  
  bar(1,mean(fc_mod_rest(:,1)),'facecolor',[.8 .8 .8])
  bar(2,mean(fc_mod_rest(:,3)),'facecolor',[1 1 1])
  bar(3,mean(fc_mod_task(:,1)),'facecolor',[1 0 0])
  bar(4,mean(fc_mod_task(:,3)),'facecolor',[1 .5 .5])
  
  axis square; tp_editplots
  axis([0 5 0 0.03]);
  set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Rest (Pbo)';'Rest (Dpz)';'Task (Pbo)';'Task (Dpz)'})
  xtickangle(gca,45)
  
% 
% for igain = 1 : 16
% 
%   diff_rest(igain) = outp.fc_sim_fr_mean(exc_avg(1),inh_avg(1),:,igain)-outp.fc_sim_fr_mean(exc_avg(1),inh_avg(1),:,baseline_gain);
%   diff_task(igain) = outp.fc_sim_fr_mean(exc_avg(2),inh_avg(2),:,igain)-outp.fc_sim_fr_mean(exc_avg(2),inh_avg(2),:,baseline_gain);
% 
% end
% 

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_highvslowgain_gain%d_v%d.pdf',igain,v_sim))

% 
[h,~,~,s]=ttest(fc_all_mod(:,:,2,1),fc_all_mod(:,:,1,1),'dim',2);
atx_p(1) = sum(h>0&s.tstat>0)/2850;
atx_n(1) = sum(h>0&s.tstat<0)/2850;
[h,~,~,s]=ttest(fc_all_mod(:,:,2,2),fc_all_mod(:,:,1,2),'dim',2);
atx_p(2) = sum(h>0&s.tstat>0)/2850;
atx_n(2) = sum(h>0&s.tstat<0)/2850;
[h,~,~,s]=ttest(fc_all_mod(:,:,3,1),fc_all_mod(:,:,1,1),'dim',2);
dpz_p(1) = sum(h>0&s.tstat>0)/2850;
dpz_n(1) = sum(h>0&s.tstat<0)/2850;
[h,~,~,s]=ttest(fc_all_mod(:,:,3,2),fc_all_mod(:,:,1,2),'dim',2);
dpz_p(2) = sum(h>0&s.tstat>0)/2850;
dpz_n(2) = sum(h>0&s.tstat<0)/2850;

% figure_w;
% subplot(2,3,5); hold on
%   
% bar(1,atx_p(1),'facecolor',[.8 .3 .3])
% bar(2,atx_n(1),'facecolor',[0 0 1])
% bar(3,atx_p(2),'facecolor',[1 0 0])
% bar(4,atx_n(2),'facecolor',[.3 .3 1])
%   
% axis square; tp_editplots
% axis([0 5 0 0.3]);
% set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Pos. (Rest)';'Neg. (Rest)';'Pos. (Task)';'Neg. (Task)'})
% xtickangle(gca,45)
% 
% subplot(2,3,6); hold on
%   
% bar(1,dpz_p(1),'facecolor',[.8 .3 .3])
% bar(2,dpz_n(1),'facecolor',[0 0 1])
% bar(3,dpz_p(2),'facecolor',[1 0 0])
% bar(4,dpz_n(2),'facecolor',[.3 .3 1])
%   
% axis square; tp_editplots
% axis([0 5 0 0.3]);
% set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Pos. (Rest)';'Neg. (Rest)';'Pos. (Task)';'Neg. (Task)'})
% xtickangle(gca,45)
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
gain = 1;
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
imagesc((p_pos_atx_sim(:,:,gain,1)-p_pos_atx_sim(:,:,gain,2)),[-0.2 0.2]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Ctx - Pos (g=%.2f)',gg(gain)')); colormap(gca,redblue)

subplot(2,3,6); hold on
imagesc(p_neg_atx_sim(:,:,gain,1)-p_neg_atx_sim(:,:,gain,2),[-0.2 0.2]); set(gca,'ydir','normal'); axis tight square; tp_editplots
set(gca,'xtick',[1:6:37],'xticklabel',num2cell(es([1:6:37])*mean(diff(Ies))))
set(gca,'ytick',[1:6:37],'yticklabel',num2cell(is([1:6:37])*mean(diff(Ies))))
title(sprintf('Ctx - Neg (g=%.2f)',gg(gain)'));colormap(gca,redblue)
% 
% for i = 1 : length(es)
%   neg(i) = p_neg_atx_sim(i,i,gain,1)-p_neg_atx_sim(i,i,gain,2);
%   pos(i) = p_pos_atx_sim(i,i,gain,1)-p_pos_atx_sim(i,i,gain,2);
% end
% 
% figure; set(gcf,'color','w');
% subplot(2,3,3); hold on
% plot(pos,'r')
% plot(neg,'b')
% axis square; tp_editplots
% axis([-1 38 -0.2 0.2]);
% xlabel('\Delta(Background input)'); ylabel('Fraction of altered correlations [%]')
% % axis([])



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

%% PLOT KS DISTANCE FROM FCD MATRICES (FIRING RATES AND ENVELOPES)
osc = osc1(1:2:end,1:2:end,3:2:31);

for iG = 1 : length(Gg)
  
  par1=nanmean(outp.fc(:,:,iG,:),4);
  par1(osc(:,:,iG)==1)=nan;
  
  figure_w;
  
  subplot(1,4,1);
  b= imagesc(par1,[0 0.1]); 
  set(gca,'ydir','normal'); axis square;  colormap(plasma);  set(b,'AlphaData',~isnan(par1))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  set(gca,'YTick',1:20:length(Iis),'YTickLabels',num2cell(Iis(1:20:end)))
  title('Mean FC (Envelopes)'); tp_editplots
  
  par1=nanmean(outp.ks_rest_env(:,:,iG,:),4);
  par2=nanmean(outp.ks_rest_fr(:,:,iG,:),4);
  
  par1((osc(:,:,iG)==1))=nan;
  par2((osc(:,:,iG)==1))=nan;
  
  k=subplot(1,4,3);
  b=imagesc(par1,[0.4 0.8]); set(gca,'ydir','normal'); set(b,'AlphaData',~isnan(par1))
  set(gca,'ydir','normal'); axis square;  colormap(k,plasma);
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  set(gca,'YTick',1:20:length(Iis),'YTickLabels',num2cell(Iis(1:20:end)))
  title('KS-D (Envelopes)'); tp_editplots
  
  k=subplot(1,4,4);
  b=imagesc(par2,[0.4 0.8]); set(gca,'ydir','normal');
  set(gca,'ydir','normal'); axis square;  colormap(k,plasma);  set(b,'AlphaData',~isnan(par2))
  set(gca,'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
  set(gca,'YTick',1:20:length(Iis),'YTickLabels',num2cell(Iis(1:20:end)))
  title('KS-D (Firing rates)'); tp_editplots
end
  
  
 
 
  
  
  
