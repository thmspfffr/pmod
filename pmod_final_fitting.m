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
% VERSION 3: 27/01/2020
% %-------------------------------------------------------------------------
v           = 3;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = [1.2:-0.01:1.10];
Gains       = [-0.1:0.02:0.12]; 
nTrials     = 1;
tmax        = 6500;  % in units of tauE
dt          = 0.01;
%-------------------------------------------------------------------------

v_sim = v;
% connectivity, AAL
v_conn =  33;

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

load ~/M.mat

clear outp
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
addpath ~/pupmod/matlab
cleandat = pupmod_loadpowcorr(v_conn,SUBJLIST,1);

mask = logical(tril(ones(76,76),-1));

% transform avg fc matrices to AAL BCN
k = 1 : 90;

% exclude subcortical regions
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

% load SC matrix, exclude subcortical regions
load ~/sc90.mat
SC = SC(include_bcn,include_bcn);

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
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,1,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); 
  tmp = tmp(include_bcn,include_bcn);
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,2,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  fc_task_indiv(:,isubj) = tmp(mask);
  
  for im = 1 : 3
    for icond = 1 : 2
      tmp = squeeze(nanmean(cleandat(:,:,isubj,im,icond,[3 4 5 6]),6));
      tmp = tp_match_aal(para,tmp); 
      tmp = tmp(include_bcn,include_bcn);
      fc(:,isubj,im,icond)=tmp(mask);
    end
  end
end

% end
k = [0 0];

%% CORRELATE FC_SIM WITH FC_EMP

if ~exist(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  
  for iies = 1 : length(Ies)
    iies
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        for iiis = 1 :length(Iis)
          
          load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,iies,iiis,iG,igain,v_sim))
          
          outp.fc_sim_fr_tmp = out.fc_FR;
          
          [outp.r_fr_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest(mask));          
          [outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest_indiv);
          outp.dist_fr_rest(iies,iiis,iG,igain) = 1-(outp.r_fr_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_fr_tmp(mask))).^2);  
          outp.dist_fr_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          outp.fc_sim_fr_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_fr_tmp(mask));
          
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          if isempty( out.peakfreq)
             out.peakfreq = nan;
          end
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          outp.alphapow(iies,iiis,iG,igain) = mean(out.alphapow);   

        end
      end
    end
  end
   
  save(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim),'outp')
 
else
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_all_v%d.mat',v_sim,v_sim))
  load(sprintf('~/pmod/proc/detosc/v%d/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim,v_sim))
end

error('!')

%% PEAK FREQUENCY

igain = find(Gains==0);
iG = find(Gg==1.15);

figure_w
ax{1} = subplot(2,3,1);
par=outp.peakfreq(:,:,iG,igain);
par(osc1(:,:,iG,igain)>0)=nan;
b=imagesc(par,[5 20]); colormap(plasma);
set(b,'AlphaData',~isnan(par)); set(gca,'ydir','normal')

ax{2} = subplot(2,3,2);
b=imagesc(par>7.8814 & par<18.7452)
set(b,'AlphaData',~isnan(par)); set(gca,'ydir','normal')

for iax = [1 2]
  colormap(ax{iax},plasma)
  xlabel(ax{iax},'Background input to I')
  ylabel(ax{iax},'Background input to E');
  set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
  set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  tp_editplots(ax{iax})
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  axis(ax{iax},'square')
end

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_peakfreq_ig%d_igain%d_v%d.eps',iG,igain,v_sim))

%%
lim = 121;

oscthres = 0;
nancol = [0.97 0.97 0.97];

clear par ax
figure; set(gcf,'color','w')

igain = find(Gains==0);
iG = find(Gg==1.15);
osc = osc1(1:lim,1:lim,iG,igain);

% plot peak freq model
ax{1} = subplot(3,2,1); hold on
par = outp.fc_sim_fr_mean(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.02]); set(b,'AlphaData',~isnan(par))

% plot corr fc_model with fc_meg
ax{2} = subplot(3,2,3); hold on
par = outp.r_fr_rest_corr(1:lim,1:lim,iG,igain);
par(osc==1)=nan;
b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isnan(par))

for iax = 1 :2 
  colormap(ax{iax},plasma)
  xlabel(ax{iax},'Background input to I')
  ylabel(ax{iax},'Background input to E');
  set(ax{iax},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
  set(ax{iax},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  tp_editplots(ax{iax})
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  axis(ax{iax},'square')
end

set(gcf,'Renderer','Painters')
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_fc_frandenv_iG%d_v%d.eps',iG,v_sim))

%% FIT RESTING STATE PLACEBO 
% (THROUGH DISTANCE BETWEEN FC_SIM And FC_EMP)
% ----------------------------------

clear idx par ax cnt idx_rest par
lim = 121;
iG         = find(Gg==1.15); 
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
    par(nanmean(nanmean(osc1(:,:,iG,:),3),4)>0)  = nan;
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

% find idx of delta_gain = 0 and coupling (Gg) = 1.15
baseline_gain=find(Gains==0);
baseline_coupling = find(Gg==1.15);

SUBJ = 1:28;
clear fc_mod_rest fc_mod_task

igain = find(Gains==0);

h = figure_w; 

subplot(2,3,1); hold on
par = outp.fc_sim_fr_mean(:,:,baseline_coupling,baseline_gain);
par(osc1(:,:,baseline_coupling,baseline_gain)==1)=nan; 

for isubj = 1:size(Mu,1)
  scatter(Mu(isubj,1),Mu(isubj,2),10,'o','markerfacecolor',cols(c(isubj),:),'markeredgecolor','none');
end

line([34 34],[37 45],'color','k')
line([34 46],[45 45],'color','k')
line([46 46],[37 45],'color','k')
line([34 46],[37 37],'color','k')

axis(gca,[1 121 1 121])
tp_editplots; axis square
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel('Background input to I'); ylabel('Background input to E'); 

subplot(2,3,4); hold on
for isubj = 1:size(Mu,1)
  scatter(Mu(isubj,1),Mu(isubj,2),20,'o','markerfacecolor',cols(c(isubj),:),'markeredgecolor','none');
end

axis(gca,[34 46 37 45])
tp_editplots; %axis square
set(gca,'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
set(gca,'XTick',1:5:length(Iis),'XTickLabels',num2cell(Iis(1:5:end)))
xlabel('Background input to I'); ylabel('Background input to E'); 

load redblue.mat

% ATX: 
igain = 11; iG = 6;
% DPZ: 
% igain = 8; iG = 10
  
for isubj = SUBJ
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),baseline_coupling,baseline_gain,v_sim))
  fc_mod_rest(isubj,1) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,1,1) = out.fc_FR(mask);
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,igain,v_sim))
  fc_mod_rest(isubj,2) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,2,1) = out.fc_FR(mask);
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+9,idx_rest.inh(isubj)+15,baseline_coupling,baseline_gain,v_sim))
  fc_mod_task(isubj,1) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,1,2) = out.fc_FR(mask);
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+9,idx_rest.inh(isubj)+15,iG,igain,v_sim))
  fc_mod_task(isubj,2) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,2,2) = out.fc_FR(mask);
  
end

par = outp.fc_sim_fr_mean(:,:,iG,igain)-outp.fc_sim_fr_mean(:,:,baseline_coupling,baseline_gain);
par(osc1(:,:,baseline_coupling,baseline_gain)>0)=nan;
par(osc1(:,:,iG,igain)>0)=nan;

subplot(2,3,2); hold on
b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');
text(10,110,['\DeltaGain' sprintf('=%.2f',Gains(igain))],'fontsize',7)
text(10,100,['\DeltaCoupl.' sprintf('=%.2f',Gg(iG)-Gg(baseline_coupling))],'fontsize',7)

scatter(round(mean(idx_rest.inh(SUBJ))),round(mean(idx_rest.exc(SUBJ))),20,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');
scatter(round(mean(idx_rest.inh(SUBJ)))+15,round(mean(idx_rest.exc(SUBJ)))+9,20,'o','markerfacecolor','y','markeredgecolor','k');

set(b,'AlphaData',~isnan(par));
colormap(gca,redblue)
ylabel(gca,'Background input to E');
xlabel(gca,'Background input to I');
set(gca,'YTick',1:40:121,'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:121,'XTickLabels',num2cell(Iis(1:40:end)))

tp_editplots; axis square; axis([1 121 1 121])
  

% DPZ
igain = 8; iG = 11

for isubj = SUBJ
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,igain,v_sim))
  fc_mod_rest(isubj,3) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,3,1) = out.fc_FR(mask);
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+9,idx_rest.inh(isubj)+14,iG,igain,v_sim))
  fc_mod_task(isubj,3) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,3,2) = out.fc_FR(mask);
  
end

par = outp.fc_sim_fr_mean(:,:,iG,igain)-outp.fc_sim_fr_mean(:,:,baseline_coupling,baseline_gain);
par(osc1(:,:,baseline_coupling,baseline_gain)>0)=nan;
par(osc1(:,:,iG,igain)==1)=nan;

subplot(2,3,5); hold on
b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');
text(10,110,['\DeltaGain' sprintf('=%.2f',Gains(igain))],'fontsize',7)
text(10,100,['\DeltaCoupl.' sprintf('=%.2f',Gg(iG)-Gg(baseline_coupling))],'fontsize',7)

scatter(round(mean(idx_rest.inh(SUBJ))),round(mean(idx_rest.exc(SUBJ))),20,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');
scatter(round(mean(idx_rest.inh(SUBJ)))+14,round(mean(idx_rest.exc(SUBJ)))+9,20,'o','markerfacecolor','y','markeredgecolor','k');

set(b,'AlphaData',~isnan(par));
colormap(gca,redblue);
ylabel(gca,'Background input to E');
xlabel(gca,'Background input to I');
set(gca,'YTick',1:40:121,'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:121,'XTickLabels',num2cell(Iis(1:40:end)))
tp_editplots; axis square; axis([1 121 1 121])

subplot(2,3,3); hold on

prc_change_atx_rest =  100*(fc_mod_rest(:,2)-fc_mod_rest(:,1))./fc_mod_rest(:,1);
prc_change_atx_task =  100*(fc_mod_task(:,2)-fc_mod_task(:,1))./fc_mod_task(:,1);
bar(2,mean(prc_change_atx_task),'edgecolor',[0 0 0],'facecolor',[1 1 0])
bar(1,mean(prc_change_atx_rest),'edgecolor',[0 0 0],'facecolor',[.95 .95 .95])
ylabel('Change in correlation [%]')

[~,p_task_atx]=ttest(prc_change_atx_task);
[~,p_rest_atx]=ttest(prc_change_atx_rest);

if p_rest_atx<0.001
  text(1,25,'***','horizontalalignment','center','verticalalignment','bottom')
elseif p_rest_atx <0.01
  text(1,25,'**','horizontalalignment','center','verticalalignment','bottom')
elseif p_rest_atx < 0.05
  text(1,25,'*','horizontalalignment','center','verticalalignment','bottom')
else 
  text(1,25,'n.s.','horizontalalignment','center','verticalalignment','bottom')
end
  
if p_task_atx<0.001
  text(2,25,'***','horizontalalignment','center','verticalalignment','bottom')
elseif p_task_atx <0.01
  text(2,25,'**','horizontalalignment','center','verticalalignment','bottom')
elseif p_task_atx < 0.05
  text(2,25,'*','horizontalalignment','center','verticalalignment','bottom')
else 
  text(2,25,'n.s.','horizontalalignment','center','verticalalignment','bottom')
end
  
axis square; tp_editplots
axis([0 3 -30 30]);
set(gca,'XTick',[1 2],'XTickLabels',{'Rest';'Task'})
set(gca,'YTick',[-30:10:30],'YTickLabels',[-30:10:30])

xtickangle(gca,45)

subplot(2,3,6); hold on

prc_change_dpz_rest =  100*(fc_mod_rest(:,3)-fc_mod_rest(:,1))./fc_mod_rest(:,1);
prc_change_dpz_task =  100*(fc_mod_task(:,3)-fc_mod_task(:,1))./fc_mod_task(:,1);
bar(1,mean(prc_change_dpz_rest),'edgecolor',[0 0 0],'facecolor',[0.95 0.95 0.95])
bar(2,mean(prc_change_dpz_task),'edgecolor',[0 0 0],'facecolor',[1 1 0])

[~,p_task_dpz]=ttest(prc_change_dpz_task);
[~,p_rest_dpz]=ttest(prc_change_dpz_rest);

if p_rest_dpz<0.001
  text(1,25,'***','horizontalalignment','center','verticalalignment','bottom')
elseif p_rest_dpz <0.01
  text(1,25,'**','horizontalalignment','center','verticalalignment','bottom')
elseif p_rest_dpz < 0.05
  text(1,25,'*','horizontalalignment','center','verticalalignment','bottom')
else 
  text(1,25,'n.s.','horizontalalignment','center','verticalalignment','bottom')
end
  
if p_task_dpz<0.001
  text(2,25,'***','horizontalalignment','center','verticalalignment','bottom')
elseif p_task_dpz <0.01
  text(2,25,'**','horizontalalignment','center','verticalalignment','bottom')
elseif p_task_dpz < 0.05
  text(2,25,'*','horizontalalignment','center','verticalalignment','bottom')
else 
  text(2,25,'n.s.','horizontalalignment','center','verticalalignment','bottom')
end

axis square; tp_editplots
axis([0 3 -30 30]);
set(gca,'XTick',[1 2],'XTickLabels',{'Rest';'Task'})
xtickangle(gca,45)
ylabel('Change in correlation [%]')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_highvslowgain_ig%d_gain%d_v%d.pdf',iG,igain,v_sim))

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

