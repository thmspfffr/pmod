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

fc_rest = squeeze(mean(mean(cleandat(1:90,1:90,:,1,1,[3 4 5 6]),3),6));
fc_task = squeeze(mean(mean(cleandat(1:90,1:90,:,1,2,[3 4 5 6]),3),6));

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
  tmp = squeeze(mean(cleandat(:,:,isubj,1,1,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); 
  tmp = tmp(include_bcn,include_bcn);
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(mean(cleandat(:,:,isubj,1,2,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  fc_task_indiv(:,isubj) = tmp(mask);
  
  for im = 1 : 3
    for icond = 1 : 2
      tmp = squeeze(mean(cleandat(:,:,isubj,im,icond,[3 4 5 6]),6));
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
             out.peakfreq = gla;
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
par(osc1(:,:,iG,igain)>0)=gla;
b=imagesc(par,[5 20]); colormap(plasma);
set(b,'AlphaData',~isgla(par)); set(gca,'ydir','normal')

ax{2} = subplot(2,3,2);
b=imagesc(par>7.8814 & par<18.7452)
set(b,'AlphaData',~isgla(par)); set(gca,'ydir','normal')

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
glacol = [0.97 0.97 0.97];

clear par ax
figure; set(gcf,'color','w')

% ---------
% Basline: delta_gain = 0 
% ---------
igain = find(Gains==0);
% ---------
% Set global coupling to 1.15 (see pmod_fitting_globalcoupling.m)
% ---------
iG = find(Gg==1.15);
% ---------
osc = osc1(1:lim,1:lim,iG,igain);

% plot peak freq model
ax{1} = subplot(3,2,1); hold on
par = outp.fc_sim_fr_mean(1:lim,1:lim,iG,igain);
par(osc==1)=gla;
b=imagesc(par,[0 0.02]); set(b,'AlphaData',~isgla(par))

% plot corr fc_model with fc_meg
ax{2} = subplot(3,2,3); hold on
par = outp.r_fr_rest_corr(1:lim,1:lim,iG,igain);
par(osc==1)=gla;
b=imagesc(par,[0 0.5]); set(b,'AlphaData',~isgla(par))

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
% ---------
% Global coupling = 1.15 (baseline)
% ---------
iG         = find(Gg==1.15); 
% ---------
% Basline: delta_gain = 0 
% ---------
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
    par(mean(mean(osc1(:,:,iG,:),3),4)>0)  = gla;
    binmap              = par<prctile(par(:),prc_thresh);
    bw                  = bwlabel(binmap,8);
    cnt                 = [];
    
    for iclust = 1 : max(bw(:))
      cnt(iclust) = sum(bw(:)==iclust);
    end
    
    % if no cluster(s), omit subject from further analysis
    if isempty(cnt)
      idx_rest.inh(isubj) = gla;
      idx_rest.exc(isubj) = gla;
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
      par1(par1==0)=gla;
      imagesc(par1)
      b=imagesc(par1,[0 0.4]); set(b,'AlphaData',~isgla(par1))
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
  par(osc>oscthresh)  = gla;
  b=imagesc(par,[0.6 1.4]); colormap(plasma)
  set(b,'AlphaData',~isgla(par))
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

% task_exc = 8;
% task_inh = 15;

task_exc = 10;
task_inh = 17;

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
par(osc1(:,:,baseline_coupling,baseline_gain)==1)=gla; 

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
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+task_exc,idx_rest.inh(isubj)+task_inh,baseline_coupling,baseline_gain,v_sim))
  fc_mod_task(isubj,1) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,1,2) = out.fc_FR(mask);
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+task_exc,idx_rest.inh(isubj)+task_inh,iG,igain,v_sim))
  fc_mod_task(isubj,2) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,2,2) = out.fc_FR(mask);
  
end

par = outp.fc_sim_fr_mean(:,:,iG,igain)-outp.fc_sim_fr_mean(:,:,baseline_coupling,baseline_gain);
par(osc1(:,:,baseline_coupling,baseline_gain)>0)=gla;
par(osc1(:,:,iG,igain)>0)=gla;

subplot(2,3,2); hold on
b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');
text(10,110,['\DeltaGain' sprintf('=%.2f',Gains(igain))],'fontsize',7)
text(10,100,['\DeltaCoupl.' sprintf('=%.2f',Gg(iG)-Gg(baseline_coupling))],'fontsize',7)

scatter(round(mean(idx_rest.inh(SUBJ))),round(mean(idx_rest.exc(SUBJ))),20,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');
scatter(round(mean(idx_rest.inh(SUBJ)))+task_inh,round(mean(idx_rest.exc(SUBJ)))+task_exc,20,'o','markerfacecolor','y','markeredgecolor','k');

set(b,'AlphaData',~isgla(par));
colormap(gca,redblue)
ylabel(gca,'Background input to E');
xlabel(gca,'Background input to I');
set(gca,'YTick',1:40:121,'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:121,'XTickLabels',num2cell(Iis(1:40:end)))

tp_editplots; axis square; axis([1 121 1 121])

% DPZ
igain =8; iG = 10;
% igain = 11; iG = 10

for isubj = SUBJ
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,igain,v_sim))
  fc_mod_rest(isubj,3) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,3,1) = out.fc_FR(mask);
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+task_exc,idx_rest.inh(isubj)+task_inh,iG,igain,v_sim))
  fc_mod_task(isubj,3) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,3,2) = out.fc_FR(mask);
  
end

par = outp.fc_sim_fr_mean(:,:,iG,igain)-outp.fc_sim_fr_mean(:,:,baseline_coupling,baseline_gain);
par(osc1(:,:,baseline_coupling,baseline_gain)>0)=gla;
par(osc1(:,:,iG,igain)==1)=gla;

subplot(2,3,5); hold on
b=imagesc(par,[-0.01 0.01]); set(gca,'ydir','normal');
text(10,110,['\DeltaGain' sprintf('=%.2f',Gains(igain))],'fontsize',7)
text(10,100,['\DeltaCoupl.' sprintf('=%.2f',Gg(iG)-Gg(baseline_coupling))],'fontsize',7)

scatter(round(mean(idx_rest.inh(SUBJ))),round(mean(idx_rest.exc(SUBJ))),20,'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k');
scatter(round(mean(idx_rest.inh(SUBJ)))+task_inh,round(mean(idx_rest.exc(SUBJ)))+task_exc,20,'o','markerfacecolor','y','markeredgecolor','k');

set(b,'AlphaData',~isgla(par));
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

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_highvslowgain_ig%d_gain%d_tE%dtI%d_v%d.pdf',iG,igain,task_exc,task_inh,v_sim))

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

%% CHECK ALTERNATIVE PARAMETERS AS REQUESTED BY REV2
% --------------------------
% 
% --------------------------

% ATX: 
% igain = 6; iG = 6;
% DPZ: 
% igain = 8; iG = 10

task_mod = [-9 -6 -3 0 3 6 9]
  
for isubj = SUBJ
  for iG = 1 : length(Gg)
    for igain = 1 : length(Gains)
      for itask = 1 :length(task_mod)
  isubj
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),baseline_coupling,baseline_gain,v_sim))
  fc_mod_rest(isubj,iG,igain,itask,1) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,iG,igain,itask,1,1) = out.fc_FR(mask);
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj),idx_rest.inh(isubj),iG,igain,v_sim))
  fc_mod_rest(isubj,iG,igain,itask,2) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,iG,igain,itask,2,1) = out.fc_FR(mask);
  
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+task_exc+task_mod(itask),idx_rest.inh(isubj)+task_inh+task_mod(itask),baseline_coupling,baseline_gain,v_sim))
  fc_mod_task(isubj,iG,igain,itask,1) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,iG,igain,itask,1,2) = out.fc_FR(mask);
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,idx_rest.exc(isubj)+task_exc+task_mod(itask),idx_rest.inh(isubj)+task_inh+task_mod(itask),iG,igain,v_sim))
  fc_mod_task(isubj,iG,igain,itask,2) = mean(out.fc_FR(mask));
  fc_all_mod(:,isubj,iG,igain,itask,2,2) = out.fc_FR(mask);
      end
    end
  end
end


%%
[h,~,~,s]=ttest(fc_all_mod(:,:,:,:,2,1),fc_all_mod(:,:,:,:,1,1),'dim',2);
pos_rest = 100*squeeze(sum(squeeze(h) > 0 & squeeze(s.tstat) > 0))./sum(mask(:));
neg_rest = 100*squeeze(sum(squeeze(h) > 0 & squeeze(s.tstat) < 0))./sum(mask(:));

[h,~,~,s]=ttest(fc_all_mod(:,:,:,:,2,2),fc_all_mod(:,:,:,:,1,2),'dim',2);
pos_task = 100*squeeze(sum(squeeze(h) > 0 & squeeze(s.tstat) > 0))./sum(mask(:));
neg_task = 100*squeeze(sum(squeeze(h) > 0 & squeeze(s.tstat) < 0))./sum(mask(:));

cmap1 = cmap(1:125,:);
cmap2 = cmap(127:end,:);

figure_w
k1 = subplot(2,2,1); 
imagesc(neg_rest,[0 30]);
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6)
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6)
colormap(k1,cmap1)

k2=subplot(2,2,2);
imagesc(pos_rest,[0 30])
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6)
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6)
colormap(k2,cmap2)

k3=subplot(2,2,3); 
imagesc(neg_task,[0 30]);
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6)
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6)
colormap(k3,cmap1)

k4= subplot(2,2,4);
imagesc(pos_task,[0 30])
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6)
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6)
colormap(k4,cmap2)

figure_w
k1 = subplot(2,2,1); 
imagesc(neg_rest-neg_task,[-20 20]);
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); ylabel('Change in coupling')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); xlabel('Change in gain')
colormap(k1,cmap)

k2=subplot(2,2,2);
imagesc(pos_rest-pos_task,[-20 20])
tp_editplots; axis square
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')
colormap(k2,cmap)

%%

itask = 8;

igain_atx = 11; iG_atx = 6
igain_dpz = 8; iG_dpz = 10

figure_w
subplot(2,2,1); hold on
h=squeeze(ttest(fc_mod_rest(:,:,:,itask,2),fc_mod_rest(:,:,:,itask,1),'dim',1));
par_rest = squeeze(mean(fc_mod_rest(:,:,:,itask,2))-mean(fc_mod_rest(:,:,:,itask,1))); par_rest(isgla(par_rest))=0; h(isgla(h))=0;
imagesc(par_rest.*squeeze(h),[-0.01 0.01])
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')


subplot(2,2,2); hold on
h=squeeze(ttest(fc_mod_task(:,:,:,3,2),fc_mod_task(:,:,:,itask,1),'dim',1));
par_task = squeeze(mean(fc_mod_task(:,:,:,itask,2))-mean(fc_mod_task(:,:,:,itask,1))); par_task(isgla(par_task))=0; h(isgla(h))=0;
imagesc(par_task,[-0.01 0.01]); 
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')


subplot(2,2,3); hold on
imagesc(par_rest-par_task,[-0.01 0.01]); par(isgla(par))=0; h(isgla(h))=0;
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')


atx_bin = (par_task > 0.0015) & abs(par_rest)<0.0015 & (par_rest-par_task)<-0.0015
dpz_bin = (abs(par_task) < 0.0015) & par_rest<-0.0015 & (par_rest-par_task)<-0.0015
dpz_bin = -1*dpz_bin;
par = atx_bin+dpz_bin;

subplot(2,2,4); hold on
imagesc(par,[-0.01 0.01]); par(isgla(par))=0; h(isgla(h))=0;
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_alternativeparams_task%d_v%d.pdf',itask,v))


% a1=squeeze(mean(mean(mean(cleandat(:,:,:,1,2,3:6),6)-mean(cleandat(:,:,:,1,1,3:6),6),1),2))

% a=mean(mean(all_pow(:,:,1,6:9,2))-mean(all_pow(:,:,1,6:9,1)),4)
% b1=squeeze(mean(mean(mean(cleandat(:,:,:,2,2,3:6),6)-mean(cleandat(:,:,:,1,2,3:6),6),1),2))

%%

for i = 1 : 24
  
  fc1(:,:,i)= corr(squeeze(M(i,3,:,:)),squeeze(M(i,3,:,:)));
  fc2(:,:,i)= corr(squeeze(M(i,4,:,:)),squeeze(M(i,4,:,:)));
  
end

mask = logical(tril(ones(90,90),-1));

[h,~,~,s]=ttest(fc2,fc1,'dim',3,'alpha',0.05);

pos = sum(h(mask) > 0 & s.tstat(mask) > 0)/sum(mask(:));
neg = sum(h(mask) > 0 & s.tstat(mask) < 0)/sum(mask(:));


%% PLOT KURAMOTO ORDER PARAMETER
    





