%% FITTING
% pmod_final_fitting

clear
% %-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 0:0.1:3;
% Gains       = 0; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
% %-------------------------------------------------------------------------
% VERSION 2: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.025:0;
Iis         = -5:0.025:-1;
Gg          = 1.7;
Gains       = [0 0.025:0.025:0.4 -0.025:-0.025:-0.4 0.425:0.025:0.6  0.625:0.025:0.7]; 
nTrials     = 1;
tmax        = 6500;  % in units of tauE
EC          = 0;
%--------------------------------------------------------------------------
% VERSION 22: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
% v           = 22;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 0.85;
% Gains       = [0 0.025:0.025:0.4 -0.025:-0.025:-0.4]; 
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
%--------------------------------------------------------------------------


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
cleandat = pupmod_loadpowcorr(v_conn,1);

% load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat

% peakfreq_rest = m_res(1);
% peakfreq_task = m_tsk(1);

load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

mask = logical(tril(ones(76,76),-1));

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

% LOAD SC MATRIX
load ~/sc90.mat
SC = SC(include_bcn,include_bcn);

ifoi = 6;
fc_rest = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,1,ifoi),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,2,ifoi),3),6));

para = [];
para.transfer = 'to_bcn';
para.N = 90;

fc_rest = tp_match_aal(para,fc_rest);
fc_task = tp_match_aal(para,fc_task);

fc_rest = fc_rest(include_bcn,include_bcn);
fc_task = fc_task(include_bcn,include_bcn);

% transform indiv. subj. matrices to AAL BCN 
for isubj = 1 : 28
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,1,ifoi),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn); 
  [corrwithfc_rest(isubj), p_corrwithfc_rest(isubj)]  = corr(tmp(mask),SC(mask));
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,2,ifoi),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn); 
  fc_task_indiv(:,isubj) = tmp(mask);
  [corrwithfc_task(isubj), p_corrwithfc_task(isubj)] = corr(tmp(mask),SC(mask));
end

k = [0 0];

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
% 
  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))
          
          outp.fc_sim_env_tmp = out.FC_env;

          [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
          [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
          
          [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv);
          [outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_task_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task_indiv);
             
          outp.dist_rest(iies,iiis,iG,igain) = 1-(outp.r_env_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);
          outp.dist_task(iies,iiis,iG,igain) = 1-(outp.r_env_task_corr(iies,iiis,iG,igain)-(mean(fc_task(mask))-mean(outp.fc_sim_env_tmp(mask))).^2);

          outp.dist_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
          outp.dist_task_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_env_task_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_task_indiv))-mean(outp.fc_sim_env_tmp(mask))).^2);
          
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          outp.alphapow(iies,iiis,iG,igain) = mean(out.alphapow);

          fclose all;
          
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
end

clear cleandat


error('!')
%% PLOT GAIN FUNCTION

aE = 1;
x = -10 : 0.01 : 10;
y =  1./(1 + exp(-x/aE) );

figure; set(gcf,'color','white');
subplot(2,2,1); hold on
plot(x,y,'k'); 
xlabel('Input'); ylabel('Output'); tp_editplots
axis([-12 12 -0.05 1.05]);

aE = 2;
y =  1./(1 + exp(-x/aE) );
plot(x,y,'b'); 

aE = 0.5;
y =  1./(1 + exp(-x/aE) );
plot(x,y,'r'); 

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_gainfun_v%d.eps',v_sim))

%% DETERMINE GLOBAL COUPLING PARAMETER FIRST

osc = osc1;
oscthresh = 0;

figure; set (gcf,'color','w');

subplot(2,4,1); hold on
par = outp.dist;
par1 = squeeze(nanmean(nanmean(par)));
plot(par1,'linewidth',2);
[i,j]=min(par1);
line([j j],[0.75 1.25],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[1 1],'linestyle','-','color',[0 0 0])
box on
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0.8:0.1:1.2,'YTickLabel',0.8:0.1:1.2);
xlabel('Global coupling'); ylabel('Distance')
axis([1 31 0.8 1.2])
axis(gca,'square')
tp_editplots

subplot(2,4,2); hold on
par = outp.r_env_rest_corr;
par1 = squeeze(nanmean(nanmean(par)));
plot(par1,'linewidth',2);
[i,j]=max(par1);
line([j j],[-0.05 0.2],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
box on
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',-0.05:0.05:0.2,'YTickLabel',-0.05:0.05:0.2);
xlabel('Global coupling'); ylabel('Mean correlation')
axis([1 31 -0.05 0.2])
axis(gca,'square')
tp_editplots


subplot(2,4,3)
par = outp.r_env_rest_corr;
par1 = squeeze(nansum(nansum(par)))
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-500 2500],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',-500:500:2500,'YTickLabel',-500:500:2500);
xlabel('Global coupling'); ylabel('Summed correlation')
axis([1 31 -500 2500])
axis(gca,'square')
tp_editplots

subplot(2,4,4)
par = outp.r_env_rest_corr;
par1 = squeeze(max(max(par)))
% par(osc1>oscthresh) = nan;
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-50 200],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0:0.2:1,'YTickLabel',0:.2:1);
xlabel('Global coupling'); ylabel('Max correlation')
axis([1 31 0 1])
axis(gca,'square')
tp_editplots
% -------------------
% SAME PLOTS AS 1-3, BUT SUSTAINED OSCILLATIONS EXCLUDED
% -------------------

subplot(2,4,5)
par = outp.dist;
par(osc1>oscthresh) = nan;
par1 = squeeze(nanmean(nanmean(par)))
[~,idx]=min(par1);
plot(par1,'linewidth',2);
line([idx idx],[0.975 1.025],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[1 1],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0.975:0.0125:1.025,'YTickLabel',0.975:0.0125:1.025);
xlabel('Global coupling'); ylabel('Distance)')
axis([1 31 0.975 1.025])
axis(gca,'square')
tp_editplots

subplot(2,4,6)
par = outp.r_env_rest_corr;
par(osc1>oscthresh) = nan;
par1 = squeeze(nanmean(nanmean(par)))
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-0.05 0.5],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',-0.01:0.01:0.04,'YTickLabel',-0.01:0.01:0.04);
xlabel('Global coupling'); ylabel('Mean Correlation')
axis([1 31 -0.01 0.04])
axis(gca,'square')
tp_editplots

subplot(2,4,7)
par = outp.r_env_rest_corr;
par(osc1>oscthresh) = nan;
par1 = squeeze(nansum(nansum(par)))
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-50 200],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',-50:50:200,'YTickLabel',-50:50:200);
xlabel('Global coupling'); ylabel('Summed correlation')
axis([1 31 -50 200])
axis(gca,'square')
tp_editplots

subplot(2,4,8)
par = outp.r_env_rest_corr;
par(osc1>oscthresh) = nan;
par1 = squeeze(max(max(par)))
[~,idx]=max(par1);
plot(par1,'linewidth',2);
line([idx idx],[-50 200],'linestyle',':','color',[0.3 0.3 0.3])
line([1 31],[0 0],'linestyle','-','color',[0 0 0])
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
set(gca,'YTick',0:0.2:1,'YTickLabel',0:.2:1);
xlabel('Global coupling'); ylabel('Max correlation')
axis([1 31 0 1])
axis(gca,'square')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_fitG_v%d.pdf',v_sim))


%%
osc = osc1(:,:,:,1);
oscthres = 0;
nancol = [0.97 0.97 0.97];

clear par ax
figure; set(gcf,'color','w')

igain = 1;
iG = 1;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = outp.fc_sim_env_mean(:,:,iG,igain);
par(osc>oscthres)=nan;
imagescnan(par,[0 0.02],'NanColor',nancol)
title('Mean FC (Envelopes)');
line([find(Iis==-4) find(Iis==-4)],[find(Ies==-3.2) find(Ies==-2.2)],'color','k')
line([find(Iis==-4) find(Iis==-2.8)],[find(Ies==-2.2) find(Ies==-2.2)],'color','k')
line([find(Iis==-2.8) find(Iis==-2.8)],[find(Ies==-3.2) find(Ies==-2.2)],'color','k')
line([find(Iis==-2.8) find(Iis==-4)],[find(Ies==-3.2) find(Ies==-3.2)],'color','k')

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par= outp.fc_sim_env_mean(:,:,iG,igain);
par(osc>oscthres)=nan;
par = par([find(Ies==-3.2):find(Ies==-2.2)],[find(Iis==-4):find(Iis==-2.8)]);
imagescnan(par,[0 0.02],'NanColor',nancol)
title('Mean FC (Envelopes)');
axis off; colormap(plasma); axis equal tight

% plot peak freq model
% ax{2} = subplot(2,2,2); hold on
% par = peakfreq_rest-squeeze(outp.peakfreq(:,:,iG,igain));
% par(abs(par)<1)=2; par(abs(par)>2)=0; 
% par(osc>oscthres)=0;
% m1 = par>0;
% par = double(par); par(par<1)=NaN;
% imagescnan(par,[-1 1],'NanColor',nancol)
% title('Peak freq: Masked');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = -log10(squeeze(outp.p_env_rest_corr(:,:,iG,igain)));
par(osc>oscthres)=nan;
imagescnan(par,[0 3],'NanColor',nancol)
title('r(FC_{sim},FC_{MEG})');

% plot peak freq model
% find significant correlations
fdr_p = fdr1(reshape(squeeze(outp.p_env_rest_corr(:,:,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.05);
% find highest correlations (95%)
% prc = prctile(reshape(outp.r_env_task_corr(:,:,1),[121*121 1]),95);

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
imagesc(par,[0 1]); %colormap([nancol; 1 0 0 ])
title('Correlation: p-values (log10)');

[i,j]= find(m2);
idx = round([mean(i) mean(j)]);

for iax = [1 3 4]
 scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
%   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
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
    c.TickLabels = {'0'; sprintf('%.3f',10^-max(c.Limits))};
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

dist      = 0;
iG        = 1; igain = 1;
osc       = osc1(:,:,iG,igain);
oscthresh = 0;

h=figure; set(gcf,'color','w')

% Loop through all subjects
for isubj = 1 : 28
  clear par bw
  
  if dist
    par                 = squeeze(outp.dist_rest_indiv(isubj,:,:,1));
    par(osc>oscthresh)  = nan;  
    par                 = par<prctile(par(:),1);
    bw                  = bwlabel(double(par),8);
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
    
    [i,kk]=find(par==1);
    idx(isubj,:) = [round(mean(i)) round(mean(kk))];
%     [i,j]               = find(par==min(par(:)));
%     idx(isubj,:)        = [i j];
  else
    
    p_corr  = fdr1(reshape(squeeze(outp.p_env_rest_indiv_corr(isubj,:,:,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.05);
    par     = squeeze(outp.p_env_rest_indiv_corr(isubj,:,:,iG,igain)) < p_corr; 
    c       = squeeze(outp.r_env_rest_indiv_corr(isubj,:,:,iG,igain)) > corrwithfc_rest(isubj);
    par     = par&c; 
    
    par(osc>oscthresh)  = 0;
    bw                  = bwlabel(par,8);
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
    r = squeeze(outp.r_env_rest_indiv_corr(isubj,:,:,iG,igain));
    [i,kk]=find(par==1);
    idx(isubj,:) = [round(sum(r(par).*i)/sum(r(par))) round(sum(r(par).*kk)/sum(r(par)))];
    
  end
    % PLOT RESULT
    ax{1}   = subplot(6,5,isubj); hold on
    par = squeeze(-log10(outp.p_env_rest_indiv_corr(isubj,:,:,iG,igain)));
    par(osc>oscthresh)  = nan;
    imagesc(par,[0 3]); colormap(plasma)   
    scatter(gca,idx(isubj,2),idx(isubj,1),20,'markerfacecolor','w','markeredgecolor','k');   
    tp_editplots; axis square
    
    set(h,'Renderer','Painters')
    set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
    set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))

    idx_rest.inh(isubj) = idx(isubj,2);
    idx_rest.exc(isubj) = idx(isubj,1);
end

if ~dist
  print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv.pdf'))
end

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',v_sim),'idx_rest')

clear idx

%% GET INDIV TASK PARAMETERS

clear idx_task idx

indiv_change_prc = 100*(mean(fc_task_indiv,1)-mean(fc_rest_indiv,1))./mean(fc_rest_indiv,1);
 
for isubj = 1:28
  
  if isnan(idx_rest.exc(isubj))
    idx_task.exc(isubj) = nan;
    idx_task.inh(isubj) = nan;
    continue
  end
  
  p_corr  = fdr1(reshape(squeeze(outp.p_env_task_indiv_corr(isubj,:,:,iG,igain)),[size(osc,1)*size(osc,2) 1]),0.05);
  par     = squeeze(outp.p_env_task_indiv_corr(isubj,:,:,iG,igain)) < 0.05; 
  par(osc>oscthresh) = 0;
  
  fc_tvr_sim = 100*(outp.fc_sim_env_mean(:,:,1,1)-outp.fc_sim_env_mean(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1))./outp.fc_sim_env_mean(idx_rest.exc(isubj),idx_rest.inh(isubj),1,1);
  d = indiv_change_prc(isubj)-fc_tvr_sim;
  d(osc>oscthresh) = inf;
  d(par==0) = inf;
  
  % DISTANCE METRIC FROM DEMIRTAS??
  mask = zeros(121,121);
  mask(idx_rest.exc(isubj):end,idx_rest.inh(isubj):end) = 1;
  d(~logical(mask)) = Inf;
  
%   d = abs(d.*(d<0));
%   d(d==0)=Inf;
  m = min(reshape(abs(d),[121*121 1]));
  [i,j]=find(abs(d)==m);
  idx = round([mean(i) mean(j)]);
  idx_task.inh(isubj) = idx(2);
  idx_task.exc(isubj) = idx(1);

end

indiv_idx.rest = idx_rest;
indiv_idx.task = idx_task;
%% PLOT REST AND TASK

figure; set(gcf,'color','white')
load redblue.mat

ax{1} = subplot(2,2,1); hold on
imagesc(zeros(121,121))
set(gca,'ydir','normal')
scatter(gca,idx_rest.inh,idx_rest.exc,25,'markerfacecolor','r','markeredgecolor','k');
axis(ax{1},[1 121 1 121])
set(ax{1},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{1},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
title('Indiv subjects')
axis(ax{1},'square')
tp_editplots; colormap([0.9 0.9 0.9])

ax{2}= subplot(2,2,2); hold on

imagesc(zeros(121,121))
set(gca,'ydir','normal')
scatter(gca,idx_task.inh,idx_task.exc,25,'markerfacecolor','y','markeredgecolor','k');
axis(ax{2},[1 121 1 121])
set(ax{2},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{2},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{2},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
title('Indiv subjects')
axis(ax{2},'square')
tp_editplots

ax{3}= subplot(2,2,3); hold on
imagesc(zeros(121,121))
scatter(gca,round(nanmean(indiv_idx.rest.inh)),round(nanmean(indiv_idx.rest.exc)),35,'markerfacecolor','r','markeredgecolor','k');
scatter(gca,round(nanmean(indiv_idx.task.inh)),round(nanmean(indiv_idx.task.exc)),35,'markerfacecolor','y','markeredgecolor','k');
axis(ax{3},[1 121 1 121])
set(ax{3},'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(ax{3},'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
xlabel(ax{3},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
title('Indiv subjects')
axis(ax{3},'square')
tp_editplots

% set(gca,'ydir','normal')
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_final_fitting_indivfits_taskandrest_v%d.eps',v_sim))

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim),'indiv_idx')

%% LOAD FRACTIION OF ALTERED CORRELATIONS FROM DATA
v = 1;
para.nfreq = 1:13;
para.alpha = 0.05;

if ~exist('emp','var')
  if ~exist('cleandat','var')
    cleandat = pupmod_loadpowcorr(v,1);
  end
  % transform indiv. subj. matrices to AAL BCN
  for isubj = 1 : 28
    isubj
    for im = 1 : 3
      for icont = 1 : 2
        for ifoi = 1:13
          tmp = squeeze(nanmean(cleandat(:,:,isubj,im,icont,ifoi),6));
          tmp = tp_match_aal(para,tmp);
          fc_indiv(:,:,isubj,im,icont,ifoi) = tmp(include_bcn,include_bcn);
        end
      end
    end
  end
  emp = pupmod_compute_altered_correlations(fc_indiv,para);
end

clear fc_indiv cleandat

%% PLOT ALTERED CORRELATIONS AS FUNCTION OF GAIN
clear FC_simrest FC_simtask pp_pos_atx_sim pp_neg_atx_sim

SUBJLIST = find(idx_rest.inh>15);

mask = logical(tril(ones(76,76),-1));
[~,idx_gain]=sort(Gains)

idx_gain1 = idx_gain(1:2:end);

pp_pos_atx_sim = nan(length(idx_gain1),2);
pp_neg_atx_sim = nan(length(idx_gain1),2);

for igain = 1:length(idx_gain(1:2:end))
    
  if Gains(idx_gain1(igain)) == 0
    pp_pos_atx_sim(igain,:) = [nan nan];
    pp_neg_atx_sim(igain,:) = [nan nan];
    continue
  end
  
  gain = idx_gain1(igain)
  
  for isubj = 1:length(SUBJLIST)
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simrest(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simrest(:,:,isubj,2) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simtask(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simtask(:,:,isubj,2) = out.FC_env; clear out
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

% PLOT
% -----------------------
figure; set(gcf,'color','w'); 
subplot(4,2,1); hold on
plot(pp_pos_atx_sim(:,1),'r:');
plot(pp_pos_atx_sim(:,2),'r');
plot(length(idx_gain1)+2,emp.n_p_atx(6,1),'o','markerfacecolor','r')
plot(length(idx_gain1)+2,emp.n_p_atx(6,2),'o','markerfacecolor','r')

axis([0 length(idx_gain1)+2  -0.05 1]); 
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[0 0.25 0.50 0.75 1],'yTickLabels',num2cell([0 0.25 0.50 0.75 1]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

subplot(4,2,3); hold on
plot(pp_neg_atx_sim(:,1),'b:');
plot(pp_neg_atx_sim(:,2),'b');
plot(length(idx_gain1)+2,emp.n_n_atx(6,1),'o','markerfacecolor','b')
plot(length(idx_gain1)+2,emp.n_n_atx(6,2),'o','markerfacecolor','b')
axis([0 length(idx_gain1)+2 -0.05 1]); 
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain1)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[0 0.25 0.50 0.75 1],'yTickLabels',num2cell([0 0.25 0.50 0.75 1]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

% PLOT CONTEXT DEPENDENCE
subplot(4,2,5); hold on
plot(pp_pos_atx_sim(:,1)-pp_pos_atx_sim(:,2),'r');
plot(pp_neg_atx_sim(:,1)-pp_neg_atx_sim(:,2),'b');
plot(length(idx_gain1)+2,emp.n_n_atx(6,1)-emp.n_n_atx(6,2),'o','markerfacecolor','b')
plot(length(idx_gain1)+2,emp.n_p_atx(6,1)-emp.n_p_atx(6,2),'o','markerfacecolor','r')
axis([0 length(idx_gain1)+2  -1 1]); 
set(gca,'XTick',[1:4:length(idx_gain1) length(idx_gain1)+2],'XTickLabels',[num2cell(Gains(idx_gain1(1:4:length(idx_gain1)))) 'Emp.'])
set(gca,'yTick',[-1 -0.5 0 0.5 1],'yTickLabels',num2cell([-1 -0.5 0 0.5 1]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_fraction_as_gain_v%d.pdf',v_sim))
%%
figure; set(gcf,'color','w'); hold on
subplot(2,2,1); hold on
bar([1,2,3], [pp_pos_atx_sim(18,1) pp_pos_atx_sim(18,2) pp_pos_atx_sim(18,1)-pp_pos_atx_sim(18,2)])
axis square; tp_editplots

subplot(2,2,2); hold on
bar([1,2,3], [pp_neg_atx_sim(18,1) pp_neg_atx_sim(18,2) pp_neg_atx_sim(18,1)-pp_neg_atx_sim(18,2)])
axis square; tp_editplots

subplot(2,2,3);
bar([1,2], [pp_pos_atx_sim(18,1)-pp_pos_atx_sim(18,2) pp_neg_atx_sim(18,1)-pp_neg_atx_sim(18,2)])
axis square; tp_editplots

% bar([1,2], )
axis(gca,[0.5 2.5 -0.75 0.75]); tp_editplots; ylabel('Fraction of significantly altered correlation')
tp_editplots

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_fraction_as_gain_barplots_v%d.pdf',v_sim))

%% PLOT ALTERED CORRELATIONS AS FUNCTION OF E-I RATIO
clear p_pos_atx_sim p_neg_atx_sim
SUBJLIST = find(idx_rest.inh>15);
% gain = 18;
mask = logical(tril(ones(76,76),-1));

es = -10:1:8;
is = -10:1:2;


p_pos_atx_sim = nan(length(es),length(is),length(idx_gain1),2);
p_neg_atx_sim = nan(length(es),length(is),length(idx_gain1),2);


for ie = 1:length(es)
  ie
  for ii = 1 : length(is)
    for igain = 1:length(idx_gain1)
      
      if es(ie) == 0 && is(ii) == 0 && Gains(idx_gain1(igain)) == 0
        p_pos_atx_sim(ie,ii,igain,:) = [nan nan];
        p_neg_atx_sim(ie,ii,igain,:) = [nan nan];
        continue
      end
      
      gain = idx_gain1(igain);
      
      exc = es(ie);
      inh = is(ii);
      
      for isubj = SUBJLIST
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),1,1,v_sim))
        FC_simrest(:,:,isubj,1) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(isubj)+exc,indiv_idx.rest.inh(isubj)+inh,1,gain,v_sim))
        FC_simrest(:,:,isubj,2) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),1,1,v_sim))
        FC_simtask(:,:,isubj,1) = out.FC_env; clear out
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(isubj)+exc,indiv_idx.task.inh(isubj)+inh,1,gain,v_sim))
        FC_simtask(:,:,isubj,2) = out.FC_env; clear out
      end
      
      [h,~,~,s]=ttest(FC_simrest(:,:,:,2),FC_simrest(:,:,:,1),'dim',3);
      p_pos_atx_sim(ie,ii,igain,1) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
      p_neg_atx_sim(ie,ii,igain,1) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
      
      [h,~,~,s]=ttest(FC_simtask(:,:,:,2),FC_simtask(:,:,:,1),'dim',3);
      p_pos_atx_sim(ie,ii,igain,2) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
      p_neg_atx_sim(ie,ii,igain,2) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
      
%       [h,~,~,s]=ttest(FC_simtask(:,:,:,2),FC_simrest(:,:,:,2),'dim',3);
%       p_pos_tvr_sim(ie,ii,igain) = sum(s.tstat(mask)>0 & h(mask)) / sum(mask(:));
%       p_neg_tvr_sim(ie,ii,igain) = sum(s.tstat(mask)<0 & h(mask)) / sum(mask(:));
%       
%       diff_tvr(ie,ii,igain) = sum(sum(triu(nanmean(FC_simtask(:,:,:,2),3)-nanmean(FC_simrest(:,:,:,2),3),1)))./sum(mask(:));
      fc_sim_rest(ie,ii,igain,:) = mean(squeeze(mean(mean(FC_simrest),2)));
      fc_sim_task(ie,ii,igain,:) = mean(squeeze(mean(mean(FC_simtask),2)));
      
    end
  end
end

%%
% PLOT
% -----------------------
figure; set(gcf,'color','w'); 
subplot(2,2,1); hold on
imagesc(p_pos_atx_sim(:,:,1),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
subplot(2,2,2); hold on
imagesc(p_neg_atx_sim(:,:,1),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
subplot(2,2,3); hold on
imagesc(p_pos_atx_sim(:,:,2),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
subplot(2,2,4); hold on
imagesc(p_neg_atx_sim(:,:,2),[0 1]); set(gca,'ydir','normal'); axis tight square; tp_editplots
% 
colormap(plasma)


figure; set(gcf,'color','w'); 
subplot(2,3,1); hold on
imagescnan(p_pos_atx_sim(:,:,Gains(idx_gain1)==0,1) - p_pos_atx_sim(:,:,Gains(idx_gain1)==0,2),[-0.85 0.85])
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))

xlabel('Background input (to I)'); ylabel('Background input (to E)')
subplot(2,3,2); hold on
imagescnan(p_pos_atx_sim(:,:,18,1) - p_pos_atx_sim(:,:,18,2),[-0.85 0.85])
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
xlabel('Background input (to I)'); ylabel('Background input (to E)')
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))

% subplot(2,2,2); hold on
a=p_pos_atx_sim(:,:,:,1)<0.1;
b=p_pos_atx_sim(:,:,:,2)>0.6;
c=p_neg_atx_sim(:,:,:,1)<0.1;
d=p_neg_atx_sim(:,:,:,2)<0.1;


subplot(2,3,3); hold on
par = a&b&c&d;
idf_gains = find(any(any(par==1)));

imagescnan(par(:,:,idf_gains(end))); 
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
xlabel('Background input (to I)'); ylabel('Background input (to E)')
title(sprintf('Gain = %.3f',Gains(idx_gain1(idf_gains(end)))))
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))

subplot(2,3,4); hold on
imagescnan(p_neg_atx_sim(:,:,Gains(idx_gain1)==0,1) - p_neg_atx_sim(:,:,Gains(idx_gain1)==0,2),[-0.4 0.4])
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
xlabel('Background input (to I)'); ylabel('Background input (to E)')
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))

subplot(2,3,5); hold on
imagescnan(p_neg_atx_sim(:,:,14,1) - p_neg_atx_sim(:,:,14,2),[-0.4 0.4])
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
xlabel('Background input (to I)'); ylabel('Background input (to E)')
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))


a=p_pos_atx_sim(:,:,:,1)<0.15;
b=p_pos_atx_sim(:,:,:,2)<0.05;
c=p_neg_atx_sim(:,:,:,1)>0.30;
d=p_neg_atx_sim(:,:,:,2)<0.075;
par = a&b&c&d;
idf_gains = find(any(any(par==1)));

subplot(2,3,6); hold on
imagesc(par(:,:,idf_gains(1))); 
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(parula)
xlabel('Background input (to I)'); ylabel('Background input (to E)')
title(sprintf('Gain = %.3f',Gains(idx_gain1(14))))
idx = find(sum(a&b&c&d,2)==max(sum(a&b&c&d,2)));
set(gca,'xtick',[1 5 9 13],'xticklabel',is([1 5 9 13])*mean(diff(Ies)))
set(gca,'ytick',[1 5 9 13 17],'yticklabel',es([1 5 9 13 17])*mean(diff(Ies)))

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_eiphaseplot_v%d.pdf',v_sim))


figure; set(gcf,'color','w'); hold on
subplot(2,2,1); hold on
bar([1,2,3], [pp_pos_atx_sim(18,1) pp_pos_atx_sim(18,2) pp_pos_atx_sim(18,1)-pp_pos_atx_sim(18,2)])
axis([0.5 2.5 -0 0.8]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots; axis square
subplot(2,2,2); hold on
bar([1,2,3], [pp_neg_atx_sim(18,1) pp_neg_atx_sim(18,2) pp_neg_atx_sim(18,1)-pp_neg_atx_sim(18,2)])
axis([0.5 2.5 0 0.8]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots; axis square

subplot(2,2,3); hold on
bar([1,2], [p_pos_atx_sim(2,10,14,1)-p_pos_atx_sim(2,10,14,2) p_neg_atx_sim(2,10,14,1)-p_neg_atx_sim(2,10,14,2)])
axis([0.5 2.5 -0.4 0.4]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots; axis square
% subplot(2,2,4); hold on
% bar([1,2,3], [])
% axis([0.5 2.5 0 0.4]); tp_editplots; ylabel('Change in correlation [in %]')
% tp_editplots; axis square

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gainmodulation_barplot_v%d.pdf',v_sim))

figure; set(gcf,'color','w'); 

% PLOT PERCENT CHANGE IN CORRELATIONS
subplot(2,2,1); hold on
bar([1,2], [100*(fc_sim_rest(2,10,14,2)-fc_sim_rest(2,10,14,1))./fc_sim_rest(2,10,14,1) (fc_sim_task(2,10,14,2)-fc_sim_task(2,10,14,1))./fc_sim_task(2,10,14,1)])
axis([0.5 2.5 -75 75]); tp_editplots; ylabel('Change in correlation [in %]')
tp_editplots; axis square

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_gain_pct_change_dpz_barplot_v%d.pdf',v_sim))


%% LOAD PARAMETERS AND PLOT GAIN VS NO GAIN
% For specific level of gain
% ----------
clear p_* FC_simrest FC_simtask p_pos_atx_sim p_neg_atx_sim permdat_atx_rest permdat_atx_task perm
isperm = 1;
SUBJLIST

mask = logical(tril(ones(76,76),-1));
gain = 45;
% 
if isperm
  load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim))

  for isubj = 1:length(SUBJLIST)
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simrest(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(SUBJLIST(isubj)),indiv_idx.rest.inh(SUBJLIST(isubj)),1,gain,v_sim))
    FC_simrest(:,:,isubj,2) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,1,v_sim))
    FC_simtask(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(SUBJLIST(isubj)),indiv_idx.task.inh(SUBJLIST(isubj)),1,gain,v_sim))
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

  nperm = 1000;
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
imagesc(par,[-0.02 0.02])
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



%% PLOT PERCENT CHANGE AS FUNCTION OF GAIN
% ------------------------------

mask = logical(tril(ones(76,76),-1));
[~,idx_gain1]=sort(Gains)

for igain = 1:length(idx_gain1)
  
  if Gains(idx_gain1(igain)) == 0
    p_pos_atx_sim(igain,:) = [nan nan];
    p_neg_atx_sim(igain,:) = [nan nan];
    continue
  end
  
  gain = idx_gain1(igain)
  
  for isubj = SUBJLIST
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),1,1,v_sim))
    FC_simrest(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.rest.exc(isubj),indiv_idx.rest.inh(isubj),1,gain,v_sim))
    FC_simrest(:,:,isubj,2) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),1,1,v_sim))
    FC_simtask(:,:,isubj,1) = out.FC_env; clear out
    load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',indiv_idx.task.exc(isubj),indiv_idx.task.inh(isubj),1,gain,v_sim))
    FC_simtask(:,:,isubj,2) = out.FC_env; clear out
  end
  
  prc_change_rest(igain) = 100*(mean(mean(mean(FC_simrest(:,:,:,2))))-mean(mean(mean(FC_simrest(:,:,:,1)))))/mean(mean(mean(FC_simrest(:,:,:,1))));
  prc_change_task(igain) = 100*(mean(mean(mean(FC_simtask(:,:,:,2))))-mean(mean(mean(FC_simtask(:,:,:,1)))))/mean(mean(mean(FC_simtask(:,:,:,1))));
  
  
  
%   diff_tvr(igain) = sum(sum(triu(nanmean(FC_simtask(:,:,:,2),3)-nanmean(FC_simrest(:,:,:,2),3),1)))./sum(mask(:));

end
%%
figure; set(gcf,'color','w')
subplot(1,2,1); hold on
bar(prc_change_rest(18)-prc_change_task(18))
axis square; tp_editplots
axis([0.5 1.5 -400 400])

print(gcf,'-depsc2',sprintf('~/pmod/plots/pupmod_rawcorr_bar_v%d.eps',v))
%%
% PLOT
% -----------------------
figure; set(gcf,'color','w'); 
subplot(4,2,1); hold on
plot(prc_change_rest,'r:');
plot(prc_change_task,'r');
% plot(length(idx_gain)+2,emp.n_p_atx(6,1),'o','markerfacecolor','r')
% plot(length(idx_gain)+2,emp.n_p_atx(6,2),'o','markerfacecolor','r')
line([17 17],[0 0.75],'color','k')
line([find(Gains(idx_gain)==Gains(12)) find(Gains(idx_gain)==Gains(12))],[0 0.5],'color','k')

axis([0 length(idx_gain)+2  -100 800]); 
set(gca,'XTick',[1:4:length(idx_gain) length(idx_gain)+2],'XTickLabels',[num2cell(Gains(idx_gain(1:4:length(idx_gain)))) 'Emp.'])
set(gca,'yTick',[0:200:800],'yTickLabels',num2cell([0:200:800]))
xlabel('Change in gain'); %ylabel('Fraction of altered correlations');
tp_editplots
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

if ~exist('emp','var')
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
end
%%
clear cleandat

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

subplot(2,2,3); hold on
imagesc(squeeze(r(1,2,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))

subplot(2,2,2); hold on
imagesc(squeeze(r(2,1,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))

subplot(2,2,4); hold on
imagesc(squeeze(r(2,2,:,:))',[0.9 1.1]);
set(gca,'ydir','normal'); axis tight square; tp_editplots; colormap(plasma)
set(gca,'xTick',1:4:45,'xTickLabels',num2cell([Gains(idx_gain(1:4:end))]))

for i = 1 : 42
  for isubj = 1 : 28
    kk = i:i+2;
    smoothed_task(i,isubj,1) = mean((r(2,1,kk,isubj)));
    smoothed_task(i,isubj,2) = mean((r(2,2,kk,isubj)));
  end
end
  
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
  
    