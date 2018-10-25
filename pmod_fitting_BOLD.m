%% FITTING
% pmod_wc_wholebrain_final_plot.m
clear
% %-------------------------------------------------------------------------
% VERSION 2: After meeting with tobi, 24-08-2018: even more fine gained
% %-------------------------------------------------------------------------
v           = 1;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0.6;
Gains       = 0;
nTrials     = 1;
tmax        = 6500; % in units of tauE
wins = [2 20];
%-------------------------------------------------------------------------
% VERSION 1: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v           = 4;
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-1;
% Gg          = 0:0.2:1;
% Gains       = 0:0.05:0.2;
% nTrials     = 3;
% tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------

v_sim = v;
% connectivity, AAL
v_conn =  1;
% simulations
% dfa, aal, lcmv
v_dfa = 2;

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
Gs = 1;

% ---------
% LOAD EMPIRICAL DFA (LCMV)
% ---------
load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],v_dfa));
% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
load ~/M.mat

pars.dim = 2;
pars.grid = 'medium';
pars.N = 90;
pars.transfer = 'to_bcn';

% dfa_emp_rest = nanmean(outp.dfa_all(:,:,1,1,1),2);
% dfa_emp_task = nanmean(outp.dfa_all(:,:,1,1,2),2);

lambda_emp_rest = nanmean(outp.lambda_all(:,:,1,1,1),2);
lambda_emp_task = nanmean(outp.lambda_all(:,:,1,1,2),2);

clear outp

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v_conn));
load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat
peakfreq_rest = m_res(1);
peakfreq_task = m_tsk(1);

load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

mask = logical(tril(ones(90,90),-1));
mask = find(triu(ones(90))-eye(90));

% transform avg fc matrices to AAL BCN
fc_rest     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,1,6),3)));
fc_task     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,2,6),3)));

% transform indiv. subj. matrices to AAL BCN
for isubj = 1 : 28
  fc_rest_indiv(:,isubj) = reshape(tp_match_aal(pars,squeeze(cleandat(:,:,isubj,1,1,6))),[90*90 1]);
  fc_task_indiv(:,isubj) = reshape(tp_match_aal(pars,squeeze(cleandat(:,:,isubj,1,2,6))),[90*90 1]);
end
% vectorize
fc_rest_indiv = fc_rest_indiv(mask,:);
fc_task_indiv = fc_task_indiv(mask,:);

% transform BOLD FC to BCN
for isubj = 1 : 24
  fc_BOLD_emp(:,isubj) = reshape(tp_match_aal(pars,corrcoef(squeeze(M(isubj,3,:,:)))),[90*90 1]);
end
fc_BOLD_emp_indiv = fc_BOLD_emp(mask,:);

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
  
  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          %         igain
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))
          
          if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
            disp('save stuff')
            FFCC = out.FC;
          end
          
          % Time scales
          outp.lambda(:,iies,iiis,iG,igain)      = mean(out.lambda,2);
          outp.lambda_env(:,iies,iiis,iG,igain)  = mean(out.lambda_env,2);
          % DFA
          %         outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
          %         outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
          
          %         [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          
          % get peak frequency from MEG data
          %         outp.peakfreq_diff_res(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_rest;
          %         outp.peakfreq_diff_tsk(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_task;
          
          pars.dim = 2;
          
          outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);
          outp.fc_sim_env_tmp = tp_match_aal(pars,out.FC_env,3);
          
          %         outp.fc_sim_var(:,iies,iiis,iG,igain) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
          
          [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
          [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
          [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
          
          [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
          [outp.r_rest_BOLD_corr(iies,iiis,iG,igain), outp.p_rest_BOLD_corr(iies,iiis,iG,igain)]=corr(outp.fc_BOLD(mask),fc_BOLD_emp_indiv);
          [outp.r_rest_BOLD_corr(iies,iiis,iG,igain), outp.p_rest_BOLD_corr(iies,iiis,iG,igain)]=corr(outp.fc_BOLD(mask),fc_BOLD(mask));
          
          %
          pars.dim = 1;
          
          %         [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
          %        	[outp.dfa_env_r_rest(iies,iiis,iG,igain), outp.dfa_env_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_env_sim(:,iies,iiis,iG,igain),[1 90]),pars));
          
          [outp.lambda_r_rest(iies,iiis,iG,igain), outp.lambda_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
          [outp.lambda_r_task(iies,iiis,iG,igain), outp.lambda_p_task(iies,iiis,iG,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
          
          [outp.lambda_env_r_rest(iies,iiis,iG,igain), outp.lambda_env_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,iies,iiis,iG,igain),[1 90]),pars));
          
          %         outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
          %         outp.dist_fc_task (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
          
          outp.fc_sim_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_tmp(mask));
          %         outp.fc_sim_all(:,:,iies,iiis,iG,igain) = outp.fc_sim_tmp;
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          %
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          
          % KURAMOTO
          outp.kuramoto_mean (iies,iiis,iG,igain) = mean(out.KOPmean);
          outp.kuramoto_std (iies,iiis,iG,igain)  = mean(out.KOPsd);
          
          %         outp.kuramoto_mean_diff (iies,iiis,iG,igain) = mean(out.KOPmean) - mean(kura_mean(:,1,1,1));
          %         outp.kuramoto_std_diff (iies,iiis,iG,igain)  = mean(out.KOPsd)- mean(kura_std(:,1,1,1));
          
          %         outp.psslp(:,iies,iiis,iG,igain)      = out.psslope;
          %         outp.psslp_env(:,iies,iiis,iG,igain)  = out.psslope_env;
          
          fclose all;
          %
        end
      end
    end
  end
  
  
  %
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim),'outp')
  % DELETE OUTOUT FILES
  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      %     iiis
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          
          fn =  sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v_sim);
          delete(sprintf('%s_processing.txt',fn))
          delete(sprintf('%s.mat',fn))
          
        end
      end
    end
  end
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
end

error('!')

%% INDIV SUBJ
iG = 1;

% osc=zeros(301,401)
oscthresh = 0
h=figure; set(gcf,'color','w')

width = 4;
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
clear idx
ax{1}=subplot(6,5,1); hold on

imagescnan(bif_mask,[-1 1])

colormap(redblue)
set(gca,'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(gca,'YTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots

for isubj = 1 : 24
  % plot lambda
  ax{1} = subplot(6,5,isubj+1); hold on
  p_corr = fdr1(reshape(squeeze(outp.p_rest_BOLD_indiv_corr(isubj,:,:,iG,1)),[size(osc,1)*size(osc,2) 1]),0.05);
  par = squeeze(outp.p_rest_BOLD_indiv_corr(isubj,:,:,iG,1))<p_corr; par=double(par);

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

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_BOLD_v%d.mat',v_sim),'idx')

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



ax{1} = subplot(2,2,4); hold on
scatter(gca,idx(:,2),idx(:,1),10,'markerfacecolor','w','markeredgecolor','k');
axis(gca,[1 size(osc1,1) 1 size(osc1,2)])

set(ax{1},'YTick',1:20:length(Ies ),'YTickLabels',num2cell(Ies(1:20:end)))
set(ax{1},'XTick',1:20:length(Iis),'XTickLabels',num2cell(Iis(1:20:end)))
tp_editplots
xlabel(ax{1},'Inhibitory input')
ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
title('Indiv subjects')
colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_sum.pdf'))


%% GET INDIV TASK PARAMETERS
indiv_change_prc = 100*(mean(fc_task_indiv,1)-mean(fc_rest_indiv,1))./mean(fc_rest_indiv,1);




for isubj = 1:28
  if isnan(idx(isubj,1)) || idx(isubj,1)==0
    continue
  end
  fc_tvr_sim = 100*(outp.fc_sim_mean(:,:,4,1)-outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),4,1))./outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),4,1)

  d = indiv_change_prc(isubj)-fc_tvr_sim;
  d(osc>0) = nan;

  mask = zeros(31,41);
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
  
  mask = mask&bif_mask

thresh = 1;
% d(~logical(bif_mask)) = Inf;
d(~logical(mask)) = Inf;
% d(1:idx(isubj,1),1:idx(isubj,2)) = Inf;

while 1
  clust = reshape(bwlabel(abs(d)<thresh,8),[31*41 1]);
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

imagesc(zeros(31,41))
scatter(ax{1},idx_task(:,2),idx_task(:,1),20,'markerfacecolor','y','markeredgecolor','k');
set(ax{1},'YTick',1:10:length(Ies),'YTickLabels',num2cell(Ies(1:10:end)))
set(ax{1},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
axis([1 41 1 31])
tp_editplots; colormap(redblue); colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_task_sum.pdf'))

clear par

indiv_idx.rest = idx;
indiv_idx.task = idx_task;

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim),'indiv_idx')

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
