%% FITTING
% pmod_final_fitting

clear
% %-------------------------------------------------------------------------
% VERSION 1:
% %-------------------------------------------------------------------------
v           = 1;
Ies         = -4:0.05:-1;
Iis         = -5:0.05:-2;
Gg          = 0:0.05:2;
Gains       = 0; 
nTrials     = 1;
tmax        = 6500;  % in units of tauE
dt          = 0.01;
%-------------------------------------------------------------------------

v_sim = v;
% connectivity, AAL
v_conn =  33;

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
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
fc_rest = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,1,1:end-1),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,1,2,1:end-1),3),6));

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

fc_rest = tp_match_aal(para,fc_rest);
fc_task = tp_match_aal(para,fc_task);

fc_rest = fc_rest(include_bcn,include_bcn);
fc_task = fc_task(include_bcn,include_bcn);

% transform indiv. subj. matrices to AAL BCN
for isubj = 1 : size(cleandat,3)
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,1,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  [corrwithfc_rest(isubj), p_corrwithfc_rest(isubj)]  = corr(tmp(mask),SC(mask));
  fc_rest_indiv(:,isubj) = tmp(mask);
  tmp = squeeze(nanmean(cleandat(:,:,isubj,1,2,[3 4 5 6]),6));
  tmp = tp_match_aal(para,tmp); tmp = tmp(include_bcn,include_bcn);
  fc_task_indiv(:,isubj) = tmp(mask);
  [corrwithfc_task(isubj), p_corrwithfc_task(isubj)] = corr(tmp(mask),SC(mask));
end

k = [0 0];
%% CORRELATION FC and SC

figure; set(gcf,'color','w');

subplot(3,4,1); hold on

r = (rand(28,1)-0.5)./5;

plot(r+ones(length(corrwithfc_rest),1),corrwithfc_rest,'o','markersize',5,'markeredgecolor','w','MarkerFaceColor','k','linewidth',0.5);
set(gca,'XTick',1:3:28,'XTickLabel',num2cell([1:3:28]))
line([0.75 1.25],[mean(corrwithfc_rest) mean(corrwithfc_rest)],'linestyle','-','color','k')

plot(r+2*ones(length(corrwithfc_task),1),corrwithfc_task,'o','markersize',5,'markeredgecolor','w','MarkerFaceColor',[0.5 0.6 0.6],'linewidth',0.5);
line([1.75 2.25],[mean(corrwithfc_rest) mean(corrwithfc_rest)],'linestyle','-','color',[0.5 0.6 0.6])
axis([0.5 2.5 0 0.4]); tp_editplots;
ylabel('Correlation FC w SC');

set(gca,'YTick',0:0.1:0.4,'YTickLabel',num2cell([0:0.1:0.4]))
tp_editplots

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_corrwithSC_v%d.pdf',v))

%%
if ~exist(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_globalcoupling_v%d.mat',v_sim,v_sim))
  
  for iies = 1 : length(Ies)
    iies
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        for iiis = 1 :length(Iis)
          
          load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v_sim,iies,iiis,iG,igain,v_sim))
          
          outp.fc_sim_fr_tmp = out.fc_FR;
          
          [outp.r_fr_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest(mask));
          [outp.r_fr_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_task(mask));          
          
          [outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_rest_indiv);
          [outp.r_fr_task_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_task_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_fr_tmp(mask),fc_task_indiv);
          
          outp.dist_fr_rest(iies,iiis,iG,igain) = 1-(outp.r_fr_rest_corr(iies,iiis,iG,igain)-(mean(fc_rest(mask))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          outp.dist_fr_task(iies,iiis,iG,igain) = 1-(outp.r_fr_task_corr(iies,iiis,iG,igain)-(mean(fc_task(mask))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          
          outp.dist_fr_rest_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_fr_rest_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_rest_indiv))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          outp.dist_fr_task_indiv(:,iies,iiis,iG,igain) = 1-(squeeze(outp.r_fr_task_indiv_corr(:,iies,iiis,iG,igain))'-(squeeze(mean(fc_task_indiv))-mean(outp.fc_sim_fr_tmp(mask))).^2);
          
          outp.fc_sim_fr_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_fr_tmp(mask));
          
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          
        end
      end
    end
  end
%   
  save(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_globalcoupling_v%d.mat',v_sim,v_sim),'outp')
% 
else
  load(sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_globalcoupling_v%d.mat',v_sim,v_sim))
  load(sprintf('~/pmod/proc/detosc/v%d/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim,v_sim))
end
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
%% PLOT SIMILARITY FC_EMP and FC_SIM AS A FUNCTION OF GLOBAL COUPLING

idxE = find(Ies == -1);
idxI = find(Iis == -2);

osc = osc1(1:idxE,1:idxI,:);
% osc = osc1;

for isubj = 1 : size(cleandat,3)
  for igg = 1 : 41
    
    par = 1-(squeeze(outp.r_fr_rest_indiv_corr(isubj,1:idxE,1:idxI,igg))'-(squeeze(mean(fc_rest_indiv(:,isubj)))-mean(outp.fc_sim_fr_tmp(mask))).^2);
    par(osc(:,:,igg)==1)=nan;
    d(isubj,igg) = nanmean(par(:));
    
    tmp = squeeze(outp.r_fr_rest_indiv_corr(isubj,:,:,igg));
    m(isubj,igg)=mean(tmp(osc(:,:,igg)>0));
    
    mmaaxx(isubj,igg)=max(tmp(osc(:,:,igg)>0));
    
  end
end

figure; set(gcf,'color','w');
subplot(3,2,1);
imagesc(m,[0 0.5]); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
% set(gca,'YTick',1:28,'YTickLabel',num2cell(1:28));
ylabel('Subjects');
colorbar; colormap(plasma)
subplot(3,2,2);
plot(mean(m)); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
[~,i]=max(mean(m));
axis([0 42 -0.05 0.5])
line([i i],[-0.05 0.5])

set(gca,'YTick',0:0.1:0.4,'YTickLabel',num2cell(0:0.1:0.4));
ylabel('Mean correlation')

subplot(3,2,3);
imagesc(d,[0.6 1.0]); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
% set(gca,'YTick',1:28,'YTickLabel',num2cell(1:28));
ylabel('Subjects');
colorbar; colormap(plasma)
subplot(3,2,4);
plot(mean(d)); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
axis([0 42 0.6 1.1])
[~,i]=min(mean(d));
set(gca,'YTick',0.6:0.1:1.1,'YTickLabel',num2cell(0.6:0.1:1.1));
ylabel('Mean distance <\delta>')
line([i i],[.7 1.1])

subplot(3,2,5);
imagesc(mmaaxx,[0 0.7]); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
% set(gca,'YTick',1:28,'YTickLabel',num2cell(1:28));
ylabel('Subjects');
colorbar; colormap(plasma)
[~,i]=min(mean(d));


subplot(3,2,6);
plot(mean(mmaaxx)); axis square; tp_editplots
set(gca,'XTick',1:10:length(Gg),'XTickLabel',num2cell(Gg(1:10:end)))
xlabel('Global coupling')
axis([0 42 0 0.8])
set(gca,'YTick',0:0.2:0.8,'YTickLabel',num2cell(0:0.2:0.8));
ylabel('Max. correlation')


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_globalcoupling_indivsubj_v%d.pdf',v_sim))
