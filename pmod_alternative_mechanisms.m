% pmod_alternative_mechanisms.m

v_sim = 3;

fc = pupmod_loadpowcorr(v_conn,SUBJLIST,1);

% transform avg fc matrices to AAL BCN
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

% ifoi = 1;
fc_rest = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,[1 2 3],1,[3 4 5 6]),3),6));
fc_task = squeeze(nanmean(nanmean(cleandat(1:90,1:90,:,[1 2 3],2,[3 4 5 6]),3),6));

para          = [];
para.transfer = 'to_bcn';
para.N        = 90;

for im = 1 : 3
fc_rest(:,:,im) = tp_match_aal(para,fc_rest(:,:,im));
fc_task(:,:,im) = tp_match_aal(para,fc_task(:,:,im));
end
fc_rest=fc_rest(include_bcn,include_bcn,:);
fc_task=fc_task(include_bcn,include_bcn,:);

diff_atx_rest = nanmean(nanmean(fc_rest(:,:,2),1),2)-nanmean(nanmean(fc_rest(:,:,1),1),2);
diff_atx_task = nanmean(nanmean(fc_task(:,:,2),1),2)-nanmean(nanmean(fc_task(:,:,1),1),2);
diff_dpz_rest = nanmean(nanmean(fc_rest(:,:,3),1),2)-nanmean(nanmean(fc_rest(:,:,1),1),2);
diff_dpz_task = nanmean(nanmean(fc_task(:,:,3),1),2)-nanmean(nanmean(fc_task(:,:,1),1),2);


if ~exist('outp','var')
  error('Load outp from final fitting first!')
end

% load idx_rest
load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',v_sim))

% loop through resting points
%%
for igain = 6

par=outp.fc_sim_fr_mean(:,:,:,igain);
par(osc1(:,:,:,igain)>0)=nan;
CeE(:,:,igain)=par;

for ie = round(mean(idx_rest.exc))
  ie
  for ii = round(mean(idx_rest.inh))
    if isnan(CeE(ie,ii,igain))
      continue
    end
    
    for task_effect = 1 : 120
%       task_effect
      if (ie+task_effect) > 121
        continue
      elseif (ii+task_effect) > 121
        continue
      end
      
      eff_E = -(ie-1):(121-(ie+task_effect));
      eff_I = -(ii-1):(121-(ii+task_effect));
      
      for i_drug_eff_on_E = 1:length(eff_E)
        for i_drug_eff_on_I = 1:length(eff_I)
          
          drug_eff_on_E = eff_E(i_drug_eff_on_E);
          drug_eff_on_I = eff_I(i_drug_eff_on_I);
          
          if ie+drug_eff_on_E<=1
            continue
          elseif ii+drug_eff_on_I<=1
            continue
          elseif ie+task_effect+drug_eff_on_E>121
            continue
          elseif ii+task_effect+drug_eff_on_I>121
            continue
          end
          
          Pbo_Rest = CeE(ie,ii,1);
          Pbo_Task = CeE(ie+task_effect,ii+task_effect,1);
          
          if ~(Pbo_Task-Pbo_Rest<-0.003)
            continue
          end
          
          Drg_Rest = CeE(ie+drug_eff_on_E,ii+drug_eff_on_I,igain);
          Drg_Task = CeE(ie+drug_eff_on_E+task_effect,ii+drug_eff_on_I+task_effect,igain);
          
          drg_eff_rest = Drg_Rest-Pbo_Rest;
          drg_eff_task = Drg_Task-Pbo_Task;
          
          corr_val(1,1) = Pbo_Rest;
          corr_val(1,2) = Pbo_Task;
          corr_val(2,1) = Drg_Rest;
          corr_val(2,2) = Drg_Task;
          
          
          
          atx_eff = drg_eff_task>0.007 & abs(drg_eff_rest)<0.001;
          
          if ~atx_eff
            continue
          else
            
            atx_eff_on_E = drug_eff_on_E;
            atx_eff_on_I = drug_eff_on_I;
 
            %            fprintf('Atx effect found...\n')
            
            for ii_drug_eff_on_E = 1:length(eff_E)
              for ii_drug_eff_on_I = 1:length(eff_I)
                
                drug_eff_on_E = eff_E(ii_drug_eff_on_E);
                drug_eff_on_I = eff_I(ii_drug_eff_on_I);
                
                if ie+drug_eff_on_E<=1
                  continue
                elseif ii+drug_eff_on_I<=1
                  continue
                elseif ie+task_effect+drug_eff_on_E>121
                  continue
                elseif ii+task_effect+drug_eff_on_I>121
                  continue
                end
                
                Pbo_Rest = CeE(ie,ii,1);
                Pbo_Task = CeE(ie+task_effect,ii+task_effect,1);
                Drg_Rest = CeE(ie+drug_eff_on_E,ii+drug_eff_on_I,igain);
                Drg_Task = CeE(ie+drug_eff_on_E+task_effect,ii+drug_eff_on_I+task_effect,igain);
                
                corr_val(3,1) = Drg_Rest;
                corr_val(3,2) = Drg_Task;
          
                drg_eff_rest = Drg_Rest-Pbo_Rest;
                drg_eff_task = Drg_Task-Pbo_Task;
                
                dpz_eff = abs(drg_eff_task)<0.003 & drg_eff_rest<-0.007;
                
                if atx_eff & dpz_eff
                  dpz_eff_on_E = drug_eff_on_E;
                  dpz_eff_on_I = drug_eff_on_I;
                  fprintf('Atx & Dpz effect found...\n')
                  error('Found!')
                end
                
              end
            end
          end 
        end
      end
    end
  end
end
end

%%
figure; set(gcf,'color','w')
subplot(2,3,[1 2 4 5])
par=CeE(:,:,igain);
c=imagesc(par,[0 0.02]); colormap(plasma); hold on

set(c,'AlphaData',~isnan(par))

scatter(ii,ie,'o','markeredgecolor','w','markerfacecolor','k')
scatter(ii+task_effect,ie+task_effect,'o','markeredgecolor','w','markerfacecolor',[.4 .4 .4])
scatter(ii+atx_eff_on_I,ie+atx_eff_on_E,'o','markeredgecolor','w','markerfacecolor',[1 .2 .0])
scatter(ii+dpz_eff_on_I,ie+dpz_eff_on_E,'o','markeredgecolor','w','markerfacecolor',[0 .2 1])
scatter(ii+task_effect+dpz_eff_on_I,ie+task_effect+dpz_eff_on_E,'o','markeredgecolor','w','markerfacecolor',[.2 .5 1])
scatter(ii+task_effect+atx_eff_on_I,ie+task_effect+atx_eff_on_E,'o','markeredgecolor','w','markerfacecolor',[1 .5 .2])

set(gca,'ydir','normal')
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
tp_editplots
axis square
% xlabel('Background input to I')
% ylabel('Background input to E')


subplot(2,3,3)

bar([1 2 3 4],[corr_val(1,1) corr_val(2,1) corr_val(1,2) corr_val(2,2)]); 
axis square; tp_editplots
axis([0 5 0 0.05])

subplot(2,3,6)

bar([1 2 3 4],[corr_val(1,1) corr_val(3,1) corr_val(1,2) corr_val(3,2)]); 
axis square; tp_editplots   
axis([0 5 0 0.05])

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_alternative_mech_v%d.pdf',v_sim))
