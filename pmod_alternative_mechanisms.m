
% loop through resting points

igain = 6;

CeE = squeeze(outp.fc_sim_fr_mean);

par=CeE(:,:,igain);
par(osc1(:,:,:,igain)>0)=nan;
CeE(:,:,igain)=par;

for ie = 1:121
  ie
  for ii = 1 : 121
    if isnan(CeE(ie,ii,igain))
      continue
    end
    
    for task_effect = 1 : 120
      
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
          Drg_Rest = CeE(ie+drug_eff_on_E,ii+drug_eff_on_I,igain);
          Drg_Task = CeE(ie+drug_eff_on_E+task_effect,ii+drug_eff_on_I+task_effect,igain);
          
          drg_eff_rest = Drg_Rest-Pbo_Rest;
          drg_eff_task = Drg_Task-Pbo_Task;
          
          
          atx_eff = drg_eff_task>0.005 & abs(drg_eff_rest)<0.001;
          
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
                
                drg_eff_rest = Drg_Rest-Pbo_Rest;
                drg_eff_task = Drg_Task-Pbo_Task;
                
                dpz_eff = abs(drg_eff_task)<0.001 & drg_eff_rest<-0.005;
                
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


%%
figure; set(gcf,'color','w')
par=CeE(:,:,6);
c=imagesc(par,[0 0.02]); colormap(plasma); hold on

set(c,'AlphaData',~isnan(par))

scatter(ii,ie,'o','markeredgecolor','w')
scatter(ii+task_effect,ie+task_effect,'o','markeredgecolor','w')
scatter(ii+atx_eff_on_I,ie+atx_eff_on_E,'o','markeredgecolor','k')
scatter(ii+dpz_eff_on_I,ie+dpz_eff_on_E,'o','markeredgecolor','r')
scatter(ii+task_effect+dpz_eff_on_I,ie+task_effect+dpz_eff_on_E,'o','markeredgecolor','r')
scatter(ii+task_effect+atx_eff_on_I,ie+task_effect+atx_eff_on_E,'o','markeredgecolor','k')

set(gca,'ydir','normal')
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
tp_editplots
axis square
xlabel('Background input to I')
ylabel('Background input to E')