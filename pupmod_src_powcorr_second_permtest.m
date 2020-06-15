%% pupmod_src_powcorr_second_permtest
% counts number of altered correlations and computes permututation test
% see hawellek et al. (2013) for details
% cleaned data is generated in pupmod_all_powcorr_periphereal

% data is used in: pupmod_plot_alteredcorr

% last update: 26-10-2018

clear 

v     = 3; 
nperm = 10000; 
alp   = 0.05;
nfoi  = 17;

par.subs     = 200;
par.allperms = nperm/par.subs;

outdir = '/home/tpfeffer/pupmod/proc/conn/';
addpath /home/tpfeffer/pconn/matlab/

SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

addpath ~/pconn/matlab/

fc = pupmod_loadpowcorr(v,SUBJLIST,1);

nvox = size(fc,1)*size(fc,1)-size(fc,1);

%%

% atx effect #1 (9.5-16 Hz), where p<0.05
fc(:,:,:,:,:,18) = tanh(nanmean(atanh(fc(:,:,:,:,:,[6 7 8 9])),6));
% dpz effect #1 (9.5-16 Hz), where p<0.05
fc(:,:,:,:,:,19) = tanh(nanmean(atanh(fc(:,:,:,:,:,[6 7 8 9 10])),6));
% task vs rest effect #1 (9.5-16 Hz), where p<0.05
fc(:,:,:,:,:,20) = tanh(nanmean(atanh(fc(:,:,:,:,:,[4 5 6])),6));

tmp = clock;

seed = ((tmp(1)+tmp(2)*tmp(3))/tmp(4)+tmp(5))*tmp(6);

rng(seed,'twister')

if ~exist(sprintf([outdir 'pupmod_src_powcorr_second_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v))
  all_idx1 = randi(2,[size(SUBJLIST,2),nperm]);
  save(sprintf([outdir 'pupmod_src_powcorr_second_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v),'all_idx1');
else
  load(sprintf([outdir 'pupmod_src_powcorr_second_permtest_perms_subs%d_nperm%d_v%d.mat'],par.subs,nperm,v));
end

dat_res1 = single(squeeze(fc(:,:,:,[1 2],1,:))); 
dat_cnt1 = single(squeeze(fc(:,:,:,[1 2],2,:)));  
dat_cnt2 = single(squeeze(fc(:,:,:,[1 3],2,:))); 
dat_res2 = single(squeeze(fc(:,:,:,[1 3],1,:))); 

taskvsrest(:,:,:,1,:) = dat_res1(:,:,:,1,:);
taskvsrest(:,:,:,2,:) = dat_cnt1(:,:,:,1,:);

for iperm = 1 : par.allperms
  
  if ~exist(sprintf([outdir 'pupmod_src_powcorr_second_permtest_iperm%d_nperm%d_v%d_processing.txt'],iperm,nperm,v))
    system(['touch ' outdir sprintf('pupmod_src_powcorr_second_permtest_iperm%d_nperm%d_v%d_processing.txt',iperm,nperm,v)]);
  else
    continue
  end
  
  for kperm = 1 : par.subs
    
    iiperm = (iperm-1)*par.subs+kperm;
    
    % within subjects permutation test 
    fprintf('Perm #%d\n',kperm);
    
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat_cnt1(:,:,i,1,:) = dat_cnt1(:,:,i,idx1(i),:);
      permdat_cnt1(:,:,i,2,:) = dat_cnt1(:,:,i,idx2(i),:);     
      permdat_res1(:,:,i,1,:) = dat_res1(:,:,i,idx1(i),:);
      permdat_res1(:,:,i,2,:) = dat_res1(:,:,i,idx2(i),:);
    end
    
    % should this be the same as for ATX??
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)
      permdat_cnt2(:,:,i,1,:) = dat_cnt2(:,:,i,idx1(i),:);
      permdat_cnt2(:,:,i,2,:) = dat_cnt2(:,:,i,idx2(i),:);
      permdat_res2(:,:,i,1,:) = dat_res2(:,:,i,idx1(i),:);
      permdat_res2(:,:,i,2,:) = dat_res2(:,:,i,idx2(i),:); 
    end

    % TASK VS REST ------------------------------------------------------
    idx1 = all_idx1(:,iiperm);
    idx2 = 3-idx1;
    
    for i = 1 : length(idx1)     
      taskvsrest_perm(:,:,i,1,:) = taskvsrest(:,:,i,idx1(i),:);
      taskvsrest_perm(:,:,i,2,:) = taskvsrest(:,:,i,idx2(i),:);       
    end
    
%     for ifoi = 1 : nfoi
      
      % -----------
      % compute ttest during task and atomoxetine
      [t_cnt1,~,~,s] = ttest(atanh(permdat_cnt1(:,:,:,2,:)),atanh(permdat_cnt1(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_cnt1 = squeeze(t_cnt1.*sign(s.tstat)); clear s
      % compute ttest during rest and atomoxetine
      [t_res1,~,~,s] = ttest(atanh(permdat_res1(:,:,:,2,:)),atanh(permdat_res1(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_res1 = squeeze(t_res1.*sign(s.tstat)); clear s
      % compute ttest during task and donepezil
      [t_cnt2,~,~,s] = ttest(atanh(permdat_cnt2(:,:,:,2,:)),atanh(permdat_cnt2(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_cnt2 = squeeze(t_cnt2.*sign(s.tstat)); clear s
      % compute ttest during rest and donepezil
      [t_res2,~,~,s] = ttest(atanh(permdat_res2(:,:,:,2,:)),atanh(permdat_res2(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_res2 = squeeze(t_res2.*sign(s.tstat)); clear s
      
      % compute double contrast, test connections for context-dependence
      [t_all1,~,~,s] = ttest(atanh(permdat_res1(:,:,:,2,:))-atanh(permdat_res1(:,:,:,1,:)),atanh(permdat_cnt1(:,:,:,2,:))-atanh(permdat_cnt1(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_all1 = squeeze(t_all1.*sign(s.tstat)); clear s
      % compute double contrast, test connections for context-dependence
      [t_all2,~,~,s] = ttest(atanh(permdat_res2(:,:,:,2,:))-atanh(permdat_res2(:,:,:,1,:)),atanh(permdat_cnt2(:,:,:,2,:))-atanh(permdat_cnt2(:,:,:,1,:)),'dim',3,'alpha',alp);
      t_all2 = squeeze(t_all2.*sign(s.tstat)); clear s
      % -----------------------
      % NUMBER OF ALTERED CORRELATONS - across space
      % -----------------------
      % count number of altered connections (atx, task)
      par.tperm_cnt1_n(kperm,:)=nansum(nansum(t_cnt1<0))./nvox;
      par.tperm_cnt1_p(kperm,:)=nansum(nansum(t_cnt1>0))./nvox;
      % count number of altered connections (atx, rest)
      par.tperm_res1_n(kperm,:)=nansum(nansum(t_res1<0))./nvox;
      par.tperm_res1_p(kperm,:)=nansum(nansum(t_res1>0))./nvox;
      % count number of altered connections (dpz, task)
      par.tperm_cnt2_n(kperm,:)=nansum(nansum(t_cnt2<0))./nvox;
      par.tperm_cnt2_p(kperm,:)=nansum(nansum(t_cnt2>0))./nvox;
      % count number of altered connections (dpz, rest)
      par.tperm_res2_n(kperm,:)=nansum(nansum(t_res2<0))./nvox;
      par.tperm_res2_p(kperm,:)=nansum(nansum(t_res2>0))./nvox;
     	% number of altered connections, irrespective of direction (atx)
      par.tperm_atx_during_task(kperm,:)=nansum(nansum(abs(t_cnt1)))./nvox;
      par.tperm_atx_during_rest(kperm,:)=nansum(nansum(abs(t_res1)))./nvox;
      % number of altered connections, irrespective of direction (dpz)
      par.tperm_dpz_during_task(kperm,:)=nansum(nansum(abs(t_cnt2)))./nvox;
      par.tperm_dpz_during_rest(kperm,:)=nansum(nansum(abs(t_res2)))./nvox;
      
      % -----------------------
      % NUMBER OF ALTERED CORRELATONS - per voxel
      % -----------------------
      siz=size(t_cnt1,1);
      % count number of altered connections (atx, task)
      par.tperm_cnt1_pervoxel_n(:,kperm,:)=nansum(t_cnt1<0)./siz;
      par.tperm_cnt1_pervoxel_p(:,kperm,:)=nansum(t_cnt1>0)./siz;
      % count number of altered connections (atx, rest)
      par.tperm_res1_pervoxel_n(:,kperm,:)=nansum(t_res1<0)./siz;
      par.tperm_res1_pervoxel_p(:,kperm,:)=nansum(t_res1>0)./siz;
      % count number of altered connections (dpz, task)
      par.tperm_cnt2_pervoxel_n(:,kperm,:)=nansum(t_cnt2<0)./siz;
      par.tperm_cnt2_pervoxel_p(:,kperm,:)=nansum(t_cnt2>0)./siz;
      % count number of altered connections (dpz, rest)
      par.tperm_res2_pervoxel_n(:,kperm,:)=nansum(t_res2<0)./siz;
      par.tperm_res2_pervoxel_p(:,kperm,:)=nansum(t_res2>0)./siz;
     	% number of altered connections, irrespective of direction (atx)
      par.tperm_atx_during_task_pervoxel(:,kperm,:)=nansum(abs(t_cnt1))./siz;
      par.tperm_atx_during_rest_pervoxel(:,kperm,:)=nansum(abs(t_res1))./siz;
      % number of altered connections, irrespective of direction (dpz)
      par.tperm_dpz_during_task_pervoxel(:,kperm,:)=nansum(abs(t_cnt2))./siz;
      par.tperm_dpz_during_rest_pervoxel(:,kperm,:)=nansum(abs(t_res2))./siz;
      
      % -----------------------
      % CONTEXT DEPENDENCE - across space
      % -----------------------
      % context dependence: diff in counts (atx) between task and rest
      par.tperm_context_diff_atx_n(kperm,:)= par.tperm_cnt1_n(kperm,:) - par.tperm_res1_n(kperm,:);
      par.tperm_context_diff_atx_p(kperm,:)= par.tperm_cnt1_p(kperm,:) - par.tperm_res1_p(kperm,:);
      % context dependence: diff in counts (dpz) between task and rest
      par.tperm_context_diff_dpz_n(kperm,:)= par.tperm_cnt2_n(kperm,:) - par.tperm_res2_n(kperm,:);
      par.tperm_context_diff_dpz_p(kperm,:)= par.tperm_cnt2_p(kperm,:) - par.tperm_res2_p(kperm,:);  
      % context-dependence: diff in counts, irrespective of direction 
      par.tperm_atx_context(kperm,:)=par.tperm_atx_during_task(kperm,:)-par.tperm_atx_during_rest(kperm,:);
      par.tperm_dpz_context(kperm,:)=par.tperm_dpz_during_task(kperm,:)-par.tperm_dpz_during_rest(kperm,:);     
       
      % -----------------------
      % CONTEXT DEPENDENCE - per voxel
      % -----------------------
      % context dependence: diff in counts (atx) between task and rest
      par.tperm_context_diff_atx_n_pervoxel(:,kperm,:) = par.tperm_cnt1_pervoxel_n(:,kperm,:) - par.tperm_res1_pervoxel_n(:,kperm,:);
      par.tperm_context_diff_atx_p_pervoxel(:,kperm,:) = par.tperm_cnt1_pervoxel_p(:,kperm,:) - par.tperm_res1_pervoxel_p(:,kperm,:);
      % context dependence: diff in counts (dpz) between task and rest
      par.tperm_context_diff_dpz_n_pervoxel(:,kperm,:) = par.tperm_cnt2_pervoxel_n(:,kperm,:) - par.tperm_res2_pervoxel_n(:,kperm,:);
      par.tperm_context_diff_dpz_p_pervoxel(:,kperm,:) = par.tperm_cnt2_pervoxel_p(:,kperm,:) - par.tperm_res2_pervoxel_p(:,kperm,:);

      % test double dissociation (atx vs. dpz, rest vs. task)
      a = par.tperm_context_diff_atx_n(kperm,:)+par.tperm_context_diff_atx_p(kperm,:);
      b = par.tperm_context_diff_dpz_n(kperm,:)+par.tperm_context_diff_dpz_p(kperm,:);
      par.tperm_doubledissociation(kperm,:) = a-b;
      
      % context-dependence: number of context-dependent connections 
      par.tperm_atx_context_test(kperm,:) = nansum(nansum(t_all1))./siz;
      par.tperm_dpz_context_test(kperm,:) = nansum(nansum(t_all2))./siz;
          
      % TASK VS REST ------------------------------------------------------
      [t_tvsr,~,~,s] = ttest(taskvsrest_perm(:,:,:,2,:),taskvsrest_perm(:,:,:,1,:),'dim',3,'alpha',alp);
      t_tvsr = squeeze(t_tvsr.*sign(s.tstat)); clear s   
      
      par.tperm_taskvsrest_p(kperm,:) = nansum(nansum(t_tvsr>0))./siz;
      par.tperm_taskvsrest_n(kperm,:) = nansum(nansum(t_tvsr<0))./siz;
      
      par.tperm_taskvsrest_pervox_p(:,kperm,:) = nansum(t_tvsr>0)./siz;
      par.tperm_taskvsrest_pervox_n(:,kperm,:) = nansum(t_tvsr<0)./siz;
      
%     end
  end
  
   save(sprintf([outdir 'pupmod_src_powcorr_second_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'par')
  
  try
    pause(randi(3))
    load(sprintf([outdir 'pupmod_src_powcorr_second_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v))
  catch me
    save(sprintf([outdir 'pupmod_src_powcorr_second_permtest_iperm%d_nperm%d_v%d.mat'],iperm,nperm,v),'par')
  end
  
end
%%
