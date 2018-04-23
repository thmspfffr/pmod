%% pupmod_src_fcmat_cortex
% PLOTS FC MATRICES, ORDERED FROM ANTERIOR TO POSTERIOR

clear

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab/
addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/

clear permdat_cnt1
clear permdat_cnt2
clear permdat_res1
clear permdat_res2
v = 12;
outdir   = '/home/tpfeffer/pupmod/proc/conn/';

addpath ~/pconn/matlab/
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

load sa_meg_template.mat

grid = select_chans(sa_meg_template.grid_cortex3000,400);

[~,ind]=sort(grid(:,2),'descend');

% SEPERATE HEMISPHERES
% left hemisphere
lh = find(grid(:,1)<0);
rh = find(grid(:,1)>0);

[~,ind_lh]=sort(grid(lh,2),'descend');
[~,ind_rh]=sort(grid(rh,2),'descend');

ind = [lh(ind_lh); rh(ind_rh)];

%%
ifoi = 6

fc_atx_rest = (squeeze(nanmean(cleandat(ind,ind,:,2,1,ifoi),3))-squeeze(nanmean(cleandat(ind,ind,:,1,1,ifoi),3)))./squeeze(nanmean(cleandat(ind,ind,:,2,1,ifoi),3));
fc_atx_task = (squeeze(nanmean(cleandat(ind,ind,:,2,2,ifoi),3))-squeeze(nanmean(cleandat(ind,ind,:,1,2,ifoi),3)))./squeeze(nanmean(cleandat(ind,ind,:,2,1,ifoi),3));
fc_dpz_rest = (squeeze(nanmean(cleandat(ind,ind,:,3,1,ifoi),3))-squeeze(nanmean(cleandat(ind,ind,:,1,1,ifoi),3)))./squeeze(nanmean(cleandat(ind,ind,:,2,1,ifoi),3));
fc_dpz_task = (squeeze(nanmean(cleandat(ind,ind,:,3,2,ifoi),3))-squeeze(nanmean(cleandat(ind,ind,:,1,2,ifoi),3)))./squeeze(nanmean(cleandat(ind,ind,:,2,1,ifoi),3));

fc_taskvsrest = squeeze(nanmean(cleandat(ind,ind,:,1,2,ifoi),3))-squeeze(nanmean(cleandat(ind,ind,:,1,1,ifoi),3));
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','white')
subplot(1,2,1)
fc_atx = rot90(triu(fc_atx_rest,1),2)+flipud(rot90(triu(fc_atx_task,1),-1));
imagesc(fc_atx,[-0.3 0.3]); colormap(cmap); axis square off

subplot(1,2,2)
fc_dpz = rot90(triu(fc_dpz_rest,1),2)+flipud(rot90(triu(fc_dpz_task,1),-1));
imagesc(fc_dpz,[-0.3 0.3]); colormap(cmap); axis square off

print(gcf,'-dpng',sprintf('~/pmod/plots/pupmod_src_fcmat_cortex_f%d.png',ifoi))



%% plot brain 

mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = 1;
% m2 = spatfiltergauss(m,g1,dd,g2);
% deg2plot(i,7) = 0;
z2 = spatfiltergauss(grid(:,2),g1,dd,g2);

para = [] ;
para.colorlimits = [min(z2) max(z2)];
% cmap = cbrewer('seq', 'YlOrRd', 100,'pchip'); %cmap(end:-1:1,:);
cmap = plasma;
cmap = cmap(end:-1:1,:);
figure; set(gcf,'color','white');
viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
% cmap = viridis; 
para.filename = '~/pmod/plots/pupmod_src_fcmat_cortex_gradient.png';
tp_showsource(z2,cmap,sa_meg_template,para);
