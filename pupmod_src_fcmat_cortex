%% pupmod_src_powcorr_degree
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

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

f = fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
aal_labels = textscan(f,'%d %s %d','headerlines',0);

load aalmask_grid_coarse.mat

if v==1
  for i = 1 : 90
    for j = 1 : 90

      idx = find(aalgrid.mask==i,1,'first');
      lab = aalgrid.labels(idx);
      idx = find(strcmp(aal_labels{2},lab));

      jdx = find(aalgrid.mask==j,1,'first');
      lab = aalgrid.labels(jdx);
      jdx = find(strcmp(aal_labels{2},lab));

      p(i,j,:,:,:,:) = cleandat(idx,jdx,:,:,:,:);

    end
  end

  cleandat = p;
end

for ifoi = 6:7
  %     ifoi
  
  s_fc(:,:,:,:,1,ifoi) = cleandat(:,:,:,:,1,ifoi);
  s_fc(:,:,:,:,2,ifoi) = cleandat(:,:,:,:,2,ifoi);
  
end

%%

clear z1 z2 th1 th2 th c deg z z1 z2 p1 p2

% reference (1) drug (1-3), (2) rest/task (1/2)
ref_cond    = [1 1];
test_cond   = [1 1];

j = 1 : size(s_fc,1);

m1 = squeeze(nanmean(s_fc(j,:,:,ref_cond(1),ref_cond(2),:),2));
s1 = squeeze(nanstd(s_fc(j,:,:,ref_cond(1),ref_cond(2),:),[],2));
m2 = squeeze(nanmean(s_fc(:,j,:,ref_cond(1),ref_cond(2),:),1));
s2 = squeeze(nanstd(s_fc(:,j,:,ref_cond(1),ref_cond(2),:),[],1));

for ifoi = 6:6
  for il = 1 : size(s_fc,1)
    
    fprintf('Computing location %d ...\n',il);
    
    jj = j(j~=il);
    
    
    for jl = 1:size(s_fc,1)
       
      x = squeeze(s_fc(il,jl,:,test_cond(1),test_cond(2),ifoi));
      
      z1 = (x' - squeeze(m1(jl,:,ifoi)))./squeeze(s1(jl,:,ifoi));
      z2 = (x' - squeeze(m2(jl,:,ifoi)))./squeeze(s2(jl,:,ifoi));
      
      [~,p1] = ttest(z1');
      [~,p2] = ttest(z2');
      
      th1(jl,:) = p1 < 0.01/2;
      th2(jl,:) = p2 < 0.01/2;
      
    end
    
    th1 = th1(jj,:);
    th2 = th2(jj,:);
    
    % any connection significant?
    th(il,jj,:) = (th1 + th2) > 0; clear th1 th2
    
  end
  c(ifoi)= sum(th(:))/(size(th,1)*size(th,1));
  deg(:,ifoi) = mean(th);
  
end


    %%
% figure;
subplot(3,2,3)
hold on
% plot(c,'color',[0.95 0.95 0.2],'linewidth',3)
area(c)
set(gca,'xtick',[1 3 5 7 9 11 13],'xticklabel',[2 4 8 16 32 64 128])
xlabel('Frequency [Hz]'); ylabel('Number of connections [%]'); title('Degree');
box off
% % set(gca,'linewidth',.5,'color','w')
% axis square
axis([1 14 0 0.4]); set(gca,'tickdir','out','ytick',[0 0.1 0.2 0.3 0.4],'yticklabel',[0 0.1 0.2 0.3 0.4],'color','w')

print(gcf,'-dpdf',sprintf('~/pupmod/plots/pupmod_src_degree.pdf'))

subplot(3,2,2); hold on; plot(1:13,c,'linewidth',3,'color',[0 0.5 1]);



%%
%% PLOT NUMBER OF ALTERED CONNECTIONS (IRRESPECTIVE OF SIGN)
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
load sa_meg_template;

if v == 1
  for ifoi = 1 : 13
    ifoi
    for i = 1 : 2113

      idx  = aalgrid.mask(i);
      if idx == 0 || idx > 90
        deg2plot(i,ifoi)=0;
      else
        deg2plot(i,ifoi) =deg(idx,ifoi);

      end

    end
  end
  
else 
  for ifoi = [6]
    deg2plot(:,ifoi) = deg(:,ifoi);
    load('~/pupmod/proc/grid_cortexlow.mat')
  end
end
  
mri   = sa_meg_template.mri;
vc    = sa_meg_template.vc;
g1    = grid;
g2    = sa_meg_template.cortex10K.vc;
dd    = 1;
% m2 = spatfiltergauss(m,g1,dd,g2);
% deg2plot(i,7) = 0;
z2 = spatfiltergauss(deg2plot(:,6),g1,dd,g2);

para = [] ;
para.colorlimits = [.30 .7];

%
figure; set(gcf,'color','white');
viewdir = [-.5 -.5 .5; .5 -.5 .5; .5 .5 .5; -.5 .5 .5];
cmap = inferno;

cmap = inferno; cmap = cmap(102:end,:);
gray = [1 1 1]; 

len = 50;

g= [linspace(gray(1), cmap(1,1),len)', linspace(gray(2), cmap(1,2),len)', linspace(gray(3), cmap(1,3),len)']

cmap = [g;cmap];


tp_showsource(z2,cmap,sa_meg_template,para);
print(gcf,'-djpeg100',sprintf('~/pupmod/plots/pupmod_src_powcorr_maps_degree_f%d_v%d.jpg',ifoi,v))

%%

for ifoi = 1 : 13

  atx = triu(nanmean(s_fc(:,:,:,2,2,ifoi),3)-nanmean(s_fc(:,:,:,1,2,ifoi),3),1);
  atx = atx(:); atx(atx==0)=[];

  dpz = triu(nanmean(s_fc(:,:,:,3,1,ifoi),3)-nanmean(s_fc(:,:,:,1,1,ifoi),3),1);
  dpz = dpz(:); dpz(dpz==0)=[];

  r(ifoi) = corr(atx,dpz);

end
%%
for i = 1 : 400
  for j = 1 : 400
    d(i,j) = sqrt([grid(i,1)-grid(j,1)].^2+[grid(i,2)-grid(j,2)].^2+[grid(i,3)-grid(j,3)].^2);
  end
  
end
    

