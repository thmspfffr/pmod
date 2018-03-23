%% pmod_WC_wholebrain_rest
% Stochastic simulation of 2*N WC nodes during "rest"
%-------------------------------------------------------------------------

clear

%-------------------------------------------------------------------------
% VERSION 1
%-------------------------------------------------------------------------
% v = 1;
% Ies = -1.5:0.2:2;
% Iis = -4:0.2:0;
% nTrials = 20;
%-------------------------------------------------------------------------
% VERSION 2
%-------------------------------------------------------------------------
% v = 2;
% Ies = -2.2:0.1:0.1;
% Iis = -3.8:0.2:-0.8;
% nTrials = 20;
%-------------------------------------------------------------------------
% VERSION 3
%-------------------------------------------------------------------------
% v = 3;
% Ies = -2.2:0.1:0.1;
% Iis = -4;
% Gg = 0.5:0.02:0.8; % 0.62
% nTrials = 20;
%-------------------------------------------------------------------------
% VERSION 4
%-------------------------------------------------------------------------
v = 4;
Ies = -3.5:0.2:0.8;
Iis = -5:0.2:-1;
Gg = 0.4:0.1:0.9; % 0.62
nTrials = 3;
%-------------------------------------------------------------------------
% VERSION 5
%-------------------------------------------------------------------------
% v = 5;
% Ies = -3:0.2:1;
% Iis = -4:0.2:-2;
% Gg = 0.5:0.05:0.8; % 0.62
% nTrials = 3;
%-------------------------------------------------------------------------


% LOAD CLEANED DATA AND COMPUTE MEAN
% v = 1;
% load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
% fc_emp_atx(:,:,1) = squeeze(nanmean(cleandat(:,:,:,2,1,6),3))-squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
% fc_emp_atx(:,:,2) = squeeze(nanmean(cleandat(:,:,:,2,2,6),3))-squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
% fc_emp_dpz(:,:,1) = squeeze(nanmean(cleandat(:,:,:,3,1,7),3))-squeeze(nanmean(cleandat(:,:,:,1,1,7),3));
% fc_emp_dpz(:,:,2) = squeeze(nanmean(cleandat(:,:,:,3,2,7),3))-squeeze(nanmean(cleandat(:,:,:,1,2,7),3));

% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));
N = size(C,1);

%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;


tauE = 1;
tauI = 2;
tau = zeros(2*N,1);
tau(1:N) = tauE;
tau(N+1:2*N) = tauI;

dt=0.01;
tmax = 10000; % in units of tauE
tspan=0:dt:tmax;
L = length(tspan);

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

sigma = 0.0005;
%Qn = (sigma*dt)^2*eye(2*N);

% Connectivity:
% dI_drug = zeros(2*N,1);
% dI_drug(1:N) = dIe_drug;
% dI_drug(N+1:2*N) = dIi_drug;

isub = find( triu(ones(N)) - eye(N) );
%%
for iies = 1 : length(Ies)
  for iiis = 1 : length(Iis)
    for iG = 1 : length(Gg)
      
      if ~exist(sprintf(['~/pmod/proc/' 'pmod_WC_wholebrain_rest_Ie%d_Ii%d_G%d_v%d_processing.txt'],iies,iiis,iG,v))
        system(['touch ' '~/pmod/proc/' sprintf('pmod_WC_wholebrain_rest_Ie%d_Ii%d_G%d_v%d_processing.txt',iies,iiis,iG,v)]);
      else
        continue
      end
      
      g = Gg(iG);
      W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
      
      FC = zeros(N,N,1);
      %     Rate = zeros(1,1);
      %     RateSD = zeros(2,1);
      
      %     Cee = zeros(1,1);
      %     CeeSD = zeros(2,1);
      
      %     FCval = cell(1,1);
      %     Epsilon = cell(1,1);
      
      KOPsd   = zeros(nTrials,1);
      KOPmean = zeros(nTrials,1);
      % Control params.
      %--------------------
      Ie = Ies(iies);
      Ii = Iis(iiis);
      
      Io=zeros(2*N,1);
      Io(1:N) = Ie;
      Io(N+1:2*N) = Ii;
      
      % transfer function:
      gE  = 1;
      gI  = 1;
      aE  = 1/gE;
      aI  = 1/gI;
      Fe  = @(x) 1./(1 + exp(-x/aE) );
      Fi  = @(x) 1./(1 + exp(-x/aI) );
      F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
      
      % Working point:
      Io(1:N)     = Ie;
      Io(N+1:2*N) = Ii;
      %
      %     FCval{1}    = zeros(length(isub),nTrials);
      %     Epsilon{1}  = zeros(N,nTrials);
      
      T       = Tds*resol; %% define time of interval
      %     freqs   = (0:Tds/2)/T; %% find the corresponding frequency in Hz
      %     nfreqs  = length(freqs);
      %     freq100 = freqs(freqs<100 & freqs>1);
      %     pp      = 1:10:length(freq100);
      %     PSD     = zeros(length(pp),N,nTrials);
      
      %     PSDstruct(1).frequencies = freq100(1:10:end)';
      
      
      for tr=1:nTrials
        fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
        r   = 0.001*rand(2*N,1);
        R   = zeros(Tds,N);
        Ri  = zeros(Tds,N);
        tt  = 0;
        % transient:
        for t = 1:50000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(2*N,1);
        end
        % simulation
        for t = 1:L
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(2*N,1);
          if mod(t,ds)==0
            tt=tt+1;
            R(tt,:)  = r(1:N);
            Ri(tt,:)  = r(N+1:end);
          end
        end
        
        %correlation:
        rE = R;
        rI = Ri;
        
        z             = rE + 1i*rI;
        ku            = sum(z,2)/N;
        KOP           = abs(ku);
        KOPsd(tr,1)   = std(KOP);
        KOPmean(tr,1) = mean(KOP);
        tmp = tp_dfa(R,[3 50],ds,0.5,15);
        dfa(tr,:,1) = tmp.exp;
        
        rc              = corrcoef(rE);
        FC(:,:,1)       = FC(:,:,1) + rc/nTrials;
        %       fc              = rc(isub);
        %       Cee(1)          = Cee(1) + mean(fc)/nTrials;
        %       CeeSD(1)        = CeeSD(1) + var(fc)/nTrials;
        %       FCval{1}(:,tr)  = fc;
        
        %       rEo       = mean(mean(rE));
        %       rEsd      = var(mean(rE));
        %       Rate(1)   = rEo-rEo;
        %       RateSD(1) = RateSD(1) + rEsd/nTrials;
        
        %       for i=1:N
        %         %autocorr
        % %         [acf,lags] = autocorr(R(:,i),300);
        % %         ii = find(acf<.2,1,'first');
        % %         Epsilon{1}(i,tr) = lags(ii)*resol;
        %         %PSD:
        %         f = R(:,i) - mean(R(:,i));
        %         xdft = fft(f);
        %         xdft = xdft(1:floor(Tds/2)+1);
        %         pw = (1/(Tds/2)) * abs(xdft).^2;
        %         psd = pw(freqs<100 & freqs>1);
        %         f = freqs(freqs<100 & freqs>1);
        %         fnew = f(1:10:end);
        %         psd  = psd(1:10:end);
        %         %       pw = gaussfilt(fnew,psd',.5);
        %         PSD(:,i,tr) = psd';
        %       end
        
      end
      
      %     CeeSD(1)  = sqrt(CeeSD(1));
      %     RateSD(1) = sqrt(RateSD(1));
      %
      %     PSDstruct(1).PSD = PSD;
      %
      out = struct('dfa',dfa,'FC',FC,'Ie',Ie,'Ii',Ii,'G',g,'KOP',KOP,'KOPsd',KOPsd,'KOPmean',KOPmean,'ku',ku);
      
      save(sprintf('~/pmod/proc/pmod_WC_wholebrain_rest_Ie%d_Ii%d_G%d_v%d.mat',iies,iiis,iG,v),'out')
      
    end
  end
end

error('!')

%%

load(sprintf(['~/pconn/proc/dfa/' 'pconn_src_dfa_aal_f%d_m%d_v%d.mat'],2,1,2),'dfa');
dfa_emp = dfa; clear dfa

clear dfa r dfa_r dist_fc fc_sim
v=1;
vv = 4;
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
mask    = logical(tril(ones(90,90),-1));
fc_rest =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));

for iies = 1 : length(Ies)
  for iiis = 1 : length(Iis)
    for iG = 1 : length(Gg)
      
%       try
      load(sprintf('~/pmod/proc/pmod_WC_wholebrain_rest_Ie%d_Ii%d_G%d_v%d.mat',iies,iiis,iG,vv))
      dfa(:,iies,iiis,iG) = squeeze(mean(out.dfa));
      pars = [];
      pars.grid = 'medium';
      pars.N = 90;
      
      fc_diff = tp_match_aal(pars,out.FC(:,:,1));
      r(iies,iiis,iG)=corr(fc_diff(mask),fc_rest(mask));
      

      
      idx = find(~isnan(dfa_emp(1:90)'));
      dfa_r (iies,iiis,iG) = corr(dfa_emp(idx)',dfa(idx,iies,iiis,iG));
      dist_fc (iies,iiis,iG)  = abs(mean(fc_diff(mask))-mean(fc_rest(mask)));
      fc_sim(iies,iiis,iG) = mean(fc_diff(mask));
      
%       catch 
%         dist_fc(iies,iiis,iG)  = nan;
%         fc_sim(iies,iiis,iG)  = nan;
%         dfa_r (iies,iiis,iG) = nan;
%         r(iies,iiis,iG) = nan;
%       end
    end
  end
end

% pars.N = 90;
% dfa_emp  =  tp_match_aal(pars,dfa_emp');

%%
j = 24;
i = 4;
h = figure; set(gcf,'color','white'); hold on
set(h,'Position',[10,10,1200,400])

gg = 1;

subplot(1,4,1);
% imagesc(flipud(squeeze(r)'),[0 0.25]); axis square tight
imagesc(flipud(squeeze(r(:,:,gg))'),[0 0.4]); axis square tight
% contourc(Ies,Gg,squeeze(r)',[0 0.25]); axis square tight


if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:2:length(Iis),'yticklabel',Iis(length(Iis):-2:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Inh. input I_I');
end

title('Correlation with FC')
colormap(inferno)

tp_editplots; hold on
scatter(j,length(Gg)-i+1,20,'markerfacecolor','w','markeredgecolor','k')

h=subplot(1,4,2);

imagesc(flipud(squeeze(mean(dfa(:,:,:,gg),1))'),[0.5 1]);  axis square tight


if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:2:length(Iis),'yticklabel',Iis(length(Iis):-2:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Inh. input I_I');
end

colormap(inferno)
title('DFA')
tp_editplots
hold on
scatter(j,length(Gg)-i+1,20,'markerfacecolor','w','markeredgecolor','k')

%
% subplot(1,4,2);
% figure; set(gcf,'color','white')
%
% imagesc(squeeze(dist_fc)',[0 0.5]); axis square
%
% if length(Iis)== 1
%   set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
%   set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(1:4:length(Gg)),'tickdir','out')
%   xlabel('Exc. input I_E'); ylabel('Global coupling G');
% else
%   set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
%   set(gca,'ytick',1:3:length(Iis),'yticklabel',Iis(1:4:length(Gg)),'tickdir','out')
%   xlabel('Exc. input I_E'); ylabel('Inh. input I_I');
% end
% title('Distance FC')
% tp_editplots
%
% figure; set(gcf,'color','white'); hold on
% 
k=subplot(1,4,3);
imagesc(flipud(squeeze(dfa_r(:,:,gg))'),[-1 0.4]);  axis square tight

if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:2:length(Iis),'yticklabel',Iis(length(Iis):-2:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Inh. input I_I');
end

colormap(inferno)
title('Correlation with DFA')
tp_editplots
hold on

scatter(j,length(Gg)-i+1,20,'markerfacecolor','w','markeredgecolor','k')
% set(k,'Position',get(k,'Position')+0.1)


l = subplot(1,4,4);
% 
idx = squeeze(mean(dfa(:,:,:,gg)))>0 & squeeze(mean(dfa(:,:,:,gg)))<2;
par = (log(abs(squeeze(r(:,:,gg)./(max(r(:))))'))+log(abs(squeeze(dfa_r(:,:,gg)./max(dfa_r(:)))))')/2;
par(~idx') = -4;

imagesc(flipud(par),[-4 0]);  axis square tight

if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:2:length(Iis),'yticklabel',Iis(length(Iis):-2:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Inh. input I_I');
end

colormap(inferno)
title('Correlation with DFA')
tp_editplots
hold on
[i,j]=find(par==max(par(:)));
scatter(j,length(Gg)-i+1,20,'markerfacecolor','w','markeredgecolor','k')
% set(l,'Position',get(l,'Position')+0.1)

% STUFF
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_WC_wholebrain_rest4_v%d.pdf',v))

%%
figure; hold on
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[5 5 800 1500])
set(gcf,'Position',[50 50 400 800],'Color','w')
imagesc(abs(0.67-squeeze(mean(dfa,1)))',[0 0.2]); axis square

set(gcf,'Position',[150 150 400 1300],'Color','w')
imagesc(abs(0.67-squeeze(mean(dfa,1)))',[0 0.2]); axis square
