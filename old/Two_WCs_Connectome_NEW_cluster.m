
% Stochastic simulation of 2*N WC nodes subject to drug at t=0

%-------------------------------------------------------------------------


clear all;

% load connectome
load ~/pmod/matlab/EC %Matt_EC
C = EC;
C = C/max(C(C>0));
N = size(C,1);


% fixed params:
nTrials = 20;
% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

g = .62;

tauE = 1;
tauI = 2;
tau = zeros(2*N,1);
tau(1:N) = tauE;
tau(N+1:2*N) = tauI;

dt=0.01;
tmax = 20000; % in units of tauE
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
W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];


% Control params.
%--------------------
Ie = -2.85;
Ii = -3.50; %-1;

% input due to task:
dIe = 1.0;
dIi = 1.2; %1.3;

Io=zeros(2*N,1);
Io(1:N) = Ie;
Io(N+1:2*N) = Ii;

dI=zeros(2*N,1);
dI(1:N) = dIe;
dI(N+1:2*N) = dIi;

for Drug = 1:2 % Loop over drugs ------------------------------------------
  
  % Drug effects:
  if Drug==1
    % inputs
    dIe_drug = 0.01;%0.08;
    dIi_drug = -.08;%0;
    % gain modulation:
    Gain_E = .20;
    Gain_I = .20;
  else
    % inputs
    dIe_drug = -0.25;
    dIi_drug = -0.25;
    % gain modulation:
    Gain_E = .10; %0.07;
    Gain_I = .10; %0.07;
  end
  dI_drug = zeros(2*N,1);
  dI_drug(1:N) = dIe_drug;
  dI_drug(N+1:2*N) = dIi_drug;
  
  Rate = zeros(4,1);
  RateSD = zeros(4,1);
  
  Cee = zeros(4,1);
  CeeSD = zeros(4,1);
  
  FC = zeros(N,N,4);
  FCval = cell(1,4);
  Epsilon = cell(1,4);
  
  KOPsd   = zeros(nTrials,4);
  KOPmean = zeros(nTrials,4);
  
  
  isub = find( triu(ones(N)) - eye(N) );
  
  display('REST ...')
  %--------------------------------------------------------------------------
  
  % transfer function:
  gE = 1;
  gI = 1;
  aE = 1/gE;
  aI = 1/gI;
  Fe = @(x) 1./(1 + exp(-x/aE) );
  Fi = @(x) 1./(1 + exp(-x/aI) );
  F = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
  
  % Working point:
  Io(1:N) = Ie;
  Io(N+1:2*N) = Ii;
  
  FCval{1} = zeros(length(isub),nTrials);
  Epsilon{1} = zeros(N,nTrials);
  
  T = Tds*resol; %% define time of interval
  freqs = (0:Tds/2)/T; %% find the corresponding frequency in Hz
  nfreqs=length(freqs);
  freq100 = freqs(freqs<100 & freqs>1);
  pp = 1:10:length(freq100);
  PSD = zeros(length(pp),N,nTrials);
  
  PSDstruct(1).frequencies = freq100(1:10:end)';
  
  
  for tr=1:nTrials
    display(sprintf('trial: %g',tr))
    r = 0.001*rand(2*N,1);
    R = zeros(Tds,N);
    Ri = zeros(Tds,N);
    tt = 0;
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
    
    z = rE + 1i*rI;
    ku=sum(z,2)/N;
    KOP=abs(ku);
    KOPsd(tr,1) = std(KOP);
    KOPmean(tr,1) = mean(KOP);
    
    rc = corrcoef(rE);
    FC(:,:,1) = FC(:,:,1) + rc/nTrials;
    fc = rc(isub);
    Cee(1) = Cee(1) + mean(fc)/nTrials;
    CeeSD(1) = CeeSD(1) + var(fc)/nTrials;
    FCval{1}(:,tr) = fc;
    
    rEo = mean(mean(rE));
    rEsd = var(mean(rE));
    Rate(1) = rEo-rEo;
    RateSD(1) = RateSD(1) + rEsd/nTrials;
    
    for i=1:N
      %autocorr
      [acf,lags] = autocorr(R(:,i),300);
      ii = find(acf<.2,1,'first');
      Epsilon{1}(i,tr) = lags(ii)*resol;
      %PSD:
      f = R(:,i) - mean(R(:,i));
      xdft = fft(f);
      xdft = xdft(1:floor(Tds/2)+1);
      pw = (1/(Tds/2)) * abs(xdft).^2;
      psd = pw(freqs<100 & freqs>1);
      f = freqs(freqs<100 & freqs>1);
      fnew = f(1:10:end);
      psd  = psd(1:10:end);
%       pw = gaussfilt(fnew,psd',.5);
%       PSD(:,i,tr) = pw;
    end
    
  end
  
  CeeSD(1)  = sqrt(CeeSD(1));
  RateSD(1) = sqrt(RateSD(1));
  
%   PSDstruct(1).PSD = PSD;
  
  
  display('REST + drug ...')
  %--------------------------------------------------------------------------
  
  % transfer function:
  gE = 1 + Gain_E;
  gI = 1 + Gain_I;
  aE = 1/gE;
  aI = 1/gI;
  Fe = @(x) 1./(1 + exp(-x/aE) );
  Fi = @(x) 1./(1 + exp(-x/aI) );
  F = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
  
  % Working point:
  Io(1:N) = Ie + dIe_drug;
  Io(N+1:2*N) = Ii+ dIi_drug;
  
  FCval{2} = zeros(length(isub),nTrials);
  
  T = Tds*resol; %% define time of interval
  freqs = (0:Tds/2)/T; %% find the corresponding frequency in Hz
  nfreqs=length(freqs);
  freq100 = freqs(freqs<100 & freqs>1);
  pp = 1:10:length(freq100);
  PSD = zeros(length(pp),N,nTrials);
  
  PSDstruct(2).frequencies = freq100(1:10:end)';
  
  
  for tr=1:nTrials
    display(sprintf('trial: %g',tr))
    r = 0.001*rand(2*N,1);
    R = zeros(Tds,N);
    Ri = zeros(Tds,N);
    tt = 0;
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
    
    z = rE + 1i*rI;
    ku=sum(z,2)/N;
    KOP=abs(ku);
    KOPsd(tr,2) = std(KOP);
    KOPmean(tr,2) = mean(KOP);
    
    rc = corrcoef(rE);
    FC(:,:,2) = FC(:,:,2) + rc/nTrials;
    fc = rc(isub);
    Cee(2) = Cee(2) + mean(fc)/nTrials;
    CeeSD(2) = CeeSD(2) + var(fc)/nTrials;
    FCval{2}(:,tr) = fc;
    
    r = mean(mean(rE));
    rEsd = var(mean(rE));
    Rate(2) = Rate(2) + (r-rEo)/nTrials;
    RateSD(2) = RateSD(2) + rEsd/nTrials;
    
    for i=1:N
      %autocorr
      [acf,lags] = autocorr(R(:,i),300);
      ii = find(acf<.2,1,'first');
      Epsilon{2}(i,tr) = lags(ii)*resol;
      %PSD:
      f = R(:,i) - mean(R(:,i));
      xdft = fft(f);
      xdft = xdft(1:floor(Tds/2)+1);
      pw = (1/(Tds/2)) * abs(xdft).^2;
      psd = pw(freqs<100 & freqs>1);
      f = freqs(freqs<100 & freqs>1);
      fnew = f(1:10:end);
      psd  = psd(1:10:end);
%       pw = gaussfilt(fnew,psd',.5);
%       PSD(:,i,tr) = pw;
    end
    
  end
  
  CeeSD(2)  = sqrt(CeeSD(2));
  RateSD(2) = sqrt(RateSD(2));
  
%   PSDstruct(2).PSD = PSD;
  
  
  display('TASK ...')
  %--------------------------------------------------------------------------
  
  % transfer function:
  gE = 1;
  gI = 1;
  aE = 1/gE;
  aI = 1/gI;
  Fe = @(x) 1./(1 + exp(-x/aE) );
  Fi = @(x) 1./(1 + exp(-x/aI) );
  F = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
  
  % Working point:
  Io(1:N) = Ie + dIe;
  Io(N+1:2*N) = Ii + dIi;
  
  FCval{3} = zeros(length(isub),nTrials);
  
  T = Tds*resol; %% define time of interval
  freqs = (0:Tds/2)/T; %% find the corresponding frequency in Hz
  nfreqs=length(freqs);
  freq100 = freqs(freqs<100 & freqs>1);
  pp = 1:10:length(freq100);
  PSD = zeros(length(pp),N,nTrials);
  
  PSDstruct(3).frequencies = freq100(1:10:end)';
  
  
  for tr=1:nTrials
    display(sprintf('trial: %g',tr))
    r = 0.001*rand(2*N,1);
    R = zeros(Tds,N);
    Ri = zeros(Tds,N);
    tt = 0;
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
    
    z = rE + 1i*rI;
    ku=sum(z,2)/N;
    KOP=abs(ku);
    KOPsd(tr,3) = std(KOP);
    KOPmean(tr,3) = mean(KOP);
    
    rc = corrcoef(rE);
    FC(:,:,3) = FC(:,:,3) + rc/nTrials;
    fc = rc(isub);
    Cee(3) = Cee(3) + mean(fc)/nTrials;
    CeeSD(3) = CeeSD(3) + var(fc)/nTrials;
    FCval{3}(:,tr) = fc;
    
    r = mean(mean(rE));
    rEsd = var(mean(rE));
    Rate(3) = Rate(3) + (r-rEo)/nTrials;
    RateSD(3) = RateSD(3) + rEsd/nTrials;
    
    for i=1:N
      %autocorr
      [acf,lags] = autocorr(R(:,i),300);
      ii = find(acf<.2,1,'first');
      Epsilon{3}(i,tr) = lags(ii)*resol;
      %PSD:
      f = R(:,i) - mean(R(:,i));
      xdft = fft(f);
      xdft = xdft(1:floor(Tds/2)+1);
      pw = (1/(Tds/2)) * abs(xdft).^2;
      psd = pw(freqs<100 & freqs>1);
      f = freqs(freqs<100 & freqs>1);
      fnew = f(1:10:end);
      psd  = psd(1:10:end);
%       pw = gaussfilt(fnew,psd',.5);
%       PSD(:,i,tr) = pw;
    end
    
  end
  
  CeeSD(3)  = sqrt(CeeSD(3));
  RateSD(3) = sqrt(RateSD(3));
  
%   PSDstruct(3).PSD = PSD;
  
  display('TASK + drug ...')
  %--------------------------------------------------------------------------
  
  % transfer function:
  gE = 1 + Gain_E;
  gI = 1 + Gain_I;
  aE = 1/gE;
  aI = 1/gI;
  Fe = @(x) 1./(1 + exp(-x/aE) );
  Fi = @(x) 1./(1 + exp(-x/aI) );
  F = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
  
  % Working point:
  Io(1:N) = Ie + dIe_drug + dIe;
  Io(N+1:2*N) = Ii + dIi_drug + dIi;
  
  FCval{4} = zeros(length(isub),nTrials);
  
  T = Tds*resol; %% define time of interval
  freqs = (0:Tds/2)/T; %% find the corresponding frequency in Hz
  nfreqs=length(freqs);
  freq100 = freqs(freqs<100 & freqs>1);
  pp = 1:10:length(freq100);
  PSD = zeros(length(pp),N,nTrials);
  
  PSDstruct(4).frequencies = freq100(1:10:end)';
  
  for tr=1:nTrials
    display(sprintf('trial: %g',tr))
    r = 0.001*rand(2*N,1);
    R = zeros(Tds,N);
    Ri = zeros(Tds,N);
    tt = 0;
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
    
    z = rE + 1i*rI;
    ku=sum(z,2)/N;
    KOP=abs(ku);
    KOPsd(tr,4) = std(KOP);
    KOPmean(tr,4) = mean(KOP);
%     dfa = tp_dfa(KOP,[3 50],1/dt,0.5,15);
%     KOPdfa(tr,4) = dfa.exp;
    
    rc = corrcoef(rE);
    FC(:,:,4) = FC(:,:,4) + rc/nTrials;
    fc = rc(isub);
    Cee(4) = Cee(4) + mean(fc)/nTrials;
    CeeSD(4) = CeeSD(4) + var(fc)/nTrials;
    FCval{4}(:,tr) = fc;
    
    r = mean(mean(rE));
    rEsd = var(mean(rE));
    Rate(4) = Rate(4) + (r-rEo)/nTrials;
    RateSD(4) = RateSD(4) + rEsd/nTrials;
    
    
    for i=1:N
      %autocorr
      [acf,lags] = autocorr(R(:,i),300);
      ii = find(acf<.2,1,'first');
      Epsilon{4}(i,tr) = lags(ii)*resol;
      %PSD:
      f = R(:,i) - mean(R(:,i));
      xdft = fft(f);
      xdft = xdft(1:floor(Tds/2)+1);
      pw = (1/(Tds/2)) * abs(xdft).^2;
      psd = pw(freqs<100 & freqs>1);
      f = freqs(freqs<100 & freqs>1);
      fnew = f(1:10:end);
      psd  = psd(1:10:end);
%       pw = gaussfilt(fnew,psd',.5);
%       PSD(:,i,tr) = pw;
    end   
  end % -------- END OF TRIAL LOOP --------
  % ---------------------------------------
  
  CeeSD(4)  = sqrt(CeeSD(4));
  RateSD(4) = sqrt(RateSD(4));
  
%   PSDstruct(4).PSD = PSD;
  
  %--------------------------------------------------------------------------
  % Save:
  
  if Drug ==1
    save('WC_connectome_ATX_AAL','Rate','RateSD','Cee','CeeSD','FC','FCval',...
      'Ie','Ii','dIe','dIi','dIe_drug','dIi_drug','Gain_E','Gain_I','Epsilon','KOPsd','KOPmean')
  else
    save('WC_connectome_DPZ_AAL','Rate','RateSD','Cee','CeeSD','FC','FCval',...
      'Ie','Ii','dIe','dIi','dIe_drug','dIi_drug','Gain_E','Gain_I','Epsilon','KOPsd','KOPmean')
  end
  
end % End loop over drugs -------------------------------------------------
%%
nTrials = 20;
figure; set(gcf,'color','white')

subplot(1,2,1); hold on
load ~/pmod/matlab/old/WC_connectome_ATX_AAL

bar(Cee); axis([0 5 0 0.02]); 
line([1 1],[Cee(1)-CeeSD(1)/sqrt(nTrials) Cee(1)+CeeSD(1)/sqrt(nTrials)])
line([2 2],[Cee(2)-CeeSD(2)/sqrt(nTrials) Cee(2)+CeeSD(2)/sqrt(nTrials)])
line([3 3],[Cee(3)-CeeSD(3)/sqrt(nTrials) Cee(3)+CeeSD(3)/sqrt(nTrials)])
line([4 4],[Cee(4)-CeeSD(4)/sqrt(nTrials) Cee(4)+CeeSD(4)/sqrt(nTrials)])
axis square
tp_editplots

fprintf('Altered connections: %.2f%%\n',100*sum(ttest(FCval{2},FCval{1},'dim',2,'alpha',0.01))/size(FCval{1},1))
fprintf('Corr > 0: %.2f%%\n',100*sum(ttest(FCval{2},FCval{1},'dim',2,'alpha',0.01).*(mean(FCval{2},2)-mean(FCval{1},2))>0)./sum(ttest(FCval{1},FCval{2},'dim',2,'alpha',0.01)))

fprintf('Altered connections: %.2f%%\n',100*sum(ttest(FCval{4},FCval{3},'dim',2,'alpha',0.01))/size(FCval{1},1))
fprintf('Corr > 0: %.2f%%\n',100*sum(ttest(FCval{4},FCval{3},'dim',2,'alpha',0.01).*(mean(FCval{4},2)-mean(FCval{3},2))>0)./sum(ttest(FCval{4},FCval{3},'dim',2,'alpha',0.01)))


load ~/pmod/matlab/old/WC_connectome_DPZ_AAL
subplot(1,2,2); hold on

bar(Cee); axis([0 5 0 0.02]); 
line([1 1],[Cee(1)-CeeSD(1)/sqrt(nTrials) Cee(1)+CeeSD(1)/sqrt(nTrials)])
line([2 2],[Cee(2)-CeeSD(2)/sqrt(nTrials) Cee(2)+CeeSD(2)/sqrt(nTrials)])
line([3 3],[Cee(3)-CeeSD(3)/sqrt(nTrials) Cee(3)+CeeSD(3)/sqrt(nTrials)])
line([4 4],[Cee(4)-CeeSD(4)/sqrt(nTrials) Cee(4)+CeeSD(4)/sqrt(nTrials)])
axis square
tp_editplots

fprintf('Altered connections: %.2f%%\n',100*sum(ttest(FCval{1},FCval{2},'dim',2,'alpha',0.01))/size(FCval{1},1))
fprintf('Corr > 0: %.2f%%\n',100*sum(ttest(FCval{2},FCval{1},'dim',2,'alpha',0.01).*(mean(FCval{2},2)-mean(FCval{1},2))>0)./sum(ttest(FCval{1},FCval{2},'dim',2,'alpha',0.01)))

fprintf('Altered connections: %.2f%%\n',100*sum(ttest(FCval{3},FCval{4},'dim',2,'alpha',0.01))/size(FCval{1},1))
fprintf('Corr > 0: %.2f%%\n',100*sum(ttest(FCval{4},FCval{3},'dim',2,'alpha',0.01).*(mean(FCval{4},2)-mean(FCval{3},2))>0)./sum(ttest(FCval{4},FCval{3},'dim',2,'alpha',0.01)))

print(gcf,'-dpdf',sprintf('~/pmod/plots/whole_brain_barplots.pdf'))

%% TASK VS REST AUTOCORRELATION

figure; set (gcf,'color','w');

subplot(1,2,1); hold on


bar([mean(mean(Epsilon{1})) mean(mean(Epsilon{3}))]); axis([0 5 0 0.03]); 
line([1 1],[mean(mean(Epsilon{1}))-std(mean(Epsilon{1})) mean(mean(Epsilon{1}))+std(mean(Epsilon{1}))])
line([2 2],[mean(mean(Epsilon{3}))-std(mean(Epsilon{3})) mean(mean(Epsilon{3}))+std(mean(Epsilon{3}))])

tp_editplots
subplot(1,2,2); hold on

% figure; set (gcf,'color','w');

bar([mean(mean(FCval{1})) mean(mean(FCval{3}))]); axis([0 5 0 0.03]); 
line([1 1],[mean(mean(FCval{1}))-std(mean(FCval{1})) mean(mean(FCval{1}))+std(mean(FCval{1}))])
line([2 2],[mean(mean(FCval{3}))-std(mean(FCval{3})) mean(mean(FCval{3}))+std(mean(FCval{3}))])
tp_editplots
print(gcf,'-dpdf',sprintf('~/pmod/plots/whole_brain_barplots_corr.pdf'))

%% FC MATRICES
clear FC_atx FC_dpz

pars = [];
pars.grid = 'medium';
pars.N = 90;


load ~/pmod/matlab/WC_connectome_ATX_AAL
for i = 1 : 4; FC_atx(:,:,i) = tp_match_aal(pars,FC(:,:,i)); end 
load ~/pmod/matlab/WC_connectome_DPZ_AAL
for i = 1 : 4; FC_dpz(:,:,i) = tp_match_aal(pars,FC(:,:,i)); end 

FC1 = tril(FC_atx(:,:,1),-1)+rot90(tril(FC_atx(:,:,3),-1),2);

figure; set(gcf,'color','w'); tp_editplots

h = subplot(2,3,1)
imagesc(FC1,[0 0.1]); 
axis square; colormap(h,plasma); axis off

g = subplot(2,3,2)
FC2 = tril(FC_atx(:,:,2)-FC_atx(:,:,1),-1)+rot90(tril(FC_atx(:,:,4)-FC_atx(:,:,3),-1),2);
imagesc(FC2,[-0.02 0.02]); 
axis square; colormap(g,cmap); axis off

g = subplot(2,3,3)
FC2 = tril(FC_dpz(:,:,2)-FC_dpz(:,:,1),-1)+rot90(tril(FC_dpz(:,:,4)-FC_dpz(:,:,3),-1),2);
imagesc(FC2,[-0.02 0.02]); 
axis square; colormap(g,cmap); axis off

h = subplot(2,3,4)
imagesc(tril(nanmean(cleandat(:,:,:,1,1,6),3),-1)+rot90(tril(nanmean(cleandat(:,:,:,1,2,6),3),-1),2),[0.02 0.06]); 
axis square; colormap(h,plasma); axis off

g = subplot(2,3,5)
FC2 = tril(nanmean(cleandat(:,:,:,2,1,6),3),-1)-tril(nanmean(cleandat(:,:,:,1,1,6),3),-1)+rot90(tril(nanmean(cleandat(:,:,:,2,2,6),3),-1)-tril(nanmean(cleandat(:,:,:,1,2,6),3),-1),2);
imagesc(FC2,[-0.01 0.01]); 
axis square; colormap(g,cmap); axis off

g = subplot(2,3,6)
FC2 = tril(nanmean(cleandat(:,:,:,3,1,7),3),-1)-tril(nanmean(cleandat(:,:,:,1,1,7),3),-1)+rot90(tril(nanmean(cleandat(:,:,:,3,2,7),3),-1)-tril(nanmean(cleandat(:,:,:,1,2,7),3),-1),2);
imagesc(FC2,[-0.01 0.01]); 
axis square; colormap(g,cmap); axis offs



print(gcf,'-dpdf',sprintf('~/pmod/plots/whole_brain_sim_fcmat.pdf'))

%% CORRELATION WITH EMPIRICAL DATA

cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

v=1;
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
mask    = logical(tril(ones(90,90),-1));

fc_emp_atx(:,:,1:2) =  nanmean(cleandat(:,:,:,1:2,1,6),3);
fc_emp_atx(:,:,3:4) =  nanmean(cleandat(:,:,:,1:2,2,6),3);

fc_emp_dpz(:,:,1:2) =  nanmean(cleandat(:,:,:,[1 3],1,6),3);
fc_emp_dpz(:,:,3:4) =  nanmean(cleandat(:,:,:,[1 3],2,6),3);

FC_emp_atx = reshape((fc_emp_atx(repmat(mask,[1 1 4]))),[4005 4]);
FC_emp_dpz = reshape((fc_emp_dpz(repmat(mask,[1 1 4]))),[4005 4]);

FC_sim_atx = reshape((FC_atx(repmat(mask,[1 1 4]))),[4005 4]);
FC_sim_dpz = reshape((FC_dpz(repmat(mask,[1 1 4]))),[4005 4]);

r_rest = corr(FC_emp_atx(:,1),FC_sim_atx(:,1));
r_task = corr(FC_emp_atx(:,2),FC_sim_atx(:,2));

r_atx_rest = corr(FC_emp_atx(:,2)-FC_emp_atx(:,1),FC_sim_atx(:,2)-FC_sim_atx(:,1));
r_atx_task = corr(FC_emp_atx(:,4)-FC_emp_atx(:,3),FC_sim_atx(:,4)-FC_sim_atx(:,3));

r_dpz_rest = corr(FC_emp_dpz(:,2)-FC_emp_atx(:,1),FC_sim_dpz(:,2)-FC_sim_dpz(:,1));
r_dpz_task = corr(FC_emp_dpz(:,4)-FC_emp_atx(:,3),FC_sim_dpz(:,4)-FC_sim_dpz(:,3));



