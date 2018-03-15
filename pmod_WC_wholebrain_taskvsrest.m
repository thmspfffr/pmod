
% Stochastic simulation of 2*N WC nodes subject to drug at t=0
%-------------------------------------------------------------------------

clear

% LOAD CLEANED DATA AND COMPUTE MEAN
v = 1;
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
fc_emp_atx(:,:,1) = squeeze(nanmean(cleandat(:,:,:,2,1,6),3))-squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_emp_atx(:,:,2) = squeeze(nanmean(cleandat(:,:,:,2,2,6),3))-squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
fc_emp_dpz(:,:,1) = squeeze(nanmean(cleandat(:,:,:,3,1,7),3))-squeeze(nanmean(cleandat(:,:,:,1,1,7),3));
fc_emp_dpz(:,:,2) = squeeze(nanmean(cleandat(:,:,:,3,2,7),3))-squeeze(nanmean(cleandat(:,:,:,1,2,7),3));

% load connectome
load EC %Matt_EC
C = EC;
C = C/max(C(C>0));
N = size(C,1);


%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------

% fixed params:
nTrials = 10;
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
W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];


% Control params.
%--------------------
Ie = -2.85;
Ii = -3.50; %-1;

Io=zeros(2*N,1);
Io(1:N) = Ie;
Io(N+1:2*N) = Ii;


%
% dI_drug = zeros(2*N,1);
% dI_drug(1:N) = dIe_drug;
% dI_drug(N+1:2*N) = dIi_drug;

Rate = zeros(2,1);
RateSD = zeros(2,1);

Cee = zeros(2,1);
CeeSD = zeros(2,1);

FC = zeros(N,N,2);
FCval = cell(1,2);
Epsilon = cell(1,2);

KOPsd   = zeros(nTrials,2);
KOPmean = zeros(nTrials,2);

isub = find( triu(ones(N)) - eye(N) );

dIE = 0.5:0.1:1.5;
dII = 0.7:0.1:1.7;

for idIE = 1 : length(dIE)
  for idII = 1 : length(dII)
    
    if ~exist(sprintf(['~/pmod/proc/' 'pmod_WC_wholebrain_taskvsrest_idIE%d_idII%d_processing.txt'],idIE,idII))
      system(['touch ' '~/pmod/proc/' sprintf('pmod_WC_wholebrain_taskvsrest_idIE%d_idII%d_processing.txt',idIE,idII)]);
    else
      continue
    end
    % input due to task
    %--------------------
    %   dIe = 1.0;
    dIe = dIE(idIE);
    %   dIi = 1.2; %1.3;
    dIi = dII(idII);
    %--------------------
    
    dI=zeros(2*N,1);
    dI(1:N) = dIe;
    dI(N+1:2*N) = dIi;

    % SIMULATE RESTING STATE
    %--------------------------------------------------------------------------
    display('REST ...')
    %--------------------------------------------------------------------------
    
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
    
    FCval{1}    = zeros(length(isub),nTrials);
    Epsilon{1}  = zeros(N,nTrials);
    
    T       = Tds*resol; %% define time of interval
    freqs   = (0:Tds/2)/T; %% find the corresponding frequency in Hz
    nfreqs  = length(freqs);
    freq100 = freqs(freqs<100 & freqs>1);
    pp      = 1:10:length(freq100);
    PSD     = zeros(length(pp),N,nTrials);
    
    PSDstruct(1).frequencies = freq100(1:10:end)';
    
    
    for tr=1:nTrials
      display(sprintf('trial: %g',tr))
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
      fc              = rc(isub);
      Cee(1)          = Cee(1) + mean(fc)/nTrials;
      CeeSD(1)        = CeeSD(1) + var(fc)/nTrials;
      FCval{1}(:,tr)  = fc;
      
      rEo       = mean(mean(rE));
      rEsd      = var(mean(rE));
      Rate(1)   = rEo-rEo;
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
        PSD(:,i,tr) = psd';
      end
      
    end
    
    CeeSD(1)  = sqrt(CeeSD(1));
    RateSD(1) = sqrt(RateSD(1));
    
    PSDstruct(1).PSD = PSD;
    
    
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
      tmp = tp_dfa(R,[3 50],ds,0.5,15);
      dfa(tr,:,2) = tmp.exp;
      
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
        PSD(:,i,tr) = psd';
      end
      
    end
    
    CeeSD(2)  = sqrt(CeeSD(2));
    RateSD(2) = sqrt(RateSD(2));
    
    PSDstruct(2).PSD = PSD;
    
    out = struct('dfa',dfa,'FC',FC,'Epsilon',Epsilon,'dIe',dIe,'dIi',dIi);
    
    save(sprintf('~/pmod/proc/pmod_WC_wholebrain_taskvsrest_idIE%d_idII%d.mat',idIE,idII),'out')
    
  end
end

%--------------------------------------------------------------------------
% Save:
%
%   if Drug ==1
%     save('WC_connectome_ATX_AAL','Rate','RateSD','Cee','CeeSD','FC','FCval',...
%       'Ie','Ii','dIe','dIi','dIe_drug','dIi_drug','Gain_E','Gain_I','Epsilon','PSDstruct','KOPsd','KOPmean')
%   else
%     save('WC_connectome_DPZ_AAL','Rate','RateSD','Cee','CeeSD','FC','FCval',...
%       'Ie','Ii','dIe','dIi','dIe_drug','dIi_drug','Gain_E','Gain_I','Epsilon','PSDstruct','KOPsd','KOPmean')
%   end

error('!')
%%
FC = FFFF;
fc_emp = squeeze(nanmean(cleandat(:,:,:,1,1:2,6),3));

FC(FC==1)=nan;

mask = logical(repmat(tril(ones(size(FC(:,:,1))),-1),[1 1 2]));


m = mean(FC(mask));
s = std(FC(mask));
m_emp = mean(fc_emp(mask));
s_emp = std(fc_emp(mask));

for i = 1 : 2
  
  pars = [];
  pars.N = 90;
  pars.grid = 'medium';
  FC(:,:,i) = tp_match_aal(pars,FC(:,:,i));
  
  FC(:,:,i) = (FC(:,:,i) - m)./s;
  %
  fc_emp(:,:,i) = (fc_emp(:,:,i) - m_emp)./s_emp;
  %
end

figure;

subplot(2,2,1)
imagesc(fc_emp(:,:,1),[-2 3]); colormap(jet); axis square
subplot(2,2,2)
imagesc(FC(:,:,1),[-2 3]); colormap(jet); axis square
subplot(2,2,3)
imagesc(fc_emp(:,:,2),[-2 3]); colormap(jet); axis square
subplot(2,2,4)
imagesc(FC(:,:,2),[-2 3]); colormap(jet); axis square

mask = logical(tril(ones(size(FC(:,:,1))),-1));

d1 = fc_emp(:,:,2)-fc_emp(:,:,1);
d2 = FC(:,:,2)-FC(:,:,1);





