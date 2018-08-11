%% pmod_wc_twonodes.m

clear

% VERSION 2
% -------------------------------
v         = 1;
GAIN      = 0;
SIG       = [0.005];
% Iis   = [-4.2 -4 -3.8];
Ies       = -8:0.1:8;
Iis       = -8:0.1:9;
gs        = 0:0.2:5;
wins = [3 50];
% -------------------------------
nTrials = 1;

numIes = length(Ies);
N = 4;

tmax        = 65000; % in units of tauE
dt=0.01;
tspan=0:dt:tmax;
L = length(tspan);
ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

tauE = 1;
tauI = 2;
% tau = zeros(2*N,1);
tau(1:2) = tauE;
tau(3:4) = tauI;

sigma = 0.0005;

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;


flp = 8;           % lowpass frequency of filter
fhi = 13;

delt  = resol;            % sampling interval
k     = 4;                  % 2nd order butterworth filter
fnq   = 1/(2*delt);       % Nyquist frequency
Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt] = butter(k,Wn);

%%
% fixed params:

for iies = 1 : length(Ies)
  for iiis = 1 : length(Iis)
    for ig = 1 : length(gs)
      
      if ~exist(['~/pmod/proc/' sprintf('pmod_wc_excbif_iis%d_ies%d_g%d_v%d_processing.txt',iiis,iies,ig,v)])
        system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_excbif_iis%d_ies%d_g%d_v%d_processing.txt',iiis,iies,ig,v)]);
      else
        continue
      end
      
      
      % Connectivity:
      W11 = [wEE -wEI; wIE -wII];
      W22 = W11;
      W12 = [gs(ig) 0; 0 0];
      W21 = [gs(ig) 0; 0 0];
      W = [W11 W12; W21 W22];
      
      Io=zeros(N,1);
      Io(2) = Iis(iiis);
      Io(4) = Iis(iiis);
      Io(1) = Ies(iies);
      Io(3) = Ies(iies);
      
      Rc     = zeros(length(Ies),4);
      RcAmpl = zeros(length(Ies),4);
      
      
      uu = -4:.1:4;
      Transf = zeros(length(uu),2,4);
      
      Gain  = 0;
      sdt   = sqrt(dt)*sigma;
      
      % transfer functions:
      % gains are given by 1/aE and 1/aI
      gE = 1+Gain;
      gI = 1+Gain;
      
      aE = 1/gE;
      aI = 1/gI;
      
      Fe = @(x) 1./(1 + exp(-x/aE) );
      Fi = @(x) 1./(1 + exp(-x/aI) );
      F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];
      
      T           = Tds*resol; %% define time of interval
      freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
      freq100     = freqs(freqs<100 & freqs>1);
      pp          = 1:10:length(freq100);
      
      for trial = 1:nTrials
        fprintf('Simulating trial %d...\n',trial)
        r = rand(N,1);
        R = zeros(Tds,N);
        tt = 0;
        
        % Warm-up:
        for t = 1:5000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau' + sdt*randn(N,1);
        end
        
        for t = 1:L
          
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau' + sdt*randn(N,1);
          if mod(t,ds)==0
            tt
            tt=tt+1;
            R(tt,:)=r;
          end
        end
        
        rE = R(:,[1 3]);
        Ampl(:,1)=abs(R(:,1) + 1i*R(:,2));
        Ampl(:,2)=abs(R(:,3) + 1i*R(:,4));
        
        outp.fc(trial) = corr(rE(:,1),rE(:,2));
        outp.fc_ampl(trial) = corr(Ampl(:,1),Ampl(:,2));
        
        tmp           = tp_dfa(rE,[3 50],1/resol,0.5,15);
        out.dfa(:,trial) = tmp.exp;
        
        for i=1:N
          fprintf('Autocorr reg %d ...\n',i)
          outp.osc(trial,:,i) = tp_detect_osc(rE(:,i));
          %autocorr
          lags = 1:round(2*(1/resol));
          if license('checkout','econometrics_toolbox')
            acorr = autocorr(rE(:,i),'NumLags',round(2*(1/resol)));
          else
            acorr = acf(rE(:,i),round(2*(1/resol)));
          end
          
          % get exp decay
          outp.lambda(i,trial) = tp_fitexpdecay(acorr(1:size(lags,2)),lags,0.01);
          
          ii = find(acorr<.2,1,'first');
          
          if isempty(ii)
            outp.lags(i,trial) = nan;
          else
            outp.lags(i,trial) = lags(ii);
          end
          
          % COMPUTE POWER SPECTRUM
          % ---------------------------
          f = rE(:,i) - mean(rE(:,i));
          xdft = fft(f);
          xdft = xdft(1:floor(Tds/2)+1);
          pw = (1/(Tds/2)) * abs(xdft).^2;
          psd = pw(freqs<100 & freqs>1);
          f = freqs(freqs<100 & freqs>1);
          fnew = f(1:10:end);
          psd  = psd(1:10:end);
          outp.PSD(:,i,trial) = psd';
          outp.f = fnew;
          
        end
        
        [~,peak_idx]=max(smooth(mean(out.PSD(outp.f>3,:),2),20));
        outp.peakfreq = outp.f(peak_idx+find(outp.f<4,1,'last'));
        
        
        % COMPUTE DFA, EXTRACT HURST
        % ---------------------------
        env = abs(hilbert(filtfilt(bfilt,afilt,rE)));
        tmp                 = tp_dfa(env,wins,1/resol,0.5,15);
        outp.dfa_env(trial,:)   = tmp.exp;
        
        outp.fc_env(trial) =corr(env(:,1),env(:,2));
        
      end
      
      save(sprintf('~/pmod/proc/pmod_wc_excbif_iis%d_ies%d_g%d_v%d.mat',iiis,iies,ig,v),'outp')
      
    end
  end
end

