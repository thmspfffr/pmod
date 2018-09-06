%% pmod_wc_wholebrain_final
% Stochastic simulation of 2*N WC nodes during "rest"
%-------------------------------------------------------------------------

clear

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.
%-------------------------------------------------------------------------
% VERSION 2: After meeting with tobi, 24-08-2018: even more fine gained
% %-------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0.6;
Gains       = 0;
nTrials     = 1;
tmax        = 6500; % in units of tauE
wins = [2 20]; 
%-------------------------------------------------------------------------
% VERSION 3: After meeting with tobi, 24-08-2018
%-------------------------------------------------------------------------
% v           = 3;
% Ies         = -4:0.01:-1;
% Iis         = -5:0.01:-1;
% Gg          = 0.6;
% Gains       = 0;
% nTrials     = 1;
% tmax        = 6500; % in units of tauE
% wins = [2 20]; 
%-------------------------------------------------------------------------
% VERSION 4: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v           = 4;
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-1;
% Gg          = 0:0.2:1;
% Gains       = 0:0.05:0.2;
% nTrials     = 3;
% tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------
% VERSION 5: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v           = 5;
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-1;
% Gg          = 0:0.2:1;
% Gains       = 0:0.05:0.2;
% nTrials     = 1;
% tmax        = 65000; % in units of tauE
%-------------------------------------------------------------------------

% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));
N = size(C,1);

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/pconn/matlab
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
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
tspan=0:dt:tmax;
L = length(tspan);
clear tspan

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

sigma = 0.0005;
%Qn = (sigma*dt)^2*eye(2*N);

isub = find( triu(ones(N)) - eye(N) );
%%
for iies = 1: length(Ies)
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
%         
        if ~exist(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
          system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt',iies,iiis,iG,igain,v)]);
        else
          continue
        end
        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        
        out.FC = zeros(N,N,1);
        out.Cee = zeros(1,1);
        out.CeeSD = zeros(2,1);
        out.FC_env = zeros(N,N,1);
        out.Cee_env = zeros(1,1);
        out.CeeSD_env = zeros(2,1);
        
        out.KOPsd   = zeros(nTrials,1);
        out.KOPmean = zeros(nTrials,1);
        
        % Control params.
        %--------------------
        out.Ie = Ies(iies);
        out.Ii = Iis(iiis);
        out.Gain = Gains(igain);
        
        % Working point:
        Io=zeros(2*N,1);
        Io(1:N) = out.Ie;
        Io(N+1:2*N) = out.Ii;
        
        % transfer function:
        gE = 1 + Gains(igain);
        gI = 1 + Gains(igain);
        aE  = 1/gE;
        aI  = 1/gI;
        Fe  = @(x) 1./(1 + exp(-x/aE) );
        Fi  = @(x) 1./(1 + exp(-x/aI) );
        F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
        
        T           = Tds*resol; %% define time of interval
        freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
        freq100     = freqs(freqs<100 & freqs>1);
        pp          = 1:10:length(freq100);
        PSD     = zeros(length(pp),N,nTrials);  
        out.PSDstruct(1).frequencies = freq100(1:10:end)';
        
        % RUN SIMULATION
        % ---------------------

        for tr=1:nTrials
          fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
          r   = 0.001*rand(2*N,1);
          R   = zeros(Tds,N);
          Ri  = zeros(Tds,N);
          tt  = 0;
          % transient:
          for t = 1:25000
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(2*N,1);
          end
          % simulation
          for t = 1:L
            %           100*t/L
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(2*N,1);
            if mod(t,ds)==0
              tt=tt+1;
              R(tt,:)  = r(1:N);
              Ri(tt,:)  = r(N+1:end);
            end
          end
          
          rE = R;
          rI = Ri;  
         	z  = rE + 1i*rI;

          clear R Ri rI 
 
          % KURAMOTO PARAMETERS
          % ---------------------
          ku                = sum(z,2)/N;   
          KOP           = abs(ku);
          out.KOPsd(tr,1)   = std(KOP);
          out.KOPmean(tr,1) = mean(KOP);
          
          clear ku KOP z

          % ---------------------
          % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
          % On E time course
          % ---------------------
%           tmp           = tp_dfa(rE,wins,1/resol,0.5,15);
%           out.dfa(tr,:) = single(tmp.exp);
 
          % FC matrix
          % ---------------------
          rc       	= corrcoef(rE);
          out.FC  	= single(out.FC) + single(rc/nTrials);
          fc      	= rc(isub);
%           out.Cee  	= out.Cee + mean(fc)/nTrials;
%           out.CeeSD	= out.CeeSD + var(fc)/nTrials;
                
          for i=1:N
%             fprintf('Autocorr reg %d ...\n',i)
% %             out.osc(tr,:,i) = tp_detect_osc(rE(:,i));
%             %autocorr
%             lags = 1:round(2*(1/resol));
%             if license('checkout','econometrics_toolbox')
%               acorr = autocorr(rE(:,i),'NumLags',round(2*(1/resol)));
%             else
%               acorr = acf(rE(:,i),round(2*(1/resol)));
%             end
            
%             % get exp decay 
%             out.lambda(i,tr) = tp_fitexpdecay(acorr(1:size(lags,2)),lags,0.01);
%             
%             ii = find(acorr<.2,1,'first');
%             
%             if isempty(ii)
%               out.lags(i,tr) = nan;
%             else
%               out.lags(i,tr) = lags(ii);           
%             end    
%   
%             % COMPUTE POWER SPECTRUM
%             % ---------------------------
            f = rE(:,i) - mean(rE(:,i));
            xdft = fft(f);
            xdft = xdft(1:floor(Tds/2)+1);
            pw = (1/(Tds/2)) * abs(xdft).^2;
            psd = pw(freqs<100 & freqs>1);
            f = freqs(freqs<100 & freqs>1);
            fnew = f(1:10:end);
            psd  = psd(1:10:end);
            PSD(:,i,tr) = psd';
            f = fnew;
            
%              % POWER SPECTRUM FIT
%             idx= find(log10(out.f)>1.5,1,'first');
%             X = [ones(1,length(out.f(idx:end)))' log10(out.f(idx:end))'];
%             Y = log10(PSD(idx:end,i,tr));
%             tmp = X\Y;   
%             out.psslope(i,tr)= tmp(2);
%         
% 
          end
          
          % EXTRACT PEAK FREQ
          % ---------------------------
          [~,peak_idx]=max(smooth(mean(PSD(f>3,:),2),20));
          out.peakfreq = f(peak_idx+find(f<4,1,'last')); 
          
          flp = 9;           % lowpass frequency of filter
          fhi = 13;
          
          delt  = resol;            % sampling interval
          k     = 4;                  % 2nd order butterworth filter
          fnq   = 1/(2*delt);       % Nyquist frequency
          Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
          [bfilt,afilt] = butter(k,Wn);
          
          % COMPUTE DFA, EXTRACT HURST
          % ---------------------------
          env = abs(hilbert(filtfilt(bfilt,afilt,rE)));
          
          flp = 12;           % lowpass frequency of filter
          fhi = 20;
          
          delt  = resol;            % sampling interval
          k     = 4;                  % 2nd order butterworth filter
          fnq   = 1/(2*delt);       % Nyquist frequency
          Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
          [bfilt,afilt] = butter(k,Wn);
          env_beta = abs(hilbert(filtfilt(bfilt,afilt,rE)));
       
          % COMPUTE CORRELATIONS BASED ON ENV
          % ---------------------------
          % Alpha range
          % ---------------------------
          rc                  = corrcoef(env);         
          out.FC_env          = single(out.FC_env) + single(rc/nTrials);
          fc_env              = rc(isub);
          out.Cee_env         = out.Cee_env + mean(fc)/nTrials;
          out.CeeSD_env       = out.CeeSD_env + var(fc)/nTrials;
%           tmp                 = tp_dfa(env,wins,1/resol,0.5,15);
%           out.dfa_env         = tmp.exp;
          % ---------------------------
          clear rc fc_env
          % ---------------------------
          % Beta filtered
          % ---------------------------
          rc                  = corrcoef(env_beta);         
          out.FC_env_beta     = single(out.FC_env) + single(rc/nTrials);
          fc_env              = rc(isub);
          out.Cee_env_beta    = out.Cee_env + mean(fc)/nTrials;
          out.CeeSD_env_beta  = out.CeeSD_env + var(fc)/nTrials;
          % ---------------------------
%           
%           for i = 1 : N
%             fprintf('Autocorr reg %d ...\n',i)
%             if license('checkout','econometrics_toolbox')
%               acorr_env = autocorr(env(:,i),'NumLags',round(2*(1/resol)));
%             else
%               acorr_env = acf(env(:,i),round(2*(1/resol)));
%             end
%             
%             jj = find(acorr_env<.2,1,'first');
% 
%             if isempty(jj) || jj > length(lags)
%               out.lags_env(i,tr) = nan;
%             else
%               out.lags_env(i,tr) = lags(jj);
%             end
%             
%             %fit exp decay
%            	out.lambda_env(i,tr) = tp_fitexpdecay(acorr_env(1:size(lags,2)),lags,0.01);
            
            % COMPUTE POWER SPECTRUM
            % ---------------------------
%             f = env(:,i) - mean(env(:,i));
%             xdft = fft(f);
%             xdft = xdft(1:floor(Tds/2)+1);
%             pw = (1/(Tds/2)) * abs(xdft).^2;
%             psd = pw(freqs<100 & freqs>1);
%             f = freqs(freqs<100 & freqs>1);
%             fnew = f(1:10:end);
%             psd  = psd(1:10:end);
%             PSD_env(:,i,tr) = psd';
%             out.f_env = fnew;
%             
%              % POWER SPECTRUM FIT
%             idx= find(log10(out.f_env)>0.5,1,'first');
%             X = [ones(1,length(out.f_env(idx:end)))' log10(out.f_env(idx:end))'];
%             Y = log10(PSD_env(idx:end,i,tr));
%             tmp = X\Y;   
%             out.psslope_env(i,tr)= tmp(2);
            
%           end
          
          

          clear rE rI env
          toc
        end
        
%         out.lambda_env = single(mean(out.lambda_env,2));
%        	out.lags_env = single(mean(out.lags_env,2));
%         out.lambda = single(mean(out.lambda,2));
%        	out.lags = single(mean(out.lags,2));

        
        save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v),'out')
        while 1
          try 
            load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v))
            break
          catch me
            save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v),'out')
          end
        end
        
      end
    end
  end
end
error('!')
%%
v = 2
for iies = 1: length(Ies)
  iies
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
%             
        try 
        load(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,v))
          catch me
          delete(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
          warning('!')
        end
        
      end
    end
  end
end
