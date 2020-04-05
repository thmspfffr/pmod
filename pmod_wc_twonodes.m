%% pmod_wc_twonodes

clear

% VERSION 2
% -------------------------------
% v         = 1;
% Gains     = 0:0.05:0.6;
% SIG       = [0.005];
% Ies         = -6:0.05:1;
% Iis         = -7:0.05:-0;
% -------------------------------
% VERSION 2
% -------------------------------
v         = 2;
Gains     = 0:0.05:0.6;
Ies       = -4:0.025:-1;
Iis       = -5:0.025:-2;
tmax      = 6500; % in units of tauE
% -------------------------------
nTrials = 1;

numIes = length(Ies);
N = 4;

dt=0.01;
tspan=0:dt:tmax;
L = length(tspan);
ds = 5;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

tauE = 1;
tauI = 2;
% tau = zeros(2*N,1);
tau([1 3]) = tauE;
tau([2 4]) = tauI;

sigma = 0.0005;

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

% Connectivity:
W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [1 0; 0 0];
W21 = [1 0; 0 0];
W = [W11 W12; W21 W22];

% ----------------
% WAVELETS
% ----------------
[wavelet,f]=tp_mkwavelet(11.3137,0.5,(1/resol));
delta_time = 6/pi./(f(2)-f(1));
delta_time = round(delta_time*1000)/1000;
t_shift    = delta_time;
n_win      = round(delta_time*(1/resol));
n_shift    = round(t_shift*(1/resol));
nseg=floor((L/10-n_win)/n_shift+1);
%%
% fixed params:

for iies = 1:length(Ies)
  for iiis = 1:length(Iis)
    for igain = 1 : length(Gains)
      
      if ~exist(sprintf('~/pmod/proc/twonodes/v%d/',v))
        mkdir(sprintf(['~/pmod/proc/twonodes/' 'v%d'],v))
      end

      outdir = sprintf(['~/pmod/proc/twonodes/v%d/'],v);
      
      fn = sprintf('pmod_wc_twonodes_Ie%d_Ii%d_gain%d_v%d',iies,iiis,igain,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
%       
      fprintf('Processing Ie%d Ii%d Gain%d ...\n',iies,iiis,igain)
      
      Io=zeros(N,1);
      Io(2) = Iis(iiis);
      Io(4) = Iis(iiis);
      Io(1) = Ies(iies);
      Io(3) = Ies(iies);

      uu = -4:.1:4;
      Transf = zeros(length(uu),2,4);
      sdt   = sqrt(dt)*sigma;
      
      % transfer functions:
      % gains are given by 1/aE and 1/aI
      gE = 1+Gains(igain);
      gI = 1+Gains(igain);

      aE  = 1/gE;
      aI  = 1/gI;
      
      Fe = @(x) 1./(1 + exp(-x/aE) );
      Fi = @(x) 1./(1 + exp(-x/aI) );
      F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];
      
      T           = Tds*resol; %% define time of interval
      freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
      freq100     = freqs(freqs<100 & freqs>1);
      pp          = 1:10:length(freq100);
      
      for trial = 1:nTrials
        tic
%         fprintf('Simulating trial %d...\n',trial)
        r = 0.001*rand(N,1);
        R = zeros(Tds,N);
        tt = 0;
        
        % Warm-up:
        for t = 1:5000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau + sdt*randn(N,1);
        end
        tt = 1; clear R Ri
        for t = 1:L
%           t
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau + sdt*randn(N,1);
          if mod(t,ds)==0
            tt=tt+1;
            R(tt,:)  = r([1 3]);
            Ri(tt,:) = r([2 4]);
          end
        end
        
        rE = R;
        tmp = corr(rE);
        out.rc_fr = tmp(1,2);
        
        for j=1:nseg
          dloc2=rE((j-1)*n_shift+1:(j-1)*n_shift+n_win,:)';
          dataf(j,:)=abs(dloc2*wavelet).^2;
        end
        out.rc_wl = corr(dataf);
        out.rc_wl_cov = cov(dataf);
        toc
      end
      
      clear R Ri rE 
      % ------------------------------
      % DETERMINE IF SYSTEM OSCILLATES
      % same as above but without noise
      % ------------------------------
      for trial = 1:nTrials
        tic
%         fprintf('Simulating trial %d...\n',trial)
        r = rand(N,1);
        R = zeros(Tds,N);        
        % Warm-up:
        for t = 1:5000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau;
        end
        tt = 1; clear R Ri
        for t = 1:L
%           t
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau ;
          if mod(t,ds)==0
            tt=tt+1;
            R(tt,:)  = r([1 3]);
            Ri(tt,:) = r([2 4]);
          end
        end
        
        for i=1:2
          [out.osc(trial,:,i) ] = tp_detect_osc(R(:,i));    
        end
        
      end
      
%       for i=1:N
%             
%             % COMPUTE POWER SPECTRUM
%             % ---------------------------
%             f = rE(:,i) - mean(rE(:,i));
%             xdft = fft(f);
%             xdft = xdft(1:floor(Tds/2)+1);
%             pw = (1/(Tds/2)) * abs(xdft).^2;
%             psd = pw(freqs<100 & freqs>1);
%             f = freqs(freqs<100 & freqs>1);
%             fnew = f(1:10:end);
%             psd  = psd(1:10:end);
%             PSD(:,i,trial) = psd';
%             f = fnew;
%             
%       end
%           
%       [~,peak_idx]=max(smooth(mean(PSD(f>6,:),2),20));
%       out.peakfreq = f(peak_idx+find(f<6,1,'last'));
      % ------------------------------
      
      
      save(sprintf([outdir '/%s.mat'],fn),'out')
        %
        % make sure file is saved
        while 1
          try
            load(sprintf([outdir '%s.mat'],fn))
            break
          catch me
            save(sprintf([outdir '%s.mat'],fn),'out')
          end
        end
        
        tp_parallel(fn,outdir,0,0)
    end
  end
end

error('!')

%% PLOT RESULTS

% load(sprintf(['~/pmod/proc/twonodes/v%d/output_v%d.mat'],v,v))

for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for igain = 1 : length(Gains)
      
      outdir = sprintf(['~/pmod/proc/twonodes/v%d/'],v);
      
      fn = sprintf('pmod_wc_twonodes_Ie%d_Ii%d_gain%d_v%d',iies,iiis,igain,v);
    	load(sprintf([outdir '%s.mat'],fn))      
      
      
      outp.fc(iies,iiis,igain)=out.rc_wl(1,2);
      outp.osc(iies,iiis,igain) = mean(squeeze(out.osc));
      outp.fc_fr(iies,iiis,igain)=out.rc_fr;
      outp.cov(iies,iiis,igain) = out.rc_wl_cov(1,2);
    end
  end
end

save(sprintf(['~/pmod/proc/twonodes/v%d/output_v%d.mat'],v),'outp')

%% PLOT

load redblue.mat

for gain = 6

h=figure; set(gcf,'color','w');


subplot(2,3,1);
igain = 1;

par = outp.fc(:,:,igain);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[0 0.5]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(gca,plasma)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots

subplot(2,3,2);
igain = gain;

par = outp.fc(:,:,igain);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[0 0.5]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(gca,plasma)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots

k=subplot(2,3,3);
igain = gain;

par = outp.fc(:,:,igain)-outp.fc(:,:,1);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[-0.5 0.5]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(k,redblue)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots


subplot(2,3,4);
igain = 1;

par = outp.fc_fr(:,:,igain);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[0 0.5]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(gca,plasma)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots


subplot(2,3,5);
igain = gain;

par = outp.fc_fr(:,:,igain);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[0 0.5]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(gca,plasma)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots

subplot(2,3,6);
igain = gain;

par = outp.fc_fr(:,:,gain)-outp.fc_fr(:,:,1);
par(outp.osc(:,:,igain)==1)=nan;

b=imagesc(par,[-0.2 0.2]); axis square
set(b,'AlphaData',~isnan(par))
set(gca,'ydir','normal')
colormap(gca,redblue)

set(h,'Renderer','Painters')
set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  
xlabel('Input to I'); ylabel('Input to E'); tp_editplots


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_twonodes_gain%d_v%d.pdf',igain,v))
end







