%% pmod_wc_singlenode
% -----------------------------------
% COMPUTE LONG-RANGE TEMPORAL CORRELATIONS IN ORDER TO TUNE SINGLE WC NODES
% TO INTRODUCE HETEROGENEITY IN THE LARGE-SCALE VERSION OF THIS MODEL
% -----------------------------------
% Fitted parameters should be introduced back into the model in:
% (1) pmod_WC_wholebrain_rest.m
% (2) pmod_WC_wholebrain_drugs.m
% -----------------------------------

clear

% -------------------------------
% VERSION 1
% -------------------------------
v    = 1;
Ies  = -6:0.1:0;
Iis  = -6:0.1:0;
% -------------------------------

sigma = 0.0005;
N = 2;

wII=4;
wIE=16;
wEI=12;
wEE=12;

g = 1;

tauE = 1;
tauI = 2;
tau = [tauE;tauI];

% Iis = -4;

numIis = length(Iis);
numIes = length(Ies);
Io=zeros(N,1);

dt=0.01;
tmax = 10000; % in units of tauE
tspan=0:dt:tmax;
L = length(tspan);

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;


% Connectivity:s
W11 = [wEE -wEI; wIE -wII];
W   = W11;

nTrials = 10;

uu = -4:.1:4;
Transf = zeros(length(uu),2);

%%


for k = 1:numIis
  for l = 1:numIes
    
    if ~exist(['~/pmod/proc/' sprintf('pmod_wc_singlenode_ii%d_ie%d_v%d.txt',k,l,v)])
      system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_singlenode_ii%d_ie%d_v%d.txt',k,l,v)]);
    else
      continue
    end
     
    for trial = 1:nTrials
      fprintf('Starting computation Ie%d Ii%d trial %d ...\n',k,l,trial);
      
      sdt = sqrt(dt)*sigma;
      
      Io(1) = Ies(l);
      Io(2) = Iis(k);
      
      % SAVE PARAMETERS
      out.Ie = Io(1);
      out.Ii = Io(2);
      
      % transfer functions:
      % gains are given by 1/aE and 1/aI
      aE = 1;
      aI = 1;
      Fe = @(x) 1./(1 + exp(-x/aE) );
      Fi = @(x) 1./(1 + exp(-x/aI) );
      F = @(x) [feval(Fe,x(1));feval(Fi,x(2))];
      
      r = rand(N,1);
      R = zeros(Tds,N);
      tt = 0;
      
      % Warm-up:
      for t = 1:2000
        u = W*r + Io;
        K = feval(F,u);
        r = r + dt*(-r./tau + K) + sdt*randn(N,1);
      end
      
      for t = 1:L
        u = W*r + Io;
        K = feval(F,u);
        r = r + dt*(-r./tau + K) + sdt*randn(N,1);
        if mod(t,ds)==0
          tt=tt+1;
          R(tt,:)=r;
        end
      end
      % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
      % On E time course
      % ---------------------
      tmp                   = tp_dfa(R,[5 75],1/dt,0.5,20);
      out.dfa(:,trial)      = tmp.exp;
      out.fluc(:,1,trial)  	= tmp.y{1};
      out.fluc(:,2,trial)  	= tmp.y{2};
      out.win(:,trial)      = tmp.win;
      % ---------------------
      
    end
    
    fprintf('Saving...\n');
    save(sprintf('~/pmod/proc/pmod_wc_singlenode_ii%d_ie%d_v%d.mat',k,l,v),'out')
    
  end
end

error('!')

%%
for k = 1 : length(Ies)
  for l = 1 : length(Iis)
    
    load(sprintf('~/pmod/proc/pmod_wc_singlenode_ii%d_ie%d_v%d.mat',k,l,v))
    
    dfa_all(k,l,:) = out.dfa(1,:);
  end
end



