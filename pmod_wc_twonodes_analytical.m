%% pmod_wc_twonodes_analytical

clear

% ----------------------------
% VERSION 1
Ies = -10:.1:5; %-4:.1:-1;
Iis = -10:.1:5; %-4:.1:-1;
tmax = 20; %5000;
Gains = 0:0.05:0.5;
% ----------------------------

N = 4;
% Connectivity:

g     = 1;
tauE  = 1;
tauI  = 2; %1
tau   = [tauE;tauI;tauE;tauI];
Io    = zeros(N,1);
dt    = 0.01;
tspan = 0:dt:tmax;
L     = length(tspan);
ds    = 10;
Tds   = length(0:ds*dt:tmax)-1;
sigma = 0.0005;
Qn    = (sigma*dt)^2*eye(N);

% Connectivity:
% ----------------
wII=4;
wIE=16;
wEI=12;
wEE=12;
W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [g 0; 0 0];
W21 = [g 0; 0 0];
W = [W11 W12; W21 W22];

%%
for iies = 1:length(Ies)
  for iiis = 1:length(Iis)
    for igain = 1 : length(Gains)
      
      if ~exist(sprintf(['~/pmod/proc/' 'pmod_wc_twonodes_analytical_Ie%d_Ii%d_igain%d_v%d_processing.txt'],iies,iiis,igain,v))
        system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_twonodes_analytical_Ie%d_Ii%d_igain%d_v%d.mat',iies,iiis,igain,v)]);
      else
        continue
      end
      
      Ii = Iis(iies);
      Io(2) = Ii;
      Io(4) = Ii;
      
      Ie = Ies(iiis);
      Io(1) = Ie;
      Io(3) = Ie;
      % transfer functions:
      % gains are given by 1/aE and 1/aI
      
      aE = 1/(1+Gains(igain));
      aI = 1/(1+Gains(igain));
      
      Fe = @(x) 1./(1 + exp(-x/aE) );
      Fi = @(x) 1./(1 + exp(-x/aI) );
      F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];
      
      r = rand(N,1);
      R = zeros(Tds,N);
      tt = 0;
      
      % get fixed points:
      for t = 1:300000
        u = W*r + Io;
        K = feval(F,u);
        r = r + dt*(-r + K)./tau;
      end
      
      for t = 1:L
        u = W*r + Io;
        K = feval(F,u);
        r = r + dt*(-r+ K)./tau;
        if mod(t,ds)==0
          tt=tt+1;
          R(tt,:)=r;
        end
      end
      
      RstatE = R(:,1);
      maxR = max(RstatE);
      minR = min(RstatE);
      
      
      if maxR-minR<10e-10
        Cov = zeros(N);
        for t = 1:10000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau;
          
          re = r(1);
          ri = r(2);
          rE = r(3);
          rI = r(4);
          % Jacobian:
          Aee = -1 + wEE*1/aE*re*(1-re);
          Aei = -wEI*1/aE*re*(1-re);
          Aie = wIE*1/aI*ri*(1-ri)*(tauE/tauI);
          Aii = -1*(tauE/tauI)- wII*1/aI*ri*(1-ri)*(tauE/tauI);
          
          AeE = g*1/aE*re*(1-re);
          AeI = 0;
          AiE = 0;
          AiI = 0;
          
          AEe = g*1/aE*rE*(1-rE);
          AIe = 0;
          AEi = 0;
          AIi = 0;
          
          AEE = -1/tauE + wEE*1/aE*rE*(1-rE);
          AEI = -wEI*1/aE*rE*(1-rE);
          AIE = wIE*1/aI*rI*(1-rI)*(tauE/tauI);
          AII = -1*(tauE/tauI) - wII*1/aI*rI*(1-rI)*(tauE/tauI);
          
          
          Jmat = [Aee Aei AeE AeI; ...
            Aie Aii AiE AiI; ...
            AEe AEi AEE AEI; ...
            AIe AIi AIE AII];
          
          delta = (Jmat*Cov + Cov*Jmat' + Qn)*dt;
          Cov = Cov + delta;
          
          tol = norm(delta);
          if tol<10e-18
            break
          end
        end
        
        Corr=zeros(N);
        for m=1:N
          for n=1:N
            Corr(m,n)=Cov(m,n)/( sqrt(Cov(m,m))*sqrt(Cov(n,n)) );
          end
        end
        
        CeE(k,l,Q) = Corr(1,3);
        
      else
        CeE(k,l,Q) = NaN;
      end
      
    end
  end
end

save(sprintf(['~/pmod/proc/pmod_wc_twonodes_analytical_Ie%d_Ii%d_igain%d_v%d'],iies,iiis,igain,v),'out')

%%



