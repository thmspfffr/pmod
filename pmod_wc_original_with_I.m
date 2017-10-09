

clear


% VERSION 1:
% -------------------------------------
v = 1
SIGMA = 0.002;
GAIN = .20;
Ies = -4:.2:0;
Iis = -4;
% -------------------------------------

N = 2;

wII=4;
wIE=16;
wEI=10;
wEE=12;

tauE = 1;
tauI = 1;
tau = [tauE;tauI];

Ie = -1;
Io=zeros(N,1);
Io(1) = Ie;


dt=0.01;
tmax = 300;
tspan=0:dt:tmax;
L = length(tspan);

ds = 2;
Tds = length(0:ds*dt:tmax)-1;

%%
  for igain = 1 : length(GAIN)
    for iinh = 1 : length(Iis)
      for inoise = 1 : length(SIGMA)
        for iasym = 1 : length(ASYM)


          if ~exist(['~/pmod/proc/' sprintf('pmod_wc_orig_gain%d_noise%d_inh%d_asym%d_v%d_proc.txt',igain,inoise,iinh,iasym,v)])
            system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_orig_gain%d_noise%d_inh%d_asym%d_v%d_proc.txt',igain,inoise,iinh,iasym,v)]);
          else
            continue
          end

% transfer functions:
aE = 1;
Fe = @(x) 1./(1 + exp(-x/aE) );

aI = 1;
Fi = @(x) 1./(1 + exp(-x/aI) );

F = @(x) [feval(Fe,x(1));feval(Fi,x(2))];

% inverse functions:
FinvE = @(x) log( aE*x./(1-aE*x) );

FinvI = @(x) log( aI*x./(1-aI*x) );


Iis = -15:.2:0;

numIis = length(Iis);

Re  = zeros(numIis,2);
Tr  = zeros(numIis,1);
De  = zeros(numIis,1);
Eig = zeros(numIis,N);

for k = 1:numIis
    
    Ii = Iis(k);
    Io(2) = Ii;
        
        W = [wEE -wEI; wIE -wII];
        
        % stability analysis:
        
        % fixed points:
        S = @(x) -x./tau + feval(F,W*x+Io);
        x0 = .5*ones(N,1);
        fp = fsolve(S,x0);
        re = fp(1);
        ri = fp(2);
        
        % Jacobian:
        Aee = -1/tauE + wEE*re*(1-re);
        Aei = -wEI*re*(1-re);
        Aie = wIE*ri*(1-ri);
        Aii = -1/tauI - wII*ri*(1-ri);
        
        A = [Aee Aei; Aie Aii];
        
        trA = trace(A);
        detA = det(A);
        
        d1 = 1/2*(trA + sqrt( trA^2 - 4*detA ));
        d2 = 1/2*(trA - sqrt( trA^2 - 4*detA ));
        
        Tr(k) = trA;
        De(k) = detA;
        Eig(k,1) = d1;
        Eig(k,2) = d2;
           
        
        r = [.6; .2];
        R = zeros(Tds,N);
        tt = 0;
        
        for t = 1:L
            
            u = W*r + Io;           
            K = feval(F,u);           
            r = r + dt*(-r./tau + K);
            
            if mod(t,ds)==0
                tt=tt+1;
                R(tt,:)=r;
            end
            
        end
        
        RstatE = R(end-500:end,1);
        Re(k,1) = max(RstatE);
        Re(k,2) = min(RstatE);
        

end


%%
figure
xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.1 .2 .22 .69])
plot(Iis,Re,'k')
hold on
Ibif = Iis(abs(Re(:,1)-Re(:,2))>0.0001); %bifurcations
plot(Ibif(1)*[1 1],[0 1],'k:')
plot(Ibif(end)*[1 1],[0 1],'k:')
xlabel('Control param. \itI_{I}','fontname','times')
ylabel('\itr_E    ','fontname','times','rotation',0)
text(-.2,1.1,'A','units','normalized','fontsize',14)

axes('position',[.42 .2 .22 .69])
plot(rE,nullE(:,1),'k')
hold on
plot(nullI(:,1),rI,'color',[.5 .5 .5])
plot(trayectory(:,1),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{I}\rm = -4','fontname','times')
text(-.2,1.1,'B','units','normalized','fontsize',14)

axes('position',[.74 .2 .22 .69])
plot(rE,nullE(:,2),'k')
hold on
plot(nullI(:,2),rI,'color',[.5 .5 .5])
plot(trayectory(:,2),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{I}\rm = -6','fontname','times')
text(-.2,1.1,'C','units','normalized','fontsize',14)

% export_fig Fig_WC_Bifurcations.tiff -r300







