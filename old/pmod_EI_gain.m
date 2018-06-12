%% BIFURCATION ANALYSIS
% Computes bifurcation diagrams

clear 

% fixed params:

N = 2;

wII=4;
wIE=16;
wEI=10;
wEE=10.5;

tauE = 1;
tauI = 2;
tau = [tauE;tauI];

Ies = -0.01:0.01:0.01;
Iis = -0.01:0.01:0.01;
Io=zeros(N,1);


dt=0.01;
tmax = 15;
tspan=0:dt:tmax;
L = length(tspan);

ds = 2;
Tds = length(0:ds*dt:tmax)-1;

% transfer functions:
aE = 1;
Fe = @(x) 1./(1 + exp(-x/aE) );

aI = 1;
Fi = @(x) 1./(1 + exp(-x/aI) );

F = @(x) [feval(Fe,x(1));feval(Fi,x(2))];

% inverse functions:
FinvE = @(x) log( aE*x./(1-aE*x) );

FinvI = @(x) log( aI*x./(1-aI*x) );


numIes = length(Ies);
numIis = length(Iis);

all_inp = 0:10:50;

Re  = zeros(numIes,2);
Tr  = zeros(numIes,1);
De  = zeros(numIes,1);
Eig = zeros(numIes,N);

for k = 1:numIes
  k
  for j = 1 :numIis
    for iinp = 1 : length(all_inp)
    fprintf('#%d ..\n',k)
    Ii=-4;
    Ie = 0;
    
        Io(1) = Ie+Ies(k);

        Io(2) = Ii+Iis(j);
        
        W = [wEE -wEI; wIE -wII];
        
        % stability analysis:
        
        % fixed points:
%         S = @(x) -x./tau + feval(F,W*x+Io);
%         x0 = .5*ones(N,1);
%         fp = fsolve(S,x0);
%         re = fp(1);
%         ri = fp(2);
        
        % Jacobian:
%         Aee = -1/tauE + wEE*re*(1-re);
%         Aei = -wEI*re*(1-re);
%         Aie = wIE*ri*(1-ri);
%         Aii = -1/tauI - wII*ri*(1-ri);
%         
%         A = [Aee Aei; Aie Aii];
%         
%         trA = trace(A);
%         detA = det(A);
%         
%         d1 = 1/2*(trA + sqrt( trA^2 - 4*detA ));
%         d2 = 1/2*(trA - sqrt( trA^2 - 4*detA ));
%         
%         Tr(k) = trA;
%         De(k) = detA;
%         Eig(k,1) = d1;
%         Eig(k,2) = d2;
%            
        outp(k,j,iinp) = 0
        ntr = 15;
        r = [.6; .2];
        R = zeros(Tds,N);
        tt = 0;
        for itrl = 1 : ntr
          
          for t = 1:L

            if t > 900 && t < 1100
              u = W*r + Io+[all_inp(iinp);0]; 
            else
              u = W*r + Io; 
            end
              K = feval(F,u);           
              r = r + dt*(-r./tau + K)+ sqrt(dt)*0.05*randn(2,1);

              R(t,:)=r;

          end

          outp(k,j,iinp) =(outp(k,j,iinp)+ (mean(R(931:1130,1))-mean(R(850:900,1)) )/mean(R(931:1130,1)))/ntr;
        end
%          Ampl(:,1)=abs(R(:,1) + 1i*R(:,2));
        
          
%         figure
%         plot(R)
%         set(gca,'ylim',[0 1])
%         %title(sprintf('wEE = %g,  wEI = %g, v1 = %2.3f + i*%2.3f',wEE,wEI,real(d1),imag(d1)))
%         title(sprintf('Ii = %2.2f, v1 = %2.3f + i*%2.3f',Ii,real(d1),imag(d1)))
%         text(Tds/2,.3,sprintf('TrA = %2.3f, D = %2.3f',trA,trA^2-4*detA))
  end
end
end
error('!')
%%
figure
xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.1 .2 .22 .69])
plot(Ies,Re,'k')
hold on
Ibif = Ies(find(abs(Re(:,1)-Re(:,2))>0.0001,1,'first')); %bifurcations
plot(Ibif(1)*[1 1],[0 1],'k')
plot(Ibif(end)*[1 1],[0 1],'k:')
xlabel('Control param. \itI_{E}','fontname','times')
ylabel('\itr_E    ','fontname','times','rotation',0)
text(-.2,1.1,'A','units','normalized','fontsize',14)
% 
axes('position',[.42 .2 .22 .69])
plot(rE,nullE(:,1),'k')
hold on
plot(nullI(:,1),rI,'color',[.5 .5 .5])
plot(trayectory(:,1),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{E}\rm = 0','fontname','times')
text(-.2,1.1,'B','units','normalized','fontsize',14)
% 
axes('position',[.74 .2 .22 .69])
plot(rE,nullE(:,2),'k')
hold on
plot(nullI(:,2),rI,'color',[.5 .5 .5])
plot(trayectory(:,2),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{E}\rm = 2.5','fontname','times')
text(-.2,1.1,'C','units','normalized','fontsize',14)

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifurcation_v%d.pdf',v))% 
% 
% 
% 
% 


