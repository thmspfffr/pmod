clear all;

% fixed params:

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


dt=0.05;
tmax = 1000;
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


Iis = -30:.1:30;
Ies = -30:.1:30;

numIis = length(Iis);
numIes = length(Ies);

Re  = zeros(numIis,2);
Tr  = zeros(numIis,1);
De  = zeros(numIis,1);
Eig = zeros(numIis,N);
for k = 1:numIis
  
  
%   if ~exist(['~/pmod/proc/' sprintf('pmod_wc_bif_k%d.txt',k)])
%     system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_bif_k%d.txt',k)]);
%   else
%     continue
%   end
  
  Ii = Iis(k);
  Io(2) = Ii;
  clear Ie
  
  for l = 1 : numIes
    
    
    Ie = Ies(l);
    Io(1) = Ie;
    
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
    %         Re(k,1) = max(RstatE);
    %         Re(k,2) = min(RstatE);
    Re(l,1) = max(RstatE);
    Re(l,2) = min(RstatE);
    %         figure
    %         plot(R)
    %         set(gca,'ylim',[0 1])
    %         %title(sprintf('wEE = %g,  wEI = %g, v1 = %2.3f + i*%2.3f',wEE,wEI,real(d1),imag(d1)))
    %         title(sprintf('Ii = %2.2f, v1 = %2.3f + i*%2.3f',Ii,real(d1),imag(d1)))
    %         text(Tds/2,.3,sprintf('TrA = %2.3f, D = %2.3f',trA,trA^2-4*detA))
    
    %     if Ii == -4
    %
    %         % null-clines (lines for which dr/dt=0)
    %         % 0 = -re/tauE + Fe( wEE*re - wEI*ri + Ie );
    %         % 0 = -ri/tauI + Fi( wIE*re - wII*ri + Ii );
    %         % FinvE( re/tauE ) = wEE*re - wEI*ri + Ie;
    %         % FinvI( ri/tauI ) = wIE*re - WII*ri + Ii;
    %         rE = 0:0.01:1;
    %         rI = 0:0.01:1;
    %         nullE(:,1) = ( wEE*rE + Ie - feval( FinvE, rE/tauE ) )/wEI; % = rI
    %         nullI(:,1) = ( wII*rI - Ii + feval( FinvI, rI/tauI ) )/wIE; % = rE
    %
    %         trayectory(:,1) = R(1:200,1) + 1i*R(1:200,2);
    %
    %     elseif Ii == -6
    %
    %         rE = 0:0.01:1;
    %         rI = 0:0.01:1;
    %         nullE(:,2) = ( wEE*rE + Ie - feval( FinvE, rE/tauE ) )/wEI; % = rI
    %         nullI(:,2) = ( wII*rI - Ii + feval( FinvI, rI/tauI ) )/wIE; % = rE
    %
    %         trayectory(:,2) = R(1:200,1) + 1i*R(1:200,2);
    %
  end
  
  try
    figure
    xSize = 15; ySize = 4.8;
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperUnits','centimeters')
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')
    
    axes('position',[.1 .2 .22 .69])
    plot(Ies,Re,'k')
    hold on
    Ibif = Ies(abs(Re(:,1)-Re(:,2))>0.0001); %bifurcations
    plot(Ibif(1)*[1 1],[0 1],'k:')
    plot(Ibif(end)*[1 1],[0 1],'k:')
    xlabel('Control param. \itI_{I}','fontname','times')
    ylabel('\itr_E    ','fontname','times','rotation',0)
    text(-.2,1.1,'A','units','normalized','fontsize',14)
    
    print(gcf,'-djpeg100',sprintf('~/pmod/plots/pmod_wc_bif_l%d.jpeg',l));
  catch me
  end
end

% axes('position',[.42 .2 .22 .69])
% plot(rE,nullE(:,1),'k')
% hold on
% plot(nullI(:,1),rI,'color',[.5 .5 .5])
% plot(trayectory(:,1),'k:')
% xlabel('\itr_E','fontname','times')
% ylabel('\itr_I    ','rotation',0,'fontname','times')
% set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
% title('\itI_{I}\rm = -4','fontname','times')
% text(-.2,1.1,'B','units','normalized','fontsize',14)
%
% axes('position',[.74 .2 .22 .69])
% plot(rE,nullE(:,2),'k')
% hold on
% plot(nullI(:,2),rI,'color',[.5 .5 .5])
% plot(trayectory(:,2),'k:')
% xlabel('\itr_E','fontname','times')
% ylabel('\itr_I    ','rotation',0,'fontname','times')
% set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
% title('\itI_{I}\rm = -6','fontname','times')
% text(-.2,1.1,'C','units','normalized','fontsize',14)

% export_fig Fig_WC_Bifurcations.tiff -r300