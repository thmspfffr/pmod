clear

% pmod_wc_original_with_i


% fixed params:

% VERSION 1:
% -------------------------------------
v = 2;
SIGMA = [0.002 0.01];
GAIN = 0:0.1:0.20;
Ies  = -2:0.1:0;
Iis  = 0:0.25:5;
ASYM = 0.9:0.1:1.1;
% -------------------------------------

%%
N = 4;
% Connectivity:

wII=4;
wIE=16;
wEI=10;
wEE=12;

g = 1;

tauE = 1;
tauI = 1; %1
tau = [tauE;tauI;tauE;tauI];

Io=zeros(N,1);

dt=0.05;
tmax = 300; %5000;
tspan=0:dt:tmax;
L = length(tspan);

ds = 1;
Tds = length(0:ds*dt:tmax)-1;

% Connectivity:
W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [g 0; 0 0];
W21 = [g 0; 0 0];

W = [W11 W12; W21 W22];

numIes = length(Ies);

nTrials = 40;

uu = -4:.1:4;
Transf = zeros(length(uu),2,4);
%%
  for igain = 1 : length(GAIN)
    for iinh = 1 : length(Iis)
      for inoise = 1 : length(SIGMA)
        for iasym = 1 : length(ASYM)


          if ~exist(['~/pmod/proc/' sprintf('pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d_proc.txt',igain,inoise,iinh,iasym,v)])
            system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d_proc.txt',igain,inoise,iinh,iasym,v)]);
          else
            continue
          end
        
          
          asym = ASYM(iasym);
          
          Gain = GAIN(igain);
          
          Ii = Iis(iinh);
          Io(2) = -Ii;
          Io(4) = -Ii;
          sdt = sqrt(dt)*SIGMA(inoise);
          
          % transfer functions:
          % gains are given by 1/aE and 1/aI
          gE = 1+asym*Gain;
          gI = 1+(1/asym)*Gain;
          
          aE = 1/gE;
          aI = 1/gI;
          
          Fe = @(x) 1./(1 + exp(-x/aE) );
          Fi = @(x) 1./(1 + exp(-x/aI) );
          F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];
          
          % store functions:
          Transf(:,1) = feval(Fe,uu);
          Transf(:,2) = feval(Fi,uu);
          
          
          for iexc = 1:numIes
            
            iexc
            Ie = Ies(iexc);
            Io(1) = Ie;
            Io(3) = Ie;
            
            Rc(iexc)     = 0;
            RcAmpl(iexc) = 0;
            
            dfa(iexc,:) = [0 0];
            
            for trial = 1:nTrials
              
              r = rand(N,1);
              R = zeros(Tds,N);
              tt = 0;
              
              % Warm-up:
              for t = 1:5000
                u = W*r + Io;
                K = feval(F,u);
                r = r + dt*(-r + K)./tau + sdt*randn(N,1);
              end
              
              for t = 1:L
                u = W*r + Io;
                K = feval(F,u);
                r = r + dt*(-r+ K)./tau + sdt*randn(N,1);
                if mod(t,ds)==0
                  tt=tt+1;
                  R(tt,:)=r;
                end
              end
              
              RstatE = R(:,[1 3]);
              Ampl(:,1)=abs(R(:,1) + 1i*R(:,2));
              Ampl(:,2)=abs(R(:,3) + 1i*R(:,4));
              
              rc = corr(RstatE);
              Rc(iexc) = Rc(iexc) + rc(2)/nTrials;
              
              rc = corr(Ampl);
              RcAmpl(iexc) = RcAmpl(iexc) + rc(2)/nTrials;
              
              d = tp_dfa(Ampl(:,1),[3 30],(1/dt),0.5,20);
              
              dfa(iexc,1) = dfa(iexc,1) + d.exp/nTrials;
              
              d = tp_dfa(Ampl(:,2),[3 30],(1/dt),0.5,20);
              dfa(iexc,2) = dfa(iexc,2) +  d.exp/nTrials;
              
            end
          end
          
          fprintf('Saving...\n');
          save(sprintf('~/pmod/proc/pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d.mat',igain,inoise,iinh,iasym,v),'Rc','RcAmpl','dfa','-v7.3')
          
          try
            load(sprintf('~/pmod/proc/pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d.mat',igain,inoise,iinh,iasym,v))
          catch me
            pause(randi(10))
            save(sprintf('~/pmod/proc/pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d.mat',igain,inoise,iinh,iasym,v),'Rc','RcAmpl','dfa','-v7.3')
          end
          
          clear RcAmpl Rc rc
          
        end
      end
    end
  end


error('!')


%%
clear alldfa alll
  for igain = 1 : length(GAIN)
    for iinh = 1 : length(Iis)
      for inoise = 1 : length(SIGMA)
        for iasym = 1 : length(ASYM)
          
          load(sprintf('~/pmod/proc/pmod_wc_orig_with_i_gain%d_noise%d_inh%d_asym%d_v%d.mat',igain,inoise,iinh,iasym,v))
          alll(:,igain,iinh,inoise,iasym) = RcAmpl;
          alldfa(:,igain,iinh,inoise,iasym) = mean(dfa,2);

        end
      end
      
    end
  end
  
  alll = squeeze(alll);
  alldfa = squeeze(alldfa);

%%
alll=alldfa;
figure;

subplot(1,4,1)
imagesc(squeeze(alll(:,1,:)),[0.5 1])
axis square
set(gca,'xtick',[1:10:length(Iis)],'xticklabel',Iis(1:10:length(Iis)))
ylabel('Excitation'); xlabel('Inhibition')

% subplot(3,4,5)
% imagesc(alll(:,2:end,3)'-alll(:,1,3)',[-0.2 0.2])
% axis square
% set(gca,'xtick',[1:10:length(Ies)],'xticklabel',Ies(1:10:length(Ies)))
% ylabel('Gain')

subplot(1,4,2)
imagesc(squeeze(alll(:,3,:)),[0.5 1])
axis square
set(gca,'xtick',[1:10:length(Iis)],'xticklabel',Iis(1:10:length(Iis)))
ylabel('Excitation'); xlabel('Inhibition')

% subplot(3,4,6)
% imagesc(alll(:,2:end,1)'-alll(:,1,1)',[-0.2 0.2])
% axis square
% set(gca,'xtick',[1:10:length(Ies)],'xticklabel',Ies(1:10:length(Ies)))
% % ylabel('Gain')

subplot(1,4,3)
imagesc(squeeze(alll(:,5,:)),[0.5 1])
axis square
set(gca,'xtick',[1:10:length(Iis)],'xticklabel',Iis(1:10:length(Iis)))
ylabel('Excitation'); xlabel('Inhibition')

% subplot(3,4,7)
% imagesc(alll(:,2:end,5)'-alll(:,1,5)',[-0.2 0.2])
% axis square
% set(gca,'xtick',[1:10:length(Ies)],'xticklabel',Ies(1:10:length(Ies)))

subplot(1,4,4)
imagesc(squeeze(alll(:,5,:))-squeeze(alll(:,1,:)),[-0.2 0.2])
% colormap(cmap);
axis square
set(gca,'xtick',[1:10:length(Iis)],'xticklabel',Iis(1:10:length(Iis)))
ylabel('Excitation'); xlabel('Inhibition')

colormap(cmap)

% subplot(3,4,5)
% imagesc(all(:,:,3)',[0, 1])
% axis square
% set(gca,'xtick',[1:6:length(Ies)],'xticklabel',Ies(1:6:25))
% ylabel('Gain')
% 
% subplot(3,4,6)
% imagesc(all(:,:,1)',[0, 1])
% axis square
% set(gca,'xtick',[1:6:length(Ies)],'xticklabel',Ies(1:6:25))
% 
% subplot(3,4,7)
% imagesc(all(:,:,5)',[0, 1])
% colormap(cmap); axis square
% set(gca,'xtick',[1:6:length(Ies)],'xticklabel',Ies(1:6:25))
% 
% 
% subplot(3,4,8)
% imagesc(all(:,:,5)'-all(:,:,1)',[-0.5 0.5])
% axis square
% set(gca,'xtick',[1:6:length(Ies)],'xticklabel',Ies(1:6:25))
%%
xlabel('Excitatory input'); ylabel('Gain')

subplot(3,3,8)
imagesc(all(:,:,3)'-all(:,:,1,1,1)',[-0.5 0.5])
axis square
set(gca,'xtick',[1:6:length(Ies)],'xticklabel',Ies(1:6:25))

xlabel('Excitatory input')

% subplot(3,3,9)
% imagesc([all(:,:,3,5,2)'-all(:,:,2,5,1)']-[all(:,:,1,5,2)'-all(:,:,2,5,1)]',[-0.5 0.5])
% colormap(cmap); axis square
% set(gca,'xtick',[1:5:length(Ies)],'xticklabel',Ies(1:5:25))


%%
figure
xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.07 .2 .22 .69])
plot([Ies(1) Ies(end)],[0 0],'k:')
hold on
plot(Ies,Rc,'.-')
xlabel('Control param. \itI_{E}','fontname','times')
ylabel('corr. EE (c_{E-E })','fontname','times')
l = {'no gain','\deltaG_{E} = \deltaG_{I}',...
  '\deltaG_{E} > \deltaG_{I}','\deltaG_{E} < \deltaG_{I}'};
col=lines(4);
set(gca,'ylim',[-.4 1.001],'fontsize',9)
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
for i=1:4
  plot(xlim(1)+.83*(xlim(2)-xlim(1)),ylim(1)+(.38-(i-1)*.1)*(ylim(2)-ylim(1)),...
    'o','color',col(i,:),'markerfacecolor',col(i,:))
  text(.9,.4-(i-1)*.11,l{i},'units','normalized','backgroundcolor','w','fontsize',8)
end
plot(-1.9*[1 1],ylim,'k:')
box off
text(-.22,1.05,'A','units','normalized','fontsize',14)


Q = 2;
Iopt = -3.2;
% Iopt = find( Rc(:,1)-Rc(:,Q)>.05 & Rc(:,1)<.7 & Rc(:,Q)>0 );

axes('position',[.405 .6 .19 .3])
plot(uu,Transf(:,1,1),'b')
hold on
plot(uu,Transf(:,1,Q),'k')
ylabel('F_{E}  ','fontname','times','rotation',0)
set(gca,'fontsize',9)

text(-.25,1.1,'B','units','normalized','fontsize',14)

axes('position',[.405 .2 .19 .3])
plot(uu,Transf(:,2,1),'b')
hold on
plot(uu,Transf(:,2,Q),'k')
set(gca,'fontsize',9)
xlabel('Input','fontname','times')
ylabel('F_{I}  ','fontname','times','rotation',0)


aE = 1;
aI = 1;

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

Ii = Iopt; %Ies(Iopt(1));
Io(1) = Ie;
Io(3) = Ie;

r = rand(N,1);
R = zeros(Tds,N);
tt = 0;
for t = 1:L
  u = W*r + Io;
  K = feval(F,u);
  r = r + dt*(-r./tau + K) + sdt*randn(N,1);
  if mod(t,ds)==0
    tt=tt+1;
    R(tt,:)=r;
  end
end

RstatE = R(tt-1999:tt,[1 3]);

axes('position',[.68 .6 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)
text(-.22,1.1,'C','units','normalized','fontsize',14)

aE = 1/(1+p*Gain);
aI = 1/(1+Gain);

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

Ie = Iopt; %Ies(Iopt(1));
Io(1) = Ie;
Io(3) = Ie;

r = rand(N,1);
R = zeros(Tds,N);
tt = 0;
for t = 1:L
  u = W*r + Io;
  K = feval(F,u);
  r = r + dt*(-r./tau + K) + sdt*randn(N,1);
  if mod(t,ds)==0
    tt=tt+1;
    R(tt,:)=r;
  end
end

RstatE = R(tt-1999:tt,[1 3]);

axes('position',[.68 .2 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)

text(.4,.8,sprintf('I_{E} = %2.2f',Ii),'units','normalized')

% export_fig Fig_TwoWC_Correlations_Ie.tiff -r300


