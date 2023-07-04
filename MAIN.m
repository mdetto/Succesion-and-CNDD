clear


JM = 100;
mp = [0 0.25];
c = 1.5;
c1 = c+1;
rmin = 0.5;
R = 500;
n1 = 1;
figure(1);clf
for j=1:2
csi =zeros(R*JM,1);
epsilon =zeros(R*JM,1);
N =zeros(R*JM,1);
TR =zeros(R*JM,1);
Ta =zeros(R*JM,1);
Tb =zeros(R*JM,1);
H = zeros(R,1);
for i=1:R
    i
[n,r,t,~,a,b]= LightCompetitionStrictCNDD_gost(JM,mp(j),rmin);
    
n2 = n1+length(r)-1;
csi(n1:n2,1) = a;
epsilon(n1:n2,1) = b;
N(n1:n2,1) = n;
Ta(n1:n2,1) = t;
Tb(n1:n2,1) = [1; t(1:end-1)];
TR(n1:n2,1) = r;
H(i) = length(n);
n1 = n2+1;
end

csi(n1:end,:)=[];
epsilon(n1:end,:)=[];
N(n1:end,:)=[];
TR(n1:end,:)=[];
Ta(n1:end,:)=[];
Tb(n1:end,:)=[];



rm = (1+rmin)/2;

subplot(2,2,j)
plot(csi,epsilon,'.')
% axis([1e-6 100 1e-6 100])
axis square
refline(1,0)
refline((c1-mp(j)*(1+rmin))/(c1+mp(j)*(1+rmin)),0)

% if j==1
%     mod2 = (sqrt(9-4.*csi.^2.*(c+2)./TR.*((c+2)./TR.^2.*csi-3))-3).^(1/2)...
%         ./(2*(c+2)./TR).^(1/2);
%     plot(csi,mod2,'r.')
%     pause
% end

xlabel('\xi','fontsize',15)
ylabel('\epsilon','fontsize',15)

title(['\alpha = ' num2str(mp(j))])

lambda = 2*(c+1)/(c+1+2*mp(j)*rm)*(1-rmin)/(JM+1);

Sc = (1-rmin)/lambda;

subplot(2,2,j+2)
histogram(H)
xline(Sc,'r-','linewidth',2)
xlabel('species richness')
legend('simulations','analytical solution')
axis square
pause(.1)
end
%% species richness as function of CNDD
clear

JM = [25 50 100];
mp = linspace(0,1,10);
c1 = 2.5;
rmin = 0.2;
R = 100;

S = zeros(length(mp),length(JM));
E = zeros(length(mp),length(JM));
S1 = zeros(length(mp),length(JM));
E1 = zeros(length(mp),length(JM));

for k=1:length(JM)
    for j=1:length(mp)
        j
        H = zeros(R,1);
        H1 = zeros(R,1);
        parfor i=1:R
            
            [n,r,t,~,a,b]= LightCompetitionStrictCNDD_gost(JM(k),mp(j),rmin);
            H(i,1) = length(n);
            
            % [r,n,t] = LightCompetitionStrictCNDD_linear_v2(JM(k)/(1-rmin),mp(j),rmin);
            % H1(i,1) = length(n);
        end
        
        S(j,k)=mean(H);
        E(j,k)=std(H)/sqrt(R);
        % S1(j,:)=mean(H1);
        % E1(j,1)=std(H1)/sqrt(R);
    end
    
end

figure(2);clf
for k=1:length(JM)
    % Sc = 1.2859 + (1+b*(1+rmin))*JM(k)/2;
    Sc = 1+(1+mp/c1*(1+rmin))*JM(k)/2;
    
    h = errorbar(mp,S(:,k),E(:,k),'.');hold all
    % h2 = errorbar(mp,S1,E1,'.');hold all
    plot(mp,Sc,'color',h.Color)
    xlabel('\beta','fontsize',15)
    ylabel('species richness')
    pause(.1)
end

h=get(gca,'children');
legend(h([1 3 5]),['{\itJ_M} = ' num2str(JM(3))],...
                  ['{\itJ_M} = ' num2str(JM(2))],...
                  ['{\itJ_M} = ' num2str(JM(1))])
legend('boxoff')    
set(gca,'xtick',0:0.2:1,'ytick',0:25:100)
xlim([-0.01 1.01])
axis square
%% graphical method
N = 3;
c = 1.5;
c1 = c+1;
d = 0.5;

r = [0.2 0.42 0.65 0.76 0.94 1.03].';
rp = (r.^-c - 1).^(-1/c);

t = linspace(0,1.2,10000);
tp = (t.^-c - 1).^(-1/c);

k = c.*r.^-c1;
LRS = 1/c*k.*t.^c1;


figure(1);clf
% subplot(211)
plot(t,LRS([1 2 4 5],:),'k-','linewidth',2);hold all
plot(t,LRS([3 6],:),'-','linewidth',2,'color',0.8*[1 1 1])
plot(r([1 2 4 5]),1,'ko','markersize',5,'MarkerFaceColor','k')
plot(r([3 6]),1,'ko','markersize',5,'MarkerFaceColor',0.8*[1 1 1])
yline(1,'k--')

xline(1,'r-')
for i=[1 2 4 5]
plot(r([i i]),[0 1],'r-')
end

axis([0 1.15 0 1.49])

% xlabel('patch age')
% ylabel('Tree height')
% xline(1,'-k','linewidth',2)
% set(gca,'xtick',[],'ytick',[])

% subplot(212)
% plot(tp,LRS,'k-','linewidth',2)
% hold all
% plot(rp,1,'ko','markersize',5,'MarkerFaceColor','k')
% yline(1,'k--')
% ylim([0 2])


%% patch-age distribution
    
    
t = linspace(0,10,1000);

t0 =[1 2 3];
p = zeros(1000,3);
m = zeros(1000,3);
for i=1:3
p(:,i) = c/(c+1)*(1+t.^c./t0(i).^c).^((-2*c+1)/c)./t0(i);
m(:,i) = (2*c+1)*t.^(c-1)./(t0(i).^c+t.^c);
end

subplot(211)
plot(t,p)
ylabel('patch-age distribution')

legend('{\itt}_0 = 1','{\itt}_0 = 3','{\itt}_0 = 3')
legend('boxoff')

subplot(212)
plot(t,m)
xlabel('patch age')
ylabel('disturbance rate')

%% break-even time and abundance simulations
clear
rmin = 0.2;
mp = [0.0 0.75];
JM = [25 100];
rep = 100000;
c = 1.5;

R = zeros(rep*JM(2),2,2);
Z = zeros(rep*JM(2),2,2);
Z2 = zeros(rep*JM(2),2,2);
for j=1:2
    for k=1:2
      
n1 = 1;
for i=1:rep
  [j k i]

[n,r,t]= LightCompetitionStrictCNDD_gost(JM(k),mp(j),rmin);

n2 = n1+length(r)-1;
t1 = [1; t(1:end-1)];
t2 = t;
Z(n1:n2,j,k) = n.*(t1+t2)/2;
Z2(n1:n2,j,k) = c/(c-1)*(t2.^(1-c)-t1.^(1-c))./(t2.^-c-t1.^-c).*n; % exact abundance
R(n1:n2,j,k) = r;
n1 = n2+1;


end

R(n1:end,j,k)=nan;
Z(n1:end,j,k)=nan;

    end
end

save('SAD2.mat','R','Z','Z2','JM','mp','rmin')





%% break-even time distributions
clear;load('SAD2.mat')
J = 20;
MP = linspace(0,0.75,J);
  
cmap = colormap('jet');
COL = zeros(length(MP),3);
COL(:,1) = interp1(1:64,cmap(:,1),linspace(1,64,J));
COL(:,2) = interp1(1:64,cmap(:,2),linspace(1,64,J));
COL(:,3) = interp1(1:64,cmap(:,3),linspace(1,64,J));

col(1,:) =  COL(1,:);
col(2,:) =  COL(J,:);

lgnd.position{1} = [0.0635  0.5835  0.0257  0.3416];
lgnd.position{2} = [0.0586  0.1100  0.0257  0.3416];

figure(1);clf  
for k=1:2
    subplot(2,2,(k-1)*2+1)
    for i=1:J
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),MP(i),rmin);
        plot(ri,pR,'r','linewidth',2,'color',COL(i,:));hold all
    end
    xlabel('break-even time')
    caxis(MP([1 J]))
    c = colorbar;
    c.Ticks = 0:0.25:0.75;
    c.Label.String = '\alpha';
    c.Label.FontSize = 14;
    c.Label.Rotation = 0;
    c.Position = lgnd.position{k};
    pos = get(gca,'Position');
    title(['analytical results {\itJ_M} = ' num2str(JM(k))])
    for j=1:2
        
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),mp(j),rmin);
        r0i = linspace(rmin,1,50);
        xi =linspace(rmin,1-(r0i(2)-r0i(1))/2,100);
        pR0 = interp1(ri,pR(1,:),xi);
        subplot(4,2,(k-1)*4+(j-1)*2+2)
        h=histogram(R(R(:,j,k)>0,j,k),r0i,'normalization','pdf');hold all
        plot(xi,pR0,'r','linewidth',2,'color',col(j,:));hold all
        set(h,'facealpha',0.2,'facecolor',col(j,:))
        set(gca,'xtick',[0.2 1],'ytick',[0 1 2],'ylim',[0 2.5])
        t = title(['\alpha = ' num2str(mp(j),2)]);
        t.Position = [0.6000    2.1917    0.0000];
        
        if k==2 && j==1
            p = get(gca,'position');
            p(2) = pos(2)+pos(4)-p(4);
            set(gca,'position',p);
        elseif k==1 && j==2
            p = get(gca,'position');
            p(2) = pos(2);
            set(gca,'position',p);
        end
        if j==1
            set(gca,'xtick',[],'ytick',[0 1 2])
        elseif j==2
            xlabel('break-even time','fontsize',12)
        end
            
        pause(.1)
    end
end

annotation(figure(1),'textbox',...
    [0.680 0.934 0.098 0.033],'String',{'simulations'},...
    'FontSize' ,14,'linestyle','none');


%% species abundance distribution
clear;load('SAD2.mat')
J = 20;
MP = linspace(0,0.75,J);
  
cmap = colormap('jet');
COL = zeros(length(MP),3);
COL(:,1) = interp1(1:64,cmap(:,1),linspace(1,64,J));
COL(:,2) = interp1(1:64,cmap(:,2),linspace(1,64,J));
COL(:,3) = interp1(1:64,cmap(:,3),linspace(1,64,J));

col(1,:) =  COL(1,:);
col(2,:) =  COL(J,:);

lgnd.position{1} = [0.0635  0.5835  0.0257  0.3416];
lgnd.position{2} = [0.0586  0.1100  0.0257  0.3416];

figure(2);clf  
for k=1:2
    subplot(2,2,(k-1)*2+1)
    for i=1:J
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),MP(i),rmin);
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',COL(i,:));hold all
    end
    xlabel('Log abundace')
    xlim([-11.5 -0.5])
    caxis(MP([1 J]))
    c = colorbar;
    c.Ticks = 0:0.25:0.75;
    c.Label.String = '\alpha';
    c.Label.FontSize = 14;
    c.Label.Rotation = 0;
    c.Position = lgnd.position{k};
    pos = get(gca,'Position');
    title(['{\itJ_M} = ' num2str(JM(k))])
    for j=1:2
        
        [M,ti,pR,ri,pZ,Lzi,SL,pN] = RenewalProcess(JM(k)/(1-rmin),mp(j),rmin);
        Lz0i = linspace(Lzi(1),Lzi(end),50);
%         xi =linspace(rmin,1-(r0i(2)-r0i(1))/2,100);
%         pZ0 = interp1(ri,pR(1,:),xi);
        subplot(4,2,(k-1)*4+(j-1)*2+2)
        h=histogram(log(Z2(Z2(:,j,k)>0,j,k)),Lz0i,'normalization','pdf');hold all
        plot(Lzi,pZ(1,:),'r','linewidth',2,'color',col(j,:));hold all
        set(h,'facealpha',0.2,'facecolor',col(j,:))
        set(gca,'xtick',[-10 -5 0],'ytick',[0 0.2 0.4],'ylim',[0 0.4],'xlim',[-11.5 -0.5])
        t = title(['\alpha = ' num2str(mp(j),2)]);
        t.Position = [-10.0755    0.3507    0.0000];
        
        if k==2 && j==1
            p = get(gca,'position');
            p(2) = pos(2)+pos(4)-p(4);
            set(gca,'position',p);
        elseif k==1 && j==2
            p = get(gca,'position');
            p(2) = pos(2);
            set(gca,'position',p);
        end
        if j==1
            set(gca,'xtick',[],'ytick',[0 1 2])
        elseif j==2
            xlabel('Log abundance','fontsize',12)
        end
            
        pause(.1)
    end
end

annotation(figure(2),'textbox',...
    [0.680 0.956 0.098 0.033],'String',{'simulations'},...
    'FontSize' ,14,'linestyle','none');

annotation(figure(2),'textbox',...
    [0.215 0.956 0.207 0.033],...
    'String','analytical results',...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');

%% find intercept of S ~ CNDD

clear


S = zeros(1,100);
S0 = zeros(1,100);
JM = round(logspace(1,3,100));
R = 1000;
mp = 0.5;
rmin = .2;

for k=1:length(JM)
    
k
H = zeros(R,1);
H0 = zeros(R,1);
parfor i=1:R

% [n,r,t]= LightCompetitionStrictCNDD_gost(JM(k),mp,rmin);
% H(i,1) = length(n);
[n,r,t]= LightCompetitionStrictCNDD_linear_v2(JM(k)/(1-rmin),mp,rmin);
H0(i,1) = length(n);

end

% S(k)=mean(H);
S0(k)=mean(H0);
end


save('S_intecept2.mat','S','S0','JM','mp','rmin')

clf
subplot(121)
load S_intecept1.mat
semilogx(JM,1/2+1/2*JM.^-1);hold all
semilogx(JM,(1+mp/2.5*(1+rmin))/2+JM.^-1)
semilogx(JM,S0./JM,'bo','markersize',4)
semilogx(JM,S./JM,'ro','markersize',4)

xlabel('\itJ_M')
ylabel('s/{\itJ_M}')
legend('linearized anlytical solution','corrected solution',...
       'linearized simulations','nonlinear simulations')
title('\alpha = 0')
axis square

subplot(122)
load S_intecept2.mat
semilogx(JM,(1/2-mp/2.5)./JM + (1/2 +1/2*mp/2.5*(1+rmin)));hold all
semilogx(JM,(1+mp/2.5*(1+rmin))/2+JM.^-1)
semilogx(JM,S0./JM,'bo','markersize',4)
semilogx(JM,S./JM,'ro','markersize',4)
xlabel('\itJ_M')
% ylabel('s/{\itJ_M}')
% legend('simuations','1/2 + 1/2{\itJ_M}^-^1','1/2 + {\itJ_M}^-^1')
title('\alpha = 0.5')
axis square

%% adatptive evolution
rmin = 0.2;
R = 100;
J = 75;
Tmax = logspace(1,3.8,J);
mp = [0 0.1 0.2];

s = ones(R,J,length(mp));
T = zeros(R,J,length(mp));
TR = zeros(R*100,length(mp));

for k=1:length(mp)
    TR0 = cell(R,1);
parfor i=1:R
    [k i]
    
    s0 = ones(1,J);
    T0 = zeros(1,J);
    [r,t,T0(1)] = AdpaptEvol(rmin,1,mp(k),0,Tmax(1),'Tmax');
    s0(1) = length(r);
    for j=2:J
    [r,t,T0(j)] = AdpaptEvol(r,t,mp(k),T0(j-1),Tmax(j),'Tmax');
    s0(j) = length(r);
    end
    T(i,:,k)=T0;
    s(i,:,k)=s0;
    
    TR0{i}=r;

end
  n1=1;
  for i=1:R
   n2 = n1+length(TR0{i})-1;
   TR(n1:n2,k)=TR0{i};
   n1 = n2+1;
  end
end

% save('AdaptiveEvolution.mat','T','TR','s','R','rmin','mp')
%% Nadaraya–Watson kernel regression
load('AdaptiveEvolution_v2.0.mat')
clf
col(1,:) = [0         0.4471    0.7412];
col(2,:) = [0.8510    0.3255    0.0980];
col(3,:) = [0.4667    0.6745    0.1882];
% col(4,:) = [0.9294    0.6941    0.1255];
lgnd = cell(3,1);
Lti = linspace(1.6,3.5,75);
ri = linspace(rmin,1,30);
for i=1:3
T0 = T(:,:,i);
s0 = s(:,:,i);
x = log2(s0(:));
y = log10(T0(:));

pdf0 = ksdensity(y,Lti);
pdf1 = ksdensity(y,Lti,'weights',x);
pdf2 = ksdensity(y,Lti,'weights',x.^2);

mi = pdf1./pdf0*mean(x);
vi = pdf2./pdf0*mean(x.^2) - mi.^2;

subplot(211)
gray_area(Lti,mi,mi-sqrt(vi/R),mi+1*sqrt(vi/R),'color',col(i,:))
hold all


lgnd{4-i} = ['\alpha = ' num2str(mp(i))];
end
subplot(211)
set(gca,'xtick',1:4,'xticklabel',10.^(1:4),...
        'ytick',0:6,'yticklabel',2.^(0:6),'xlim',[1.2 4.1])
h = get(gca,'children');
legend(h([1 3 5]),lgnd)
legend('boxoff')
xlabel('evolutionary time')
ylabel('species richness')

Tstar = (1-rmin)/2/(0.1*sqrt(2/pi))/(0.01*2);

%%  adatptive evolution distributions

rmin = 0.2;
R = 200;
J = 75;
Tmax = logspace(1,3.8,J);
mp = [0 0.1 0.2];
S0 = [20 40];
TR1 = zeros(R,S0(1),length(mp));
TR2 = zeros(R,S0(2),length(mp));

for k=1:length(mp)
parfor i=1:R
    [k i]
    
    [r,t,T0] = AdpaptEvol(rmin,1,mp(k),0,S0(1),'Smax');
    TR1(i,:,k)=r;
    
    [r,t,T0] = AdpaptEvol(r,t,mp(k),T0,S0(2),'Smax');
    TR2(i,:,k)=r;

end

end

%% plot analysis above
col(1,:) = [0         0.4471    0.7412];
col(2,:) = [0.8510    0.3255    0.0980];
col(3,:) = [0.4667    0.6745    0.1882];

subplot(223);cla
ri = linspace(rmin+.01,1-.01,100);
N = zeros(length(ri),3);
for i=1:3
    TR0 = TR1(:,:,i);
N(:,i) = ksdensity(TR0(:),ri,'support',[rmin 1]);
end

for j=1:3
plot(ri,N(:,j),'color',col(j,:),'linewidth',2);hold all
end
title(['{\its} = ' num2str(S0(1))])
ylim([0 12])
xlabel('break-even time')
ylabel('pdf')

subplot(224);cla
ri = linspace(rmin+.01,1-.01,50);
N = zeros(length(ri),3);
for i=1:3
    TR0 = TR2(:,:,i);
N(:,i) = ksdensity(TR0(:),ri,'support',[rmin 1]);
end
for j=1:3
plot(ri,N(:,j),'color',col(j,:),'linewidth',2);hold all
end
title(['{\its} = ' num2str(S0(2))])
ylim([0 12])
xlabel('break-even time')
%% forest structure
clear

figure(5);clf
t0 = 1;
S = [25 50 100];
mp = 0;
c = 1.5;
c1 = c+1;
rmin = 0.35;
R = 200;
g = 1;

b = 0.3;

xmax = g*1;
xmin = g*rmin^(1-b)*0.8;
x = linspace(xmin,xmax,1000);

Z1 = zeros(length(x),R);
H = zeros(R,1);
for m=1:3
for j=1:R
    j
    [n,r,t]= LightCompetitionStrictCNDD_gost(S(m),mp,rmin);
    f = g*r.^-b;
    ta = t;
    tb = [1;t(1:end-1)];
    if length(n)>1
        
        W1 = heaviside(f.*ta-x);
        W2 = heaviside(f.*tb-x);
        
        Z1(:,j) = sum(n./f.*W1 + ((x./f).^-c - tb.^-c)./(f.*c.*r.^-c1).*(1-W1).*W2);
        H(j)=length(n);
    else
        
        W1 = heaviside(r.*t-x);
        W2 = heaviside(r.*1-x);
        
        Z1(:,j) = (n.*W1 + ((x./f).^-c - tb.^-c)./(c*r.^-c1).*(1-W1).*W2)./f;
    end
    
    
end


semilogx(x,mean(Z1,2),'-');hold all
    pause(.1)
end
    
if b==0
Z0 = (1-x/g)/g;
Z0(x<g*rmin) = (1-rmin)/g;
elseif b==-1
Z0 = -log(sqrt(x/g))/g;
Z0(x<g*rmin^2) = -log(rmin)/g;
elseif b>0
Z0 = 1/(b+1)*(1-(x/g).^((b+1)/(1-b)))/g;
Z0(x<g*rmin^(1-b)) =  1/(b+1)*(1-rmin^(b+1))/g;
end

% semilogx(x,mean(Z1,2),'-');hold all
loglog(x,Z0,'-k','linewidth',2)

xlabel('tree size (x)')
ylabel('abundance')
legend(['S  = ' num2str(S(1))],...
       ['S  = ' num2str(S(2))],...
       ['S  = ' num2str(S(3))],...
        'analytical model');