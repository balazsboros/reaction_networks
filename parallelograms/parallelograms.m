dirname = 'C:/bboros/Dropbox/dfc1thm/3d/parallelogram_paper/jpg/';
colors = {'#EDB120','#77AC30','#0000FF'};
%         stable    unstable  semistable
%         ocher     green     blue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Section 3.1 Supercritical Andronov-Hopf bifurcation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: the parameters are (k1,k2,k3,k4,t,gamma)

% Here gamma = 1.
% This is the first plot in Figure 3.
p   = [ 8 ; 1/8 ; 1 ; 1 ; 1.1 ; 1 ];
eq  = [p(5)^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));p(5)^p(6)];
ap  = 5;
x0  = init_EP_EP(@ode_prlgrm_irrev,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',  0.1);
opt = contset(opt,'Singularities',  1);
opt = contset(opt,'Backward',       1);
[x1,~,s1] = cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x     = x1(1:2,s1(2).index);
hh    = 1e-6;
ntst  = 20;
ncol  = 4;
opt   = contset;
opt   = contset(opt,'MaxStepSize',     1);
opt   = contset(opt,'InitStepSize', 0.01);
opt   = contset(opt,'Singularities',   1);
opt   = contset(opt,'Multipliers',     1);
opt   = contset(opt,'Backward',        1);
opt   = contset(opt,'Adapt',           1);
opt   = contset(opt,'MaxNumPoints',  100);

[x0,v0] = init_H_LC(@ode_prlgrm_irrev, x, p, ap, hh, ntst, ncol);
xlc = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
eqx = t.^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));
eqy = t.^p(6);
eqz = t.^2*(p(3)*p(4)/p(1)/p(2))^(1/2/p(6));
c = p(6)*eqx + p(6)*eqy + 2*eqz;
z = (c-p(6)*x-p(6)*y)/2;

t = linspace(0,min(t),1000);
eqx2 = t.^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));
eqy2 = t;
eqz2 = t.^2*(p(3)*p(4)/p(1)/p(2))^(1/2/p(6));

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 2
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = 1 : 1 : size(xlc,2)
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
plot3(eqx,eqy,eqz,'Color',colors{2},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([170,20]);
legend({'stable','unstable'},'FontSize',12,'Location','NorthEast');
title(['L_1 < 0 for \gamma = ',num2str(p(6)),', \kappa_1 = ',num2str(p(1)),...
  ', \kappa_2 = 1/8, \kappa_3 = ',num2str(p(3)),', \kappa_4 = ',num2str(p(4))],'FontSize',12);
saveas(fig,[dirname,'1LC_gamma1.jpg']);
%%

% Here gamma = 2 (the homogeneous case).
% This is Figure 2.

p = [ 16 ; 1/16 ; 1 ; 1 ; 1 ; 2];
eq = [p(5)^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));p(5)^p(6);...
  (p(3)*p(4)/p(1)/p(2))^(1/2/p(6))*p(5)^2];
c = sum(eq.*[p(6);p(6);2]);

% First we find the stable limit cycle numerically.
odecell = ode_prlgrm_irrev;
options = odeset('RelTol',1e-8);
[time,yfwd] = ode45(odecell{2},linspace(0,100,50000),eq(1:2)*1.03,...
  options,p(1),p(2),p(3),p(4),p(5),p(6));
yfwd = [yfwd,(c-p(6)*yfwd(:,1)-p(6)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Here we find a time that is slightly larger than the period.
[time,yfwd] = ode45(odecell{2},linspace(0,1.04,1e3),...
  yfwd(end,1:2),options,p(1),p(2),p(3),p(4),p(5),p(6));
yfwd = [yfwd,(c-p(6)*yfwd(:,1)-p(6)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

x = yfwd(:,1);
y = yfwd(:,2);
z = yfwd(:,3);
t = linspace(0,p(5),101);
eqx = t.^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));
eqy = t.^p(6);
eqz = t.^2*(p(3)*p(4)/p(1)/p(2))^(1/2/p(6));

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 2
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = 1 : numel(t)
  plot3(x/p(5)*t(i),y/p(5)*t(i),z/p(5)*t(i),'Color',colors{1});
end
plot3(eqx,eqy,eqz,'Color',colors{2},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([170,20]);
legend({'stable','unstable'},'FontSize',12,'Location','NorthEast');
title(['The limit cycles for \gamma = ',num2str(p(6)),', \kappa_1 = ',num2str(p(1)),...
  ', \kappa_2 = 1/16, \kappa_3 = ',num2str(p(3)),', \kappa_4 = ',num2str(p(4))],'FontSize',12);
saveas(fig,[dirname,'1LC_gamma2.jpg']);

%%

% Here gamma = 3.
% This is the second plot in Figure 3.

p  = [ 16 ; 1/16 ; 1 ; 1 ; 1.02 ; 3];
eq = [p(5)^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));p(5)^p(6)];
ap = 5;
x0 = init_EP_EP(@ode_prlgrm_irrev,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',  0.01);
opt = contset(opt,'Singularities',   1);
opt = contset(opt,'Backward',        1);
[x1,~,s1] = cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x = x1(1:2,s1(2).index);
hh = 1e-6;
ntst = 20;
ncol = 4;
opt = contset;
opt = contset(opt,'MaxNumPoints', 1000);
opt = contset(opt,'MaxStepSize',     1);
opt = contset(opt,'InitStepSize',  0.1);
opt = contset(opt,'Singularities',   1);
opt = contset(opt,'Multipliers',     1);
opt = contset(opt,'Backward',        1);
opt = contset(opt,'Adapt',           1);

[x0,v0] = init_H_LC(@ode_prlgrm_irrev, x, p, ap, hh, ntst, ncol);
xlc1 = cont(@limitcycle,x0,v0,opt);

% We continue the last limit cycle (with increased step size).
p(5) = xlc1(end,end);
X0 = [xlc1(1,end),xlc1(2,end)];
tmax = 1.2*xlc1(end-1,end);
odecell = ode_prlgrm_irrev;
[time,yfwd] = ode45(odecell{2},linspace(0,tmax,2000),X0,options,...
  p(1),p(2),p(3),p(4),p(5),p(6));

tolerance = 1e-4;
[x0,v0]=initOrbLC(@ode_prlgrm_irrev,...
  time,yfwd,p,ap,ntst,ncol,tolerance);
opt = contset(opt,'MaxStepSize',100);
opt = contset(opt,'InitStepSize',10);
xlc2 = cont(@limitcycle,x0,v0,opt);

% We continue the last limit cycle (with even further increased step size).
p(5) = xlc2(end,end);
X0 = [xlc2(1,end),xlc2(2,end)];
tmax = 1.2*xlc2(end-1,end);
[time,yfwd] = ode45(odecell{2},linspace(0,tmax,2000),X0,options,...
  p(1),p(2),p(3),p(4),p(5),p(6));

tolerance = 1e-4;
[x0,v0]=initOrbLC(@ode_prlgrm_irrev,...
  time,yfwd,p,ap,ntst,ncol,tolerance);
opt = contset(opt,'MaxStepSize',10000);
opt = contset(opt,'InitStepSize',100);
xlc3 = cont(@limitcycle,x0,v0,opt);

xlc = [xlc1,xlc2,xlc3];
nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
eqx = t.^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));
eqy = t.^p(6);
eqz = t.^2*(p(3)*p(4)/p(1)/p(2))^(1/2/p(6));
c = p(6)*eqx + p(6)*eqy + 2*eqz;
z = (c-p(6)*x-p(6)*y)/2;
maxX = max(max(x));
maxY = max(max(y));
maxZ = max(max(z));

thr = 0.15;
latest = eqx(1)+eqy(1)+eqz(1); idx = 1;
for i = 2 : numel(t)
  if abs(eqx(i)+eqy(i)+eqz(i)-latest)>thr
    idx = [idx,i]; latest = eqx(i)+eqy(i)+eqz(i);
  end
end

t2 = linspace(max(t),1.1*max(t),100);
eqx2 = t2.^p(6)*sqrt(p(1)*p(4)/p(2)/p(3));
eqy2 = t2.^p(6);
eqz2 = t2.^2*(p(3)*p(4)/p(1)/p(2))^(1/2/p(6));

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 2
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = idx
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
plot3(eqx,eqy,eqz,'Color',colors{2},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([170,20]);
legend({'stable','unstable'},...
  'FontSize',12,'Location','NorthEast');
title(['L_1 < 0 for \gamma = ',num2str(p(6)),', \kappa_1 = ',num2str(p(1)),...
  ', \kappa_2 = 1/16, \kappa_3 = ',num2str(p(3)),', \kappa_4 = ',num2str(p(4))],'FontSize',12);
saveas(fig,[dirname,'1LC_gamma3.jpg']);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Section 3.2 Subcritical Andronov-Hopf bifurcation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: the parameters are (a,b,t,gamma)

% Here gamma = 1.
% This is the third plot in Figure 5.

p = [ 0.005 ; 0.14 ; 0.9 ; 1 ];
eq = [2*p(3)^p(4);1/2*p(3)^p(4)];
ap = 3;

x0 = init_EP_EP(@ode_prlgrm_rev,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',  0.001);
opt = contset(opt,'Singularities',    1);
opt = contset(opt,'Backward',         1);
[x1,~,s1] = cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x     = x1(1:2,s1(2).index);
hh    = 1e-6;
ntst  = 40;
ncol  = 4;
opt  = contset(opt,'MaxStepSize',1);
opt  = contset(opt,'Multipliers',1);
opt  = contset(opt,'Adapt',1);
opt  = contset(opt,'MaxNumPoints',220);

[x0,v0] = init_H_LC(@ode_prlgrm_rev, x, p, ap, hh, ntst, ncol);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
equilx = 2*t.^p(4);
equily = 1/2*t.^p(4);
equilz = t.^2;
c = p(4)*equilx + p(4)*equily + 2*equilz;
z = (c-p(4)*x-p(4)*y)/2;

t1 = linspace(0,t(1),1000);
eqx1 = 2*t1.^p(4);
eqy1 = 1/2*t1.^p(4);
eqz1 = t1.^2;
t2 = linspace(t(1),max(t),1000);
eqx2 = 2*t2.^p(4);
eqy2 = 1/2*t2.^p(4);
eqz2 = t2.^2;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
idx2 = slc(3).index;
idx1 = 1 : idx2-1;
idx3 = idx2+1 : size(x,2);
latest = t(1);
for i = idx1(1) : idx1(end)
  if abs(t(i)-latest)>0.002
    plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
    latest = t(i);
  end
end
plot3(x(:,idx2),y(:,idx2),z(:,idx2),'Color',colors{3},'LineWidth',1.5);
latest = t(idx2);
for i = idx3(1) : 1 : idx3(end)
  if abs(t(i)-latest)>0.005
    plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
    latest = t(i);
  end
end
plot3(eqx1,eqy1,eqz1,'Color',colors{1},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{2},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([-170,20]);
legend({'stable','unstable','semistable'},...
  'FontSize',12,'Location','NorthWest');
title(['L_1 > 0 for \gamma = ',num2str(p(4)),', a = ',num2str(p(1)),', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'2LC_gamma1.jpg']);

%%

% Here gamma = 2 (the homogeneous case).

% First check for what value of b the fold bifurcation
% of the limit cycles occur.
% (Denoted by LPC in MATCONT for "Limit Point Cycle".)
p = [1/200;20/200;1;2];
ap = 2;
eq = p(3)^p(4)*[2;1/2];
x0 = init_EP_EP(@ode_prlgrm_rev,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',0.01);
opt = contset(opt,'Singularities',1);
[x,v,s,h,f] = cont(@equilibrium,x0,[],opt);

x1    = x(1:2,s(2).index);
p(ap) = x(3,s(2).index);
hh    = 1e-6;
ntst  = 20;
ncol  = 4;
[x0,v0] = init_H_LC(@ode_prlgrm_rev, x1, p, ap, hh, ntst, ncol);
opt = contset(opt,'MaxStepSize',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'Adapt',1);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);
disp(['fold bifurcation of limit cycles at b = ',num2str(xlc(end,slc(3).index)*200,'%.2f'),'/200']);
%%
% This is Figure 4.
p = [1/200,22/200,1,2];
eq = [2*p(3)^p(4);1/2*p(3)^p(4);p(3)^2];
c = sum(eq.*[p(4);p(4);2]);

% First we find the stable limit cycle numerically.
odecell = ode_prlgrm_rev;
options = odeset('RelTol',1e-8);
[time,yfwd] = ode45(odecell{2},linspace(0,5000,1e5),eq(1:2)*1.2,...
  options,p(1),p(2),p(3),p(4));
yfwd = [yfwd,(c-p(4)*yfwd(:,1)-p(4)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Here we find a time that is slightly larger than the period.
[time,yfwd] = ode45(odecell{2},linspace(0,250,1e4),...
  yfwd(end,1:2),options,p(1),p(2),p(3),p(4));
yfwd = [yfwd,(c-p(4)*yfwd(:,1)-p(4)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');
x1 = yfwd(:,1);
y1 = yfwd(:,2);
z1 = yfwd(:,3);

% Next we find the unstable limit cycle numerically
% (by going backward in time).
odecell = ode_prlgrm_rev;
options = odeset('RelTol',1e-8);
[time,yfwd] = ode45(odecell{2},linspace(0,-5000,1e5),eq(1:2)*1.06,...
  options,p(1),p(2),p(3),p(4));
yfwd = [yfwd,(c-p(4)*yfwd(:,1)-p(4)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Here we find a time that is slightly larger than the period.
[time,yfwd] = ode45(odecell{2},linspace(0,-90,1e4),...
  yfwd(end,1:2),options,p(1),p(2),p(3),p(4));
yfwd = [yfwd,(c-p(4)*yfwd(:,1)-p(4)*yfwd(:,2))/2];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');
x2 = yfwd(:,1);
y2 = yfwd(:,2);
z2 = yfwd(:,3);

t = linspace(0,p(3),101);
eqx = 2*t.^p(4);
eqy = 1/2*t.^p(4);
eqz = t.^2;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 2
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = 1 : numel(t)
  plot3(x1/p(3)*t(i),y1/p(3)*t(i),z1/p(3)*t(i),'Color',colors{1});
end
for i = 1 : numel(t)
  plot3(x2/p(3)*t(i),y2/p(3)*t(i),z2/p(3)*t(i),'Color',colors{2});
end
plot3(eqx,eqy,eqz,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([-170,20]);
legend({'stable','unstable'},'FontSize',12,'Location','NorthEast');
title(['The limit cycles for \gamma = 2, a = ',num2str(p(1)),', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'2LC_gamma2.jpg']);

%%

% Here gamma = 3.
% This is the fourth plot in Figure 5.

p = [ 1/200 ; 1/20 ; 1.5 ; 3 ];
eq = [2*p(3)^p(4);1/2*p(3)^p(4)];
ap = 3;

x0 = init_EP_EP(@ode_prlgrm_rev,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',  0.0001);
opt = contset(opt,'Singularities',     1);
opt = contset(opt,'Backward',          1);
opt = contset(opt,'MaxNumPoints',    100);
[x1,~,s1] = cont(@equilibrium,x0,[],opt);

% Now we continue the limit cycle that is born at
% the Hopf bifurcation in a tricky way, changing
% the step size frequently.
maxstepsize = [1,1.2:0.2:5,5.5:0.5:8,9:16,...
  18:4:30,35:5:50,60:10:100,120:20:200,300,400,...
  500:500:4000,5000:1000:10000,12000:2000:18000,...
  20000:5000:50000];
maxnumpoints = [100,repmat(5,1,numel(maxstepsize)-1)];

p(ap) = x1(end,s1(2).index);
x = x1(1:2,s1(2).index);
hh = 1e-8;
ntst = 20;
ncol = 4;
opt  = contset;
opt  = contset(opt,'MaxStepSize',  maxstepsize(1));
opt  = contset(opt,'InitStepSize',            0.1);
opt  = contset(opt,'Singularities',             1);
opt  = contset(opt,'Multipliers',               1);
opt  = contset(opt,'Backward',                  1);
opt  = contset(opt,'Adapt',                     1);
opt  = contset(opt,'MaxNumPoints',maxnumpoints(1));

[x0,v0] = init_H_LC(@ode_prlgrm_rev, x, p, ap, hh, ntst, ncol);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

odecell = ode_prlgrm_rev;
options = odeset('RelTol',1e-8);
tolerance = 1e-3;

for i = 2 : numel(maxstepsize)
  p(ap) = xlc(end,end);
  [time,yfwd] = ode45(odecell{2},...
    linspace(0,1.2*xlc(end-1,end),5000),...
    xlc(1:2,end),options,p(1),p(2),p(3),p(4));
  [x0,v0]=initOrbLC(@ode_prlgrm_rev,...
    time,yfwd,p,ap,ntst,ncol,tolerance);
  opt = contset(opt,'MaxStepSize',maxstepsize(i));
  opt = contset(opt,'InitStepSize',maxstepsize(i-1));
  opt = contset(opt,'MaxNumPoints',maxnumpoints(i));
  xlc_new = cont(@limitcycle,x0,v0,opt);
  xlc = [xlc,xlc_new];
end

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
c = p(4)*2*t.^p(4) + p(4)*1/2*t.^p(4) + 2*t.^2;
z = (c-p(4)*x-p(4)*y)/2;

t1 = linspace(0,t(1),1000);
eqx1 = 2*t1.^p(4);
eqy1 = 1/2*t1.^p(4);
eqz1 = t1.^2;
t2 = linspace(t(1),1.1*max(t),1000);
eqx2 = 2*t2.^p(4);
eqy2 = 1/2*t2.^p(4);
eqz2 = t2.^2;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
idx2 = slc(2).index;
idx1 = 1 : idx2-1;
idx3 = idx2+1 : size(xlc,2);
latest = t(1);
for i = idx1(1) : idx1(end)
  if abs(t(i)-latest)>0.0
    plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
    latest = t(i);
  end
end
plot3(x(:,idx2),y(:,idx2),z(:,idx2),'Color',colors{3},'LineWidth',1.5);
latest = c(idx2);
tmp = nan(0);
idxs = nan(1,numel(c));
for i = idx3(1) : idx3(end)
  if c(i)<latest
    plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
    latest = c(i);
    tmp = [tmp,c(i)];
    idxs(i) = i;
  end
end
plot3(eqx1,eqy1,eqz1,'Color',colors{2},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([-170,20]);
legend({'stable','unstable','semistable'},...
  'FontSize',12,'Location','NorthWest');
title(['L_1 > 0 for \gamma = ',num2str(p(4)),...
  ', a = ',num2str(p(1)),', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'2LC_gamma3.jpg']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Section 3.3 Two Andronov-Hopf points %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note: the parameters are (a,b,t).

% Here L1 is negative at both Hopf points.
% This is the first plot in Figure 7.
p   = [ 1/100 ; 1/100 ; 1];
eq  = [p(3) ; p(3)^2/4];
ap  = 3;
x0  = init_EP_EP(@ode_prlgrm_wide,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',  0.01);
opt = contset(opt,'Singularities',   1);
opt = contset(opt,'Backward',        1);
[x1,~,s1]=cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x     = x1(1:2,s1(2).index);
hh    = 1e-6;
ntst  = 80;
ncol  = 4;
opt   = contset;
opt   = contset(opt,'Singularities',1);
opt   = contset(opt,'Multipliers',  1);
opt   = contset(opt,'Backward',     1);
opt   = contset(opt,'Adapt',        1);

[x0,v0] = init_H_LC(@ode_prlgrm_wide, x, p, ap, hh, ntst, ncol);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

idx = slc(3).index;
nphase = 2;
x = xlc(1:nphase:end-2,1:idx);
y = xlc(2:nphase:end-2,1:idx);
t = xlc(end,1:idx);
c = t + 2*1/4*t.^2 + 4*1/4*t.^4;
z = (c-x-2*y)/4;

t1 = linspace(0,0.72,1000);
eqx1 = t1; eqy1 = 1/4*t1.^2; eqz1 = 1/4*t1.^4;
t2 = linspace(min(t),max(t),1000);
eqx2 = t2; eqy2 = 1/4*t2.^2; eqz2 = 1/4*t2.^4;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 2
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = [1:20,21:2:size(x,2)]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
plot3(eqx1,eqy1,eqz1,'Color',colors{1},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{2},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([30,20]);
legend({'stable','unstable'},'FontSize',12,'Location','NorthWest');
title(['L_1^{(1)} < 0, L_1^{(2)} < 0 for a = ',num2str(p(1)),...
  ', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'wide_neg_neg.jpg']);

%%

% Here L1 is negative at the lower Hopf point,
% while it is positive at the upper Hopf point
% This is the second plot in Figure 7.

p = [ 6/1000 ; 3/100 ; 1];
eq  = [p(3) ; p(3)^2/4];
ap  = 3;
x0  = init_EP_EP(@ode_prlgrm_wide,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',0.01);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);
[x1,~,s1]=cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x     = x1(1:2,s1(2).index);
hh    = 1e-6;
ntst  = 80;
ncol  = 4;
opt   = contset;
opt   = contset(opt,'Singularities',1);
opt   = contset(opt,'Multipliers',1);
opt   = contset(opt,'Backward',1);
opt   = contset(opt,'Adapt',1);

[x0,v0] = init_H_LC(@ode_prlgrm_wide, x, p, ap, hh, ntst, ncol);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
c = t + 2*1/4*t.^2 + 4*1/4*t.^4;
z = (c-x-2*y)/4;

t1 = linspace(0,0.72,1000);
eqx1 = t1; eqy1 = 1/4*t1.^2; eqz1 = 1/4*t1.^4;
t2 = linspace(min(t),t(1),1000);
eqx2 = t2; eqy2 = 1/4*t2.^2; eqz2 = 1/4*t2.^4;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = 1 : 1 : slc(3).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
end
i = slc(3).index;
plot3(x(:,i),y(:,i),z(:,i),'Color',colors{3},'LineWidth',1.5);
[~,idx] = min(t);
for i = [slc(3).index+1:2:50,52:3:idx]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
plot3(eqx1,eqy1,eqz1,'Color',colors{1},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{2},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([30,20]);
legend({'stable','unstable','semistable'},'FontSize',12,...
  'Location','NorthWest');
title(['L_1^{(1)} < 0, L_1^{(2)} > 0 for a = ',num2str(p(1)),...
  ', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'wide_neg_pos.jpg']);

%%

% Here L1 is positive at both Hopf points.
% This is the third plot in Figure 7.
p = [ 4/1000 ; 0.046 ; 1];
eq  = [p(3) ; p(3)^2/4];
ap  = 3;
x0  = init_EP_EP(@ode_prlgrm_wide,eq,p,ap);
opt = contset;
opt = contset(opt,'MaxStepSize',0.01);
opt = contset(opt,'MaxNumPoints',200);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Backward',1);
[x1,~,s1] = cont(@equilibrium,x0,[],opt);

p(ap) = x1(end,s1(2).index);
x     = x1(1:2,s1(2).index);
hh    = 1e-6;
ntst  = 80;
ncol  = 4;
opt  = contset;
opt  = contset(opt,'MaxStepSize',0.11);
opt  = contset(opt,'Singularities',1);
opt  = contset(opt,'Multipliers',1);
opt  = contset(opt,'Backward',1);
opt  = contset(opt,'Adapt',1);
opt  = contset(opt,'MaxNumPoints',500);

[x0,v0] = init_H_LC(@ode_prlgrm_wide, x, p, ap, hh, ntst, ncol);
[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
c = t + 2*1/4*t.^2 + 4*1/4*t.^4;
z = (c-x-2*y)/4;

t1 = linspace(0,0.72,1000);
eqx1 = t1; eqy1 = 1/4*t1.^2; eqz1 = 1/4*t1.^4;
t2 = linspace(t(slc(5).index),t(1),1000);
eqx2 = t2; eqy2 = 1/4*t2.^2; eqz2 = 1/4*t2.^4;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = [1:slc(2).index-1, slc(3).index+1:slc(4).index]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
end
for i = slc(2).index+1 : 4 : slc(3).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
for i = [slc(2).index,slc(3).index]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{3},'LineWidth',1.5);
end
plot3(eqx1,eqy1,eqz1,'Color',colors{1},'LineWidth',2);
plot3(eqx2,eqy2,eqz2,'Color',colors{2},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([30,20]);
legend({'stable','unstable','semistable'},'FontSize',12,...
  'Location','NorthWest');
title(['L_1^{(1)} > 0, L_1^{(2)} > 0 for a = ',...
  num2str(p(1)),', b = ',num2str(p(2))],'FontSize',12);
saveas(fig,[dirname,'wide_pos_pos.jpg']);

%%

% This is the critial point: there is a single
% Hopf point and L1 is positive there.
% The cross section is a "figure 8".
% This is the fourth plot in Figure 7.

p = [ 1/256 ; 12/256 ; 1/2];
eq = [p(3) ; p(3)^2/4 ; p(3)^4/4];
c = sum(eq.*[1;2;4]);

% First we find the stable limit cycle numerically.
odecell = ode_prlgrm_wide;
options = odeset('RelTol',1e-8);
[~,yfwd] = ode45(odecell{2},linspace(0,1000,1e5),...
  0.97*eq(1:2)',options,p(1),p(2),p(3));
yfwd = [yfwd,(c-yfwd(:,1)-2*yfwd(:,2))/4];
[time,yfwd] = ode45(odecell{2},linspace(0,1000,1e5),...
  yfwd(end,1:2),options,p(1),p(2),p(3));
yfwd = [yfwd,(c-yfwd(:,1)-2*yfwd(:,2))/4];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Here we find a time that is larger than the period,
% but smaller than twice the period.
[time,yfwd] = ode45(odecell{2},linspace(0,20,1e5),...
  yfwd(end,1:2),options,p(1),p(2),p(3));
yfwd = [yfwd,(c-yfwd(:,1)-2*yfwd(:,2))/4];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Now follow that stable limit cycle as p(3)=t varies.
ntst = 20; ncol = 4; tolerance = 1e-4;
ap = 3;
[x0,v0]=initOrbLC(@ode_prlgrm_wide,...
  time,yfwd(:,1:2),p,ap,ntst,ncol,tolerance);
opt = contset;
opt = contset(opt,'MaxNumPoints',1100);
opt = contset(opt,'MaxStepSize',0.05);
opt = contset(opt,'MinStepSize',0.001);
opt = contset(opt,'InitStepSize',0.05);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'Backward',1);
opt = contset(opt,'Adapt',1);

[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
c = t + 2*1/4*t.^2 + 4*1/4*t.^4;
z = (c-x-2*y)/4;

t = linspace(0,0.72,1000);
eqx = t;
eqy = 1/4*t.^2;
eqz = 1/4*t.^4;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = slc(2).index+1:slc(3).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
end
for i = slc(3).index+1:6:slc(4).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
for i = [slc(2).index,slc(3).index]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{3},'LineWidth',1.5);
end
plot3(eqx,eqy,eqz,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([30,20]);
legend({'stable','unstable','semistable'},'FontSize',12,...
  'Location','NorthWest');
title('L_1 > 0 (t^{(1)}=t^{(2)}=1/2) for a = 1/256, b = 12/256',...
  'FontSize',12);
saveas(fig,[dirname,'wide_figure_8.jpg']);

%%

% Torus.
% This is the fifth plot in Figure 7.

p = [ 1/240 ; 12/256 ; 1/2];
eq = [p(3) ; p(3)^2/4 ; p(3)^4/4];
c = sum(eq.*[1;2;4]);

% First we find the stable limit cycle numerically.
odecell = ode_prlgrm_wide;
options = odeset('RelTol',1e-8);
[time,yfwd] = ode45(odecell{2},linspace(0,1000,1e5),...
  0.97*eq(1:2)',options,p(1),p(2),p(3));
yfwd = [yfwd,(c-yfwd(:,1)-2*yfwd(:,2))/4];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Here we find a time that is larger than the period,
% but smaller than twice the period.
[time,yfwd] = ode45(odecell{2},linspace(0,20,1e5),...
  yfwd(end,1:2),options,p(1),p(2),p(3));
yfwd = [yfwd,(c-yfwd(:,1)-2*yfwd(:,2))/4];
fullscreenFigure;
subplot(2,2,1); plot(yfwd(:,1),yfwd(:,2),'r'); grid on; xlabel('x'); ylabel('y');
subplot(2,2,2); plot(time,yfwd(:,1)); grid on; xlabel('time'); ylabel('x');
subplot(2,2,3); plot(time,yfwd(:,2)); grid on; xlabel('time'); ylabel('y');
subplot(2,2,4); plot(time,yfwd(:,3)); grid on; xlabel('time'); ylabel('z');

% Now follow that stable limit cycle as p(3)=t varies.
ntst = 20; ncol = 4; tolerance = 1e-4;
ap = 3;
[x0,v0]=initOrbLC(@ode_prlgrm_wide,...
  time,yfwd(:,1:2),p,ap,ntst,ncol,tolerance);
opt = contset;
opt = contset(opt,'MaxNumPoints',700);
opt = contset(opt,'MaxStepSize', 0.05);
opt = contset(opt,'MinStepSize', 0.001);
opt = contset(opt,'InitStepSize',0.05);
opt = contset(opt,'Singularities',1);
opt = contset(opt,'Multipliers',1);
opt = contset(opt,'Backward',1);
opt = contset(opt,'Adapt',1);

[xlc,~,slc] = cont(@limitcycle,x0,v0,opt);

nphase = 2;
x = xlc(1:nphase:end-2,:);
y = xlc(2:nphase:end-2,:);
t = xlc(end,:);
c = t + 2*1/4*t.^2 + 4*1/4*t.^4;
z = (c-x-2*y)/4;

t = linspace(0,0.72,1000);
eqx = t;
eqy = 1/4*t.^2;
eqz = 1/4*t.^4;

fig = figure('units','normalized','outerposition',[0,0,0.4,0.6]);
view(3);
hold on;
for i = 1 : 3
  plot(nan,nan,'LineWidth',3,'Color',colors{i});
end
for i = slc(2).index+1:2:slc(3).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{2});
end
for i = slc(3).index+1:4:slc(4).index-1
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{1});
end
for i = [slc(2).index,slc(3).index]
  plot3(x(:,i),y(:,i),z(:,i),'Color',colors{3},'LineWidth',1.5);
end
plot3(eqx,eqy,eqz,'Color',colors{1},'LineWidth',2);
hold off;
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
view([30,20]);
legend({'stable','unstable','semistable'},'FontSize',12,'Location','NorthWest');
title('Every positive equilibrium is stable for a = 1/240, b = 12/256','FontSize',12);
saveas(fig,[dirname,'wide_torus.jpg']);