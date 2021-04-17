%% 微信公众号：数学模型（MATHmodels）
%  联系方式：添加微信 KingOfModels
clc; clear
figure('position',[50,50,900, 600])

% 绘制车辆和轨迹
h = plotcar(0, 0, 0, 0); hold on
h2 = plot(0,0, ':m', 'linewidth',2);

% 绘制边界
plot([-5, 8, 8,14,14,20,20.0,25.0,25,20,20, 3,3,-5],...
     [-3,-3,10,10,-3,-3, 6.5, 6.5,10,10,16,16,3, 3], '--k', 'linewidth',2)

%% 前轮偏转角度函数 fdelta 和后轴中点速度函数 fv 可由控制函数或手动设置得到
% % [fdelta, fv] = autocontrol(bd, amax, vmax, omega, phimax, dir); 
% 前轮最大转动角速度为 max(dbeta) = 400/16 = 25 deg/s，这里取为 20 deg/s
dbeta = 20;          % [deg/s] 前轮转动角速度，可增大至 25
tmax = 470/16/dbeta; % 前轮最大偏转角为 470/16 度，则转动时间为 470/16/dbeta
% 前轮偏转角度随时间的变化（= 盘转动角度/16）：
fdelta = @(t) 0 + ...
              dbeta*min(t- 1.0,tmax*1).*(t> 1.0) - ...
              dbeta*min(t-10.0,tmax*2).*(t>10.0) + ...
              dbeta*min(t-29.0,tmax*1).*(t>29.0) + ...
              dbeta*min(t-36.0,tmax*1).*(t>36.0) - ...
              dbeta*min(t-66.8,tmax*1).*(t>66.8);

% 后轴中点速度（由油门和刹车控制）
fv = @(t) 1.0 - ...
          0.2*min(t-32.3,6.5).*(t>32.3) + ...
          0.2*min(t-72.0,1.5).*(t>72.0);

% 求解微分方程的时间间隔和范围
dt = 0.01;
t = [0:dt:75]';

% 车身角度 phi，车身中心坐标 (xc,yc) 的初始值 y0 = [phi0, xc0, yc0]
y0 = [0,0,0];

% 求解微分方程数值解，得到任意时刻车身的转角 phi(t) 和车身中心位置 (xc(t), yc(t))
[t,y] = ode45(@odecar, t, y0, [], fv, fdelta);
phi = y(:,1); xc = y(:,2); yc = y(:,3);

% 动态绘图
for i = 1:length(t)
    plotcar(xc(i),yc(i),phi(i),fdelta(t(i)),h);
    set(h2, 'XData', xc(1:i),'YData', yc(1:i));
    drawnow
end

figure
subplot(1,2,1); plot(t,fdelta(t)*16); xlabel('时间(s)'); ylabel('方向盘角度')
subplot(1,2,2); plot(t,fv(t));        xlabel('时间(s)'); ylabel('后车轮速度')

%% ------------------------------------------------------------------------

function dy = odecar(t, y, fv, fdelta)
%% 模型微分方程：输入 y 向量包括车身角度和车身中心位置坐标，输出 dy 是相应的导数
phi = y(1);                        % [deg  ] 车身角度
vx = fv(t);                        % [m/s  ] 沿车身方向的速度
l = 2.8;                           % [m    ] 轴距
delta = fdelta(t);                 % [deg  ] 前轮转角，由方向盘控制
vo = [vx; 0];                      % [m/s  ] 后轴中点速度
omega = vx/l*delta;                % [deg/s] 车身角速度
vg = vo + [0; l/2*deg2rad(omega)]; % [m/s  ] 车身中心速度
vc = [cosd(phi) -sind(phi); sind(phi) cosd(phi)]*vg; % 车身速度变换到真实坐标
dy = [omega; vc];   % dy/dt = [d(phi)/dt; d(xc)/dt; d(yc)/dt]
end

%% ------------------------------------------------------------------------

function h = plotcar(xc, yc, theta, beta, h)
%% 通过变标系旋转和平移变换绘制车身和车胎

% 车长、车宽和轴距
L = 5.0; W = 2.0; l = 2.8;  

% 车身轮廓
x0 = L/2*[1.00, 0.95, 0.85,-0.90,-1.00,-1.00,-0.90, 0.85, 0.95, 1.00]';
y0 = W/2*[0.00, 0.60, 1.00, 1.00, 0.80,-0.80,-1.00,-1.00,-0.60, 0.00]';

% 车胎轮廓
xt = L/10*[1.00, 0.98, 0.95,-0.95,-0.98,-1.00,-1.00,-0.98,-0.95, 0.95, 0.98, 1.00 1.00]';
yt = W/10*[0.60, 0.90, 1.00, 1.00, 0.90, 0.60,-0.60,-0.90,-1.00,-1.00,-0.90,-0.60 0.60]';

% 前车胎：旋转车胎并平移至前车胎位置
[xf, yf] = rotxyd(xt, yt, 0, 0, beta);
xf = xf + [1, 1]*l/2;
yf = yf + [1,-1]*(W/2-W/6);

% 后车胎：平移至前后胎位置
xb = xt - [1, 1]*l/2;
yb = yt + [1,-1]*(W/2-W/6);

% 旋转整车（包括车身、前后车胎）
[x,y]   = rotxyd(x0, y0, 0, 0, theta);
[xf,yf] = rotxyd(xf, yf, 0, 0, theta);
[xb,yb] = rotxyd(xb, yb, 0, 0, theta);

% 平移整车（包括车身、前后车胎）
x = x + xc; y = y + yc;
xb = [xb(:,1); NaN; xb(:,2)] + xc;
yb = [yb(:,1); NaN; yb(:,2)] + yc;
xf = [xf(:,1); NaN; xf(:,2)] + xc;
yf = [yf(:,1); NaN; yf(:,2)] + yc;

% 绘制整车（包括车身、前后车胎）
if nargin<=4
   h = plot(x, y,'k',xb, yb,'b', xf, yf, 'r', 'linewidth',2);
   axis image
   axis([-0.6,5,-0.6,3.2]*L)
else
   set(h(1), 'XData', x , 'YData', y );
   set(h(2), 'XData', xb, 'YData', yb);
   set(h(3), 'XData', xf, 'YData', yf);
end
end

%% ------------------------------------------------------------------------

function [x, y] = rotxyd(x0, y0, xc, yc, deg)
%% 坐标系转换：将点 (x0, y0) 绕 (xc, yc) 旋转 deg 度
x = (x0-xc)*cosd(deg) - (y0-yc)*sind(deg) + xc; 
y = (x0-xc)*sind(deg) + (y0-yc)*cosd(deg) + yc;
end