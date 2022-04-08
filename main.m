clear all
close all
%% simulation parameters 
tspan = 50; % simulation time - [s]
h = 0.1; % sampling intervel - [s]
times = 0:h:tspan; % time sequence - [s]
nmax = length(times);
U = zeros(4,nmax);
taumax = 10; % maximum torqu of wheels - [N*m]
dtaumax = taumax*0.5; % maximum torqu of wheels - [N*m/sec]
Fextmax = 20; % maximum external force - [N]
U(3:4,:) = [Fextmax*sin(0.1*pi*times*0);20*ones(size(times))];
Ur = U;
X = zeros(5,nmax);
% X(:,1) = [0,10,deg2rad(-90),0,0]';
X(:,1) = [0,0,0,0,0]';
%% generate referance trajactory
traj = 2;
switch traj
    case  1 % circle
        tc = 10; % 
        v = 1.2*pi*(1-cos(2*pi*times/tc));
        w = -v/12;
        % generate trajectory
        Xd = zeros(3,nmax);
        for n = 1:nmax-1
            xd = Xd(:,n); 
            Xd(:,n+1) = xd + h* [cos(xd(3)),0;sin(xd(3)),0;0,1]*[v(n+1),w(n+1)]';
        end
    case 2 % '8' shape
        w = zeros(1,nmax);
        tc = 10; 
        v = 1.2*pi*(1-cos(2*pi*times/tc));
        for i = 1:2
            sn = 2*(i-1)*tc/h + 1;
            en = 2*i*tc/h + 1;    
            w(sn:en) = (-1)^i * v(sn:en)/12;
        end
    case 3 % square
        tc = 10;
        v = 20/tc*ones(1,nmax);
        w = zeros(1,nmax);
        for i = 1:4
            sn = 100*i-19;
            en = 100*i+21; 
            if i ==4
                w(sn:en) = -pi/8*(sin(pi/2*(0:0.1:4)-pi/2)+1);
            else
                w(sn:en) = -pi/8*(sin(pi/2*(0:0.1:4)-pi/2)+1);
            end
        end
end
% generate trajectory
Xd = zeros(3,nmax);
for n = 1:nmax-1
    xd = Xd(:,n); 
    Xd(:,n+1) = xd + h* [cos(xd(3)),0;sin(xd(3)),0;0,1]*[v(n+1),w(n+1)]';
end
%	plot referrence trajectory
% figure
% plot(Xd(1,:),Xd(2,:),'b:','Linewidth',2)
% axis equal
%% model predictive control parameters
Np = 4; % predictive horizon
Nc = 4; % control horizon
% using SQP for nonlinear optimization 
options = optimoptions('fmincon','Display','off','Algorithm','sqp'); %,'StepTolerance',1e-4,'MaxIterations',200
%% simulation main loop
runtime = zeros(1,nmax);
for n = 1:nmax-1    
    x0 = X(:,n);
    u0 = U(:,n);
    if nmax - n<= Np
        Np = nmax-n;
        Nc = Np;
    end
    tic
    obj = @(du) J(x0,Xd(:,n:n+Np-1),u0,du,U(3:4,n:n+Np-1),h,Np,Nc);
    lb = repmat(max([-dtaumax*h,-dtaumax*h],-taumax)',Nc,1);
    ub = repmat(min([dtaumax*h,dtaumax*h],taumax)',Nc,1);
    dubst = fmincon(obj,rand(Nc*2,1),[],[],[],[],lb,ub,[],options);
    U(1:2,n+1) = u0(1:2) + dubst(1:2);
    X(:,n+1) = robot(x0,U(:,n+1),h);

    runtime(n) = toc; 
end
%% plots
fsize = [400,0,400,400];
linew = 2;
dX = [zeros(5,1),diff(X,1,2)/h]; % velocity

figure
subplot(2,1,1)
plot(times,X(1,:),'r-','Linewidth',linew);
ylabel('x [m]');
title('Position and velocity of X');
subplot(2,1,2)
plot(times,dX(1,:),'b-','Linewidth',linew);
ylabel('x^{dot} [m/sec]');
xlabel('Time [sec]');
set(gcf,'position',fsize)
saveas(gcf,'plots\x.png')

figure
subplot(2,1,1)
plot(times,X(2,:),'r-','Linewidth',linew);
ylabel('y [m]');
title('Position and velocity of Y');
subplot(2,1,2)
plot(times,dX(2,:),'b-','Linewidth',linew);
ylabel('y^{dot} [m/sec]');
xlabel('Time [sec]');
set(gcf,'position',fsize)
saveas(gcf,'plots\y.png')

figure
subplot(2,1,1)
plot(times,X(3,:),'r-','Linewidth',linew);
ylabel('\theta [rad]');
title('Orientation and angular velocity');
subplot(2,1,2)
plot(times,dX(3,:),'b-','Linewidth',linew);
ylabel('\theta^{dot} [rad/sec]');
xlabel('Time [sec]');
set(gcf,'position',fsize)
saveas(gcf,'plots\theta.png')

figure
subplot(2,1,1)
plot(times,X(4,:),'r-','Linewidth',linew);
ylabel('\omega_{1} [rad/sec]');
title('Angular velocity of wheels');
subplot(2,1,2)
plot(times,X(5,:),'b-','Linewidth',linew);
ylabel('\omega_{2} [rad/sec]');
xlabel('Time [sec]');
set(gcf,'position',fsize)
saveas(gcf,'plots\omega.png')

figure
subplot(2,1,1)
plot(times,dX(4,:),'r-','Linewidth',linew);
ylabel('\omega^{dot}_{1} [rad/sec]');
title('Angular acceleration of wheels');
subplot(2,1,2)
plot(times,dX(5,:),'b-','Linewidth',linew);
ylabel('\omega^{dot}_{2} [rad/sec]');
xlabel('Time [sec]');
set(gcf,'position',fsize)
saveas(gcf,'plots\acc.png')

figure
plot(times,U(3,:),'r-','Linewidth',linew);
hold on
plot(times,U(4,:),'b:','Linewidth',linew);
hold off
legend({'F_{ext,x}','F_{ext,y}'},'Location','best')
title('External disturbances');
set(gcf,'position',[400,300,400,300])
saveas(gcf,'plots\d.png')

figure
plot(Xd(1,:),Xd(2,:),'b:','Linewidth',2*linew);
hold on
plot(X(1,:),X(2,:),'r-','Linewidth',linew);
hold off
axis equal
legend({'reference','actual'},'Location','best')
title('Trajactory tacking');
set(gcf,'position',[400,300,400,300])
saveas(gcf,'plots\trajactory.png')