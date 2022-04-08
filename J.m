function [ f ] = J( x0,Xd,u0,dU,uext,h,Np,Nc )
%% Objective function
X = zeros(5,Np);
X(:,1) = x0;
e = zeros(3,Np);
dU = reshape(dU,2,Nc);
dU = [dU,dU(:,end).*ones(2,Np-Nc)];
U = zeros(4,Np);
utemp = cumsum([u0(1:2,1),dU],2);
U(1:2,:) = utemp(:,2:end);
U(3:4,:) = uext;
for i = 1:Np
    x = X(:,i);
    xd = Xd(:,i);
    x_ = robot(x,U(:,i),h);
    e(:,i) = [ cos(x(3)), sin(x(3)), 0; -sin(x(3)), cos(x(3)), 0; 0, 0, 1]*(xd - x_(1:3));
    X(:,i+1) = x_;
end
Qe = 100*[1,1,1]';
Qdu = 0.1*[1,1]';
Qet = 10*[1,1,1]';
Qdut = 0.1*[1,1]';
f =  sum(sum(Qe.*e.^2)) + sum(sum(Qdu.*dU.^2)) + sum(sum(Qet.*e(:,end).^2)) + sum(sum(Qdut.*dU(:,end).^2));
end

