% To solve the TSP Problem
% Improve version
% auther :alehua
% info  : E = A*sum(Vxi-1)^2+ A*sum(Vxi-1)^2+D*sum(Vxi*dist*Vy,i+1)
%                    x change         i change
% shortcoming E have no upper limit

clear;clc;
N = 10;
times = 3000;
x = 2*rand(1,N)-1;
y = 2*rand(1,N)-1;
citys = [x',y'];
save_citys = citys;
citys = [citys;citys(1,:)];
[r,c] = size(citys);
Init_L = 0;
for i=2:r
    Init_L = Init_L+dist(citys(i-1,:),citys(i,:)');
end
figure(1)
plot(citys(:,1),citys(:,2))
title('Init Route')
hold on 
plot(citys(:,1),citys(:,2),'ro')
for i=1:N
   text(citys(i,1)+0.01,citys(i,2)+0.01,num2str(i));
end
hold off

% step 1
A=1.5;
D=1.5;
u0=0.65;
step=0.01;
% step 2
DistanceCity=dist(save_citys,save_citys');
% step 3  
u=2*rand(N,N)-1;
U=0.5*u0*log(N-1)+u;
V=(1+tanh(U/u0))/2;   %sigmoid function
%V = 1./(1+exp(-u0.*U));
check_V=zeros(N,N);
% step 4
E_data = zeros(1,times);
std_eps = 0.0001;
for k=1:1:times
    % calc dU
    t1=repmat(sum(V,2)-1,1,N);
    t2=repmat(sum(V,1)-1,N,1);
    PermitV1=V(:,2:N);
    PermitV1=[PermitV1 V(:,1)];
    t3=DistanceCity*PermitV1;   %dist*Vy,i+1
    du=-1*(A*t1+A*t2+D*t3);
    U=U+du*step;
    %V = 1./(1+exp(-u0.*U));
    V=(1+tanh(U/u0))/2;
    % calc energy
    t4=sumsqr(sum(V,2)-1);
    t5=sumsqr(sum(V,1)-1);
    PermitV2=V(:,2:N);
    PermitV2=[PermitV2 V(:,1)];
    temp=DistanceCity*PermitV2;
    t6=sum(sum(V.*temp));
    E=0.5*(A*t4+A*t5+D*t6);
    E_data(k) = E;
    if((k>10&&(E_data(k)-E_data(k-1)<std_eps)&&(E_data(k-1)-E_data(k-2)<std_eps)))
        % vilid test
        check_V=zeros(N,N);
        [M,index]=max(V);
        for j=1:N
           check_V(index(j),j)=1;
        end
        C=sum(check_V);
        R=sum(check_V');
        if(~sumsqr(C-R))
            break;
        else
            continue;     
        end
    else
        continue;
    end
end

% result 
% [M,index]=max(V);
if(k==times)
    disp('no suit route')
else
finnal_city = zeros(N,2);
for i=1:N
   finnal_city(i,:)=save_citys(index(i),:);
end
finnal_city = [finnal_city;finnal_city(1,:)];
figure(2)
plot(finnal_city(:,1),finnal_city(:,2))
title('Finnal Route')
hold on 
plot(finnal_city(:,1),finnal_city(:,2),'ro')
for i=1:N
   text(citys(i,1)+0.01,citys(i,2)+0.01,num2str(i));
end
hold off
end
































