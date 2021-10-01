%solve the initial problem with ERK method. This function may take around 1
%min to provide result. Be patient.
%Change h may results in the bad behavior of the solution.
clear
format long;

V=[1;0;0];%set v_0
t=0;%set t_0=0
h=0.0005;%set time step length h
%This function may take around 1 min to provide result. Be patient pls.
%Change h may results in the bad behavior of the solution.
N=100/h;%get the step number
temp=[];

%4-Stage ERK, straightforward algorithom.
for n=2:1:N+1
    k_1=f(V(:,n-1));%calculate k1
    k_2=f(V(:,n-1)+(h/2)*k_1);
    k_3=f(V(:,n-1)+(h/2)*k_2);
    k_4=f(V(:,n-1)+h*k_3);
    temp=V(:,n-1)+h*(k_1/6+k_2/3+k_3/3+k_4/6);
    V=[V,temp];
end

%plot the result
nexttile
plot(0:h:100,V(1,:)); 
hold on;
plot(0:h:100,V(3,:));
nexttile
plot(0:h:100,V(2,:));

%Define the function for the odes system.
function f_n=f(x)
f_n=[-0.04*x(1)+10^4*x(2)*x(3);
        0.04*x(1)-10^4*x(2)*x(3)-3*10^7*x(2)^2;
        3*10^7*x(2)^2];
end


