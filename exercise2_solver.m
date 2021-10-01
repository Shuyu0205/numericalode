%This file use LMM to solve the given system of odes
clear
format long;
V=[2;3];%set v_0
f=[-1;-2];%set f_0
t=0;%set t_0=0
h=0.0005;%set time step length h
N=20/h;%get the step number
temp=[];%set an empty temp vector for later using

%use Euler Method once to find v_1 and f_1
temp=V(:,1)+h*f(:,1);%calculate y_n-1,here we get y_1
V=[V,temp];%y=[v_0,v_1]
t=t+h;%move forward
f=[f,[-2*V(1,2)+V(2,2)+2*sin(t);998*V(1,2)-999*V(2,2)+999*(cos(t)-sin(t))]];%[f0,f1]

%use AB(2) Method to find v_2 and f_2
temp=V(:,2)+(h/2)*(3*f(:,2)-f(:,1));%calculate y_n-1,here we get y_1
V=[V,temp];%y=[v0,v1,v2]
t=t+h;%move forward
f=[f,[-2*V(1,3)+V(2,3)+2*sin(t);998*V(1,3)-999*V(2,3)+999*(cos(t)-sin(t))]];%calculate f_n-1

for n=4:1:N+1
    temp=V(:,n-3)+3*h*f(:,n-3);%calculate v_n-1
    V=[V,temp];%V=[v_0,v_1,...]
    t=t+h;%move forward
    f=[f,[-2*V(1,n)+V(2,n)+2*sin(t);998*V(1,n)-999*V(2,n)+999*(cos(t)-sin(t))]];%calculate f_n-1
end

%Plot option: plot t-x_1,t-x_2 to see time dependence
plot(0:h:20,V(1,:));
xlabel('t');
hold on
plot(0:h:20,V(2,:));
ylabel('x_1,x_2');




