%PECE method in Exercise 3 for the Lorenz system
clear
format long;
V=[1;1;1];%set v_0
f=[10*(V(2,:)-V(1,:));V(1,:)*(28-V(3,:))-V(2,:);V(1,:)*V(2,:)-(8/3)*V(3,:)];%set f_0
t=0;%set t_0=0
h=0.0025;%set time step length h
N=100/h;%get the step number
temp=[];%set an empty temp vector for later using
t=0;%set t_0=0


%use Euler Method once to find v_1 and f_1
temp=V(:,1)+h*f(:,1);%calculate y_n-1,here we get y_1
V=[V,temp];%y=[v_0,v_1]
t=t+h;%move forward
f=[f,[10*(V(2,2)-V(1,2));V(1,2)*(28-V(3,2))-V(2,2);V(1,2)*V(2,2)-(8/3)*V(3,2)]];%[f0,f1]

%use AB(2) Method to find v_2 and f_2
temp=V(:,2)+(h/2)*(3*f(:,2)-f(:,1));%calculate y_n-1,here we get y_1
V=[V,temp];%y=[v0,v1,v2]
t=t+h;%move forward
f=[f,[10*(V(2,3)-V(1,3));V(1,3)*(28-V(3,3))-V(2,3);V(1,3)*V(2,3)-(8/3)*V(3,3)]];%update and save f

%We got [V0,V1,V2], [f0,f1,f2],we then start to get V4... by PECE method

for n=4:1:N+1
    %Start with a predictor
    V_pre=V(:,n-1)+(h/12)*(23*f(:,n-1)-16*f(:,n-2)+5*f(:,n-3));
    %Evaluate
    f_pre=[10*(V_pre(2,:)-V_pre(1,:));V_pre(1,:)*(28-V_pre(3,:))-V_pre(2,:);V_pre(1,:)*V_pre(2,:)-(8/3)*V_pre(3,:)];
    %Correct
    V_cor=V(:,n-1)+(h/12)*(8*f(:,n-1)+5*f_pre-f(:,n-2));
    %Evaluate
    f_cor=[10*(V_cor(2,:)-V_cor(1,:));V_cor(1,:)*(28-V_cor(3,:))-V_cor(2,:);V_cor(1,:)*V_cor(2,:)-(8/3)*V_cor(3,:)];
    f=[f,f_cor];%update and save f
    V=[V,V_cor];%update and save V
    t=t+h;%move forward
end

plot3(V(1,:),V(2,:),V(3,:));
xlabel('x');
ylabel('y');
zlabel('z');
hold on
    
    
    