%LMM in Exercise 2 for the Lorenz system
clear
format long;
V=[1;1;1];%set v_0
f=[10*(V(2,:)-V(1,:));V(1,:)*(28-V(3,:))-V(2,:);V(1,:)*V(2,:)-(8/3)*V(3,:)];%set f_0
t=0;%set t_0=0
h=0.005;%set time step length h
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
    temp=V(:,n-3)+3*h*f(:,n-3);
    %Evaluate 
    V=[V,temp];%update and save V
    f=[f,[10*(temp(2,:)-temp(1,:));temp(1,:)*(28-temp(3,:))-temp(2,:);temp(1,:)*temp(2,:)-(8/3)*temp(3,:)]];%update and save f  
    t=t+h;%move forward
end

plot3(V(1,:),V(2,:),V(3,:));
xlabel('x');
ylabel('y');
zlabel('z');
hold on
    
    
    