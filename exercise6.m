%implicit method plus the Newton's method for the Lorenz system
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

%We got [V0,V1,V2], [f0,f1,f2],we then start to get V4... by Newton's method

for n=4:1:N+1
    V_np1=V(:,n-1);%Give an initial guess, we start from a known value. It's a safer choice than starting from the origin.
    e=1; %First set an initial guessing error.
    while e>=0.001 %we choose the tolerance to be 0.001. Peform Newton's method.
        preserve=V_np1;
        %We give the jacobian explicitly in the code.
        J=[1-h*((29/48)*(-10)),-h*((29/48)*(10)),0;
           -h*((29/48)*(28)),1-h*((29/48)*(-1)),h*((29/48)*V_np1(1,:));
           -h*((29/48)*V_np1(2,:)),-h*((29/48)*V_np1(1,:)),1-h*((29/48)*(-(8/3)))];
        %We build the function F explicitly in the code
        F=V_np1-2*V(:,n-1)+(5/4)*V(:,n-2)-(1/4)*V(:,n-3)-h*(-(1/48)*f(:,n-2)-(1/3)*f(:,n-1)+(29/48)*[10*(V_np1(2,:)-V_np1(1,:));V_np1(1,:)*(28-V_np1(3,:))-V_np1(2,:);V_np1(1,:)*V_np1(2,:)-(8/3)*V_np1(3,:)]);
        %We iterate.
        V_np1=V_np1-J\F; %backslash is better then inv(), but do the same thing.
        e=norm(V_np1-preserve)/norm(V_np1);% we calculate the new error and compare it with the tolerance as a stopping criterion
    end
    %after applying the newton method at each time step, we update the final
    %solution.
    f=[f,[10*(V_np1(2,:)-V_np1(1,:));V_np1(1,:)*(28-V_np1(3,:))-V_np1(2,:);V_np1(1,:)*V_np1(2,:)-(8/3)*V_np1(3,:)]];%update and save f
    V=[V,V_np1];%update and save V
    t=t+h;%move forward
end

plot3(V(1,:),V(2,:),V(3,:));
xlabel('x');
ylabel('y');
zlabel('z');

