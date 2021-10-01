%exercise 1 ODEs by AB(5) method
clear
format long;
V=[1;1];%set v_0
s=10^17;
f=[-1;-1];%set f_0
t=0;%set t_0=0
h=10^(-20);%set time step length h
N=10000;%get the step number
temp=[];%set an empty temp vector for later using

%use Euler Method once to find v_1 and f_1
temp=V(:,1)+h*f(:,1);%calculate y_n-1,here we get v_1
V=[V,temp];%y=[v_0,v_1]
t=t+h;%move forward
f=[f,[-2*V(1,2)+V(2,2);(s-2)*V(1,2)+(1-s)*V(2,2)]];%[f0,f1]

%use AB(2) Method to find v_2 and f_2
temp=V(:,2)+(h/2)*(3*f(:,2)-f(:,1));%v_2
V=[V,temp];%y=[v0,v1,v2]
t=t+h;%move forward
f=[f,[-2*V(1,3)+V(2,3);(s-2)*V(1,3)+(1-s)*V(2,3)]];%[f0,f1,f2]

%AB(3)
temp=V(:,3)+(h/12)*(23*f(:,3)-16*f(:,2)+5*f(:,1));%v_3
V=[V,temp];%V=[v_0,v_1,v_2,v_3]
t=t+h;%move forward
f=[f,[-2*V(1,4)+V(2,4);(s-2)*V(1,4)+(1-s)*V(2,4)]];%[f0,f1,f2,f3]

%AB(4)
temp=V(:,4)+(h/24)*(55*f(:,4)-59*f(:,3)+37*f(:,2)-9*f(:,1));%v_4
V=[V,temp];%V=[v_0,v_1,v_2,v_3,v_4]
t=t+h;%move forward
f=[f,[-2*V(1,5)+V(2,5);(s-2)*V(1,5)+(1-s)*V(2,5)]];%[f0,f1,f2,f3,f4]


for n=6:1:N+1
    temp=V(:,n-1)+(h/720)*(1901*f(:,n-1)-2774*f(:,n-2)+2616*f(:,n-3)-1274*f(:,n-4)+251*f(:,n-5));%calculate v_n-1
    V=[V,temp];%V=[v_0,v_1,...]
    t=t+h;%move forward
    f=[f,[-2*V(1,n)+V(2,n);(s-2)*V(1,n)+(1-s)*V(2,n)]];%calculate f_n-1
end

%Plot option: plot t-x_1,t-x_2 to see time dependence
nexttile
plot(0:h:10000*h,V(1,:))
nexttile
plot(0:h:10000*h,V(2,:))