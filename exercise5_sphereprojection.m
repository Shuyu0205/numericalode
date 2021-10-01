%Solve the Rabinovich¨CFabrikant system on a sphere of radius r = 3 
clear
format long;
V=[-1;0;0.5];%set v_0
t=0;%set t_0=0
h=0.005;%set time step length h
N=50/h;%get the step number
m=1/4-sqrt(3)/6;%for shorter jacobian
p=1/4+sqrt(3)/6;%for shorter jacobian
q=1/4;%for shorter jacobian
temp=[];
b_1=0.5; %set b_1
b_4=0.5; %set b_1

for n=2:1:N+1
  e=1;%set a initial error, can be arbitrary, better big.
  k_np=[-1;0;0;0;0;0];%set an initial guess. It's a safe guess for good reason.
  %Newton's iteration
  while e>=0.00001 %test the tolerance
      preserve=k_np; %keep the old value so that can calculte new tolerance
      %Jacobian
      a_n=V(1,n-1)+q*h*k_np(1,:)+m*h*k_np(4,:);
      b_n=V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:);
      c_n=V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:);
      a_p=V(1,n-1)+p*h*k_np(1,:)+q*h*k_np(4,:);
      b_p=V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:);
      c_p=V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:);
      J=[1-0.5*h*b_n*a_n-(0.87/4)*h,-(c_n-1+(a_n)^2)*(h/4),-(1/4)*h*b_n,-2*m*h*b_n*a_n-0.87*m*h,-(c_n-1+(a_n)^2)*m*h,-m*h*b_n;
         -(3*c_n+1-(a_n)^2)*h*(1/4),1-(0.87/4)*h,-(3/4)*h*a_n,-(3*c_n+1-(a_n)^2)*h*m+2*(a_n)^2*m*h,-0.87*m*h,-3*m*h*a_n;
         0.5*h*c_n*b_n,0.5*h*c_n*a_n,1+0.5*h*(1.1+a_n*b_n),2*m*h*c_n*b_n,2*m*h*c_n*a_n,2*m*h*(1.1+a_n*b_n);
         2*p*h*b_p*a_p-(0.87*p)*h,-(c_p-1+(a_p)^2)*(h*p),-p*h*b_p,1-0.5*h*b_p*a_p-(0.87/4)*h,-(c_p-1+(a_p)^2)*0.25*h,-0.25*h*b_p;
         -(3*c_p+1-(a_p)^2)*h*p,(0.87*p)*h,-(3*p)*h*a_p,-(3*c_p+1-(a_p)^2)*h*0.25+2*(a_p)^2*0.25*h,1-0.87*0.25*h,-3*0.25*h*a_p;
         2*p*h*c_p*b_p,2*p*h*c_p*a_p,2*p*h*(1.1+a_p*b_p),2*0.25*h*c_p*b_p,2*0.25*h*c_p*a_p,1+2*0.25*h*(1.1+a_p*b_p)];
     %Target function F
      F_k=[k_np(1,:)-b_n*(c_n-1+(a_n)^2)-0.87*a_n;
           k_np(2,:)-a_n*(c_n+1-(a_n)^2)-0.87*b_n ;
           k_np(3,:)+2*c_n*(1.1+a_n*b_n) ;
           k_np(4,:)-b_p*(c_p-1+(a_p)^2)-0.87*a_p ;
           k_np(5,:)-a_p*(c_p+1-(a_p)^2)-0.87*b_p ;
           k_np(6,:)+2*c_p*(1.1+a_p*b_p)];
      %Newton's method
      k_np=k_np-pinv(J)*F_k;
      %New Tolerance
      e=norm(k_np-preserve)/norm(k_np);
  end
  %Update result
  k_1=k_np(1:3);%get k_1
  k_4=k_np(4:6);%get k_4
  %Perform the final step of the RK method to get an solution at one point.
  temp=V(:,n-1)+b_1*h.*k_1+b_4*h.*k_4;
  V=[V,temp];% update and n=n+1 then.
end

%Projection method
temp=[];
for n=1:1:length(V)
    L=norm(V(:,n));
    temp=(3.*V(:,n))/L;
    V(:,n)=temp;
end

[x y z] = sphere;
a=[0 0 0 3]
s1=surf(x*a(1,4)+a(1,1),y*a(1,4)+a(1,2),z*a(1,4)+a(1,3));
daspect([1 1 1])
view(30,10)
hold on 
plot3(V(1,:),V(2,:),V(3,:),'-r','linewidth',2);
