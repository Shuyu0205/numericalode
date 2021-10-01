%solve the initial problem with IRK method
clear
format long;
V=[1;0;0];%set v_0
t=0;%set t_0=0
h=0.005;%set time step length h
N=100/h;%get the step number
m=1/4-sqrt(3)/6;%for shorter jacobian
p=1/4+sqrt(3)/6;%for shorter jacobian
q=1/4;%for shorter jacobian
temp=[];
b_1=0.5; %set b_1
b_4=0.5; %set b_1

for n=2:1:N+1
  e=1;%set a initial error, can be arbitrary, better big.
  k_np=[1;0;0;0;0;0];%set an initial guess. It's a safe guess for good reason.
  %Newton's iteration
  while e>=0.00001 %test the tolerance
      preserve=k_np; %keep the old value so that can calculte new tolerance
      %Jacobian
      J=[1+0.04*q*h,-(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*q*h,-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,0.04*m*h,-(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*m*h,-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h;
         -0.04*q*h,1+(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*q*h+3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,-0.04*m*h,(10^4)*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))*m*h+3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h,(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h;
         0,-3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*q*h,1,0,-3*(10^7)*2*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*m*h,0;
         0.04*p*h,-(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*p*h,-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,1+0.04*q*h,-(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*q*h,-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h;
         -0.04*p*h,(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*p*h+3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h  ,(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,-0.04*q*h,1+(10^4)*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))*q*h+3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h,(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h;
         0,-3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*p*h,0,0,-3*(10^7)*2*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*q*h,1];
     %Target function F
      F_k=[k_np(1,:)+0.04*(V(1,n-1)+q*h*k_np(1,:)+m*h*k_np(4,:))-(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:));
           k_np(2,:)-0.04*(V(1,n-1)+q*h*k_np(1,:)+m*h*k_np(4,:))+(10^4)*(V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))*(V(3,n-1)+q*h*k_np(3,:)+m*h*k_np(6,:))+3*(10^7)*((V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))^2);
           k_np(3,:)-3*(10^7)*((V(2,n-1)+q*h*k_np(2,:)+m*h*k_np(5,:))^2);
           k_np(4,:)+0.04*(V(1,n-1)+p*h*k_np(1,:)+q*h*k_np(4,:))-(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:));
           k_np(5,:)-0.04*(V(1,n-1)+p*h*k_np(1,:)+q*h*k_np(4,:))+(10^4)*(V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))*(V(3,n-1)+p*h*k_np(3,:)+q*h*k_np(6,:))+3*(10^7)*((V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))^2);
           k_np(6,:)-3*(10^7)*((V(2,n-1)+p*h*k_np(2,:)+q*h*k_np(5,:))^2)];
      %Newton's method
      k_np=k_np-inv(J)*F_k;
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


%plot the result
nexttile
plot(0:h:100,V(1,:)); 
hold on;
plot(0:h:100,V(3,:));
nexttile
plot(0:h:100,V(2,:));

 

    

