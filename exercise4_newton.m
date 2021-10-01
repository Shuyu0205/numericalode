clear
delta_t=1;
delta_x=1;
lamb=delta_t/(delta_x)^2;
e=0.1;
r=1;
L=10;
T=10;
v=-1;
A=[];%Stiffness matrix A
U=[];%Solution matrix U
%build up the initial condition vector;
U_0=[];
A=[0.0913,0,0,0;
    0.2725,0.0913,0,0;
    0.1860,0.3473,0.0193,0;
    0.3286,0.2094,0.3695,0.0913]
%set the initial condition without the boundary
for j=1:1:(L/delta_x)-1
    U_0=[U_0;cos((2*pi*j*delta_x)/L)*sin((10*pi*delta_x*j)/L)];
end
%give the dirichlet boundary condition
U=U_0;%sloution matrix without boundary
U_B=[0;U_0;0];%set boundary

%PECE

for n=1:1:(T/delta_t)
    U_boundary_temp=[];
    for j=2:1:(L/delta_x)
        e=1;
        k_np=[0;1;0;1]
        %Use newton method to get k1,k2,k3,k4
        while e>=0.6 %test the tolerance
              preserve=k_np; %keep the old value so that can calculte new tolerance
              %Jacobian
              J=[1-r*(2*(U_B(j,n)+delta_t*(A(1,1)*k_np(1,:)))*delta_t*A(1,1)*(1-U_B(j,n)-delta_t*(A(1,1)*k_np(1,:)))+((U_B(j,n)+delta_t*(A(1,1)*k_np(1,:)))^2)*(-delta_t*A(1,1))),0,0,0;
                 1-r*(2*(U_B(j,n)+delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))*delta_t*A(2,1)*(1-U_B(j,n)-delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))+((U_B(j,n)+delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))^2)*(-delta_t*A(2,1))),1-r*(2*(U_B(j,n)+delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))*delta_t*A(2,2)*(1-U_B(j,n)-delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))+((U_B(j,n)+delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)))^2)*(-delta_t*A(2,2))),0,0;
                 1-r*(2*(U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))*delta_t*A(3,1)*(1-U_B(j,n)-delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))+((U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))^2)*(-delta_t*A(3,1))),1-r*(2*(U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))*delta_t*A(3,2)*(1-U_B(j,n)-delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))+((U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))^2)*(-delta_t*A(3,2))),1-r*(2*(U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))*delta_t*A(3,3)*(1-U_B(j,n)-delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))+((U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)))^2)*(-delta_t*A(3,3))),0;
                 1-r*(2*(U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))*delta_t*A(4,1)*(1-U_B(j,n)-delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))+((U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))^2)*(-delta_t*A(4,1))),1-r*(2*(U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))*delta_t*A(4,2)*(1-U_B(j,n)-delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))+((U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))^2)*(-delta_t*A(4,2))),1-r*(2*(U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))*delta_t*A(4,3)*(1-U_B(j,n)-delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))+((U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))^2)*(-delta_t*A(4,3))),1-r*(2*(U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))*delta_t*A(4,4)*(1-U_B(j,n)-delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))+((U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))^2)*(-delta_t*A(4,4)))];
             %Target function F
              F_k=[k_np(1,:)+v*((U_B(j+1,n)-U_B(j-1,n))/(2*delta_x))-e*((U_B(j+1,n)-2*U_B(j,n)+U_B(j-1,n))/((delta_x)^2))-r*((U_B(j,n)+delta_t*(A(1,1)*k_np(1,:)+A(1,2)*k_np(2,:)+A(1,3)*k_np(3,:)+A(1,4)*k_np(4,:)))^2)*(1-U_B(j,n)-delta_t*(A(1,1)*k_np(1,:)+A(1,2)*k_np(2,:)+A(1,3)*k_np(3,:)+A(1,4)*k_np(4,:)));
                   k_np(2,:)+v*((U_B(j+1,n)-U_B(j-1,n))/(2*delta_x))-e*((U_B(j+1,n)-2*U_B(j,n)+U_B(j-1,n))/((delta_x)^2))-r*((U_B(j,n)+delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)+A(2,3)*k_np(3,:)+A(2,4)*k_np(4,:)))^2)*(1-U_B(j,n)-delta_t*(A(2,1)*k_np(1,:)+A(2,2)*k_np(2,:)+A(2,3)*k_np(3,:)+A(2,4)*k_np(4,:)));
                   k_np(3,:)+v*((U_B(j+1,n)-U_B(j-1,n))/(2*delta_x))-e*((U_B(j+1,n)-2*U_B(j,n)+U_B(j-1,n))/((delta_x)^2))-r*((U_B(j,n)+delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)+A(3,4)*k_np(4,:)))^2)*(1-U_B(j,n)-delta_t*(A(3,1)*k_np(1,:)+A(3,2)*k_np(2,:)+A(3,3)*k_np(3,:)+A(3,4)*k_np(4,:)));
                   k_np(4,:)+v*((U_B(j+1,n)-U_B(j-1,n))/(2*delta_x))-e*((U_B(j+1,n)-2*U_B(j,n)+U_B(j-1,n))/((delta_x)^2))-r*((U_B(j,n)+delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))^2)*(1-U_B(j,n)-delta_t*(A(4,1)*k_np(1,:)+A(4,2)*k_np(2,:)+A(4,3)*k_np(3,:)+A(4,4)*k_np(4,:)))
                  ];
              %Newton's method
              k_np=k_np-J\F_k;
              %New Tolerance
              e=norm(k_np-preserve)/norm(k_np);
        end
        k_11=k_np(1,:);%get k
        k_22=k_np(2,:); 
        k_33=k_np(3,:);
        k_44=k_np(4,:);
        U_boundary_temp=[U_boundary_temp;U_B(j,n)+delta_t*(0.2202*k_11+0.2692*k_22+0.2886*k_33+0.222*k_44)]; 
    end
    U_boundary_temp=[0;U_boundary_temp;0];
    U_B=[U_B, U_boundary_temp];
end
waterfall(0:delta_x:L,0:delta_t:T,U_B);
