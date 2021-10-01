%Solve pde (1) by finite difference
clear
delta_t=0.02;
delta_x=0.02;
lamb=delta_t/(delta_x)^2;
e=0.1;
r=1/e;
L=10;
T=10;
A=[];%Stiffness matrix A
U=[];%Solution matrix U
%build up the initial condition vector;
U_0=[];
for j=1:1:(L/delta_x)-1
    U_0=[U_0;exp((-L^2)*(j*(delta_x)-(L/2))^2)];
end
%Set the starting solution
U=U_0;
U_boundary=[exp(-(L^4)/4);U_0;exp(-(L^4)/4)];%set boundary
%Assembly the stiffness matrix A
for j=1:1:(L/delta_x)-1
    C=[];%a row vector
    for n=1:1:(T/delta_t)-1
        if n==j-1%assign element
            C=[C,-e*lamb];
        elseif n==j%assign element
            C=[C,1+2*e*lamb+delta_t*r];
        elseif n==j+1%assign element
            C=[C,-e*lamb];
        else%assign element
            C=[C,0];
        end
    end%assign row
    A=[A;C];
end

for n=1:1:(T/delta_t)
    %Calculate the solution
    U_temp=A\(U(:,n)+(delta_t*r));
    %Use Euler method to approximate the boundary.
    U_temp_boundary=[U_boundary(1,n)+delta_t*r*(1-U_boundary(1,n));U_temp;U_boundary(1,n)+delta_t*r*(1-U_boundary(1,n))];%add boundary condition
    U=[U,U_temp];
    U_boundary=[U_boundary, U_temp_boundary];%The matrix after assembly the boundary conditions
end
waterfall(0:delta_x:L,0:delta_t:T,U_boundary);
    