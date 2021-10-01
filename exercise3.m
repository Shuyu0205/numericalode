%Solve pde (1) by finite difference
clear
Global_error=[];
for i=0:1:4  %Choose delta_x=delta_t=0.5, see if the global error converges
    delta_t=(2^(-i))*0.5;
    delta_x=(2^(-i))*0.5;
    lamb=delta_t/(delta_x)^2;
    e=0.1;
    r=1/e;
    L=10;
    T=10;
    A=[];%Stiffness matrix A
    U=[];%Solution matrix U
    %build up the initial condition vector;
    U_0=[];
    U_exact_T=[];
    U_temp=[];
    for j=1:1:(L/delta_x)-1
        U_0=[U_0;0];
    end
    %Set the starting solution
    U=U_0;
    U_boundary=[0;U_0;0];%set boundary
    %Assembly the stiffness matrix A
    for j=1:1:(L/delta_x)+1
        C=[];%a row vector
        for n=1:1:(T/delta_t)+1
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
    for j=1:1:(L/delta_x)+1
        x(j)=(j-1)*delta_x;
    end



    for n=1:1:(T/delta_t)
        %Calculate the solution
        modified=[];
        for j=1:1:(L/delta_x)+1
            modified=[modified;delta_t*(sin(0.2*pi*x(j))+0.2*0.2*pi^2*e*delta_t*(n-1)*sin(0.2*pi*x(j))-r*(1-delta_t*(n-1)*sin(0.2*pi*x(j))))];
        end
        U_temp_boundary=A\(U_boundary(:,n)+(delta_t*r)+modified);
        %Use Euler method to approximate the boundary.
        U_boundary=[U_boundary, U_temp_boundary];%The matrix after assembly the boundary conditions
    end
    for j=1:1:(L/delta_x)+1;
        U_exact_T=[U_exact_T;T*sin(0.2*pi*x(j))];
    end
    Global_error=[Global_error;norm(U_exact_T-U_boundary(:,end))];
end
Global_error
waterfall(0:delta_x:L,0:delta_t:T,U_boundary);
    