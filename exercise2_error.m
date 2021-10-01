%This file calculate and plot the error of the LMM in exercise2_solver.
%Change step length h in exercise2_solver to produce different results.
%Get the numerical solution
format long
exercise2_solver
%use dsolve to get the exact solution
syms x_1(t) x_2(t)
eqns = [diff(x_1,t)==-2*x_1+x_2+2*sin(t), diff(x_2,t)==998*x_1-999*x_2+999*(cos(t)-sin(t))];
cond = [x_1(0)==2, x_2(0)==3];
[x_1Sol(t),x_2Sol(t)] = dsolve(eqns,cond)
V_exact=[double(x_1Sol(0:h:20));double(x_2Sol(0:h:20))];%t=20

error=abs(V_exact-V);
nexttile
plot(0:h:20,error(1,:))
title('x_1 Global error')
nexttile
plot(0:h:20,error(2,:))
title('x_2 Global error')

