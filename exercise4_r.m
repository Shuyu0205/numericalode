%solve r by chosen h from three region and calculate their length to verify the strict root
%condition
clear
format long;
%h=0
p_1=[1,-1,0,0];
%h=-2
p_2=[1,-73/36,37/18,-25/36];
%h=-1
p_3=[1,-103/144,17/36,-25/144];
Sol_1=roots(p_1)
Length_1=abs(Sol_1)
Sol_2=roots(p_2)
Length_2=abs(Sol_2)
Sol_3=roots(p_3)
Length_3=abs(Sol_3)