%Plot the h1,h2 by the boundary locus method 
clear
format long;
s=0:0.001:2*pi;
h_1=(720*((exp(5i*s))-exp(4i*s)))./(1901*exp(4i*s)-2774*exp(3i*s)+2616*exp(2i*s)-1274*exp(i*s)+251);

%plot h
plot(h_1);

