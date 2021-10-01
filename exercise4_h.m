%Plot the h1,h2 by the boundary locus method 
clear
format long;
s=0:0.00001:2*pi;
h_1=((156*exp(2i*s)-12*exp(i*s))+12*exp(i*s).*sqrt(460*exp(3i*s)-611*exp(2i*s)+390*exp(i*s)-99))./(160*exp(i*s)-230*exp(2i*s)-50);
h_2=((156*exp(2i*s)-12*exp(i*s))-12*exp(i*s).*sqrt(460*exp(3i*s)-611*exp(2i*s)+390*exp(i*s)-99))./(160*exp(i*s)-230*exp(2i*s)-50);

%plot h
plot(h_1);
hold on
plot(h_2);


%Numericall apprximate the interval end points, which will be needed for
%boundary locus method analysis in exercise4_r.m
p_1=0;
p_2=0;
for n=1:1:length(h_1);
    %We seek h_1 such that im(h_1) very close to 0 and we discard very similar
    %results as they are helpless.
    if abs(imag(h_1(n))-0)<=10^-6 && abs(real(h_1(n))-p_1)>=10^-2;
        p_1=real(h_1(n))
    end
    %We seek h_2 such that im(h_2) very close to 0 and we discard very similar
    %results as they are helpless.
    if abs(imag(h_2(n))-0)<=10^-6 && abs(real(h_2(n))-p_2)>=10^-2;
        p_2=real(h_2(n))
    end
end

