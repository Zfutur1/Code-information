function x = SR(a,b,c,h,x1)

x=zeros(1,length(x1));

for i=1:length(x1)-1
   
    k1=h*(-a*x(i)+((1+a)/b)*x(i).^3-c*x(i)^5+x1(i));
    k2=h*(-a*(x(i)+k1/2)+((1+a)/b)*(x(i)+k1/2).^3-c*(x(i)+k1/2).^5+x1(i));
    k3=h*(-a*(x(i)+k2/2)+((1+a)/b)*(x(i)+k2/2).^3-c*(x(i)+k2/2).^5+x1(i));
    k4=h*(-a*(x(i)+k3)+((1+a)/b)*(x(i)+k3).^3-c*(x(i)+k3).^5+x1(i));
    x(i+1)=x(i)+(1/6)*(k1+2*k2+2*k3+k4);  

end

end