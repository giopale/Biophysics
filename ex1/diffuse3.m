function [x,r] = diffuse3(D,dt,T,p)

N=ceil(T/dt +.5);

sig = sqrt(2*D*dt/p);

x=zeros(N,3);
delta=zeros(N-1,3);
r=zeros(N,1);
g = rand(N,3);

for i=2:N
    for j=1:3
        if g(i,j)<p
            if g(i,j)<p/2
                x(i,j) = x(i-1,j) - sig;
            else
                x(i,j) = x(i-1,j) + sig;
            end
        else
            x(i,j)=x(i-1,j);
        end      
    end
    r(i)= sum(x(i,:).^2);
end
r=sqrt(r);