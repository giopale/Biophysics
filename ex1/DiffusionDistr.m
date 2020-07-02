%Produces histograms of the distance of diffused particles at fixed times
clear


k=30;
for l=1:k
    p=l*1/k;
Z = 1e5; %number of particles
dt = 1e-5;
T = 1.01e-2;
D=200; %um^2/s diffusion coefficient
N=ceil(T/dt +.5); %number of timesteps

r=zeros(Z,N);
for i=1:Z
    [~,rad]=diffuse3(D,dt,T,p);
    r(i,:)=rad';
end




% xth=linspace(0,12,1000);
edges=0:.08:12;
xhis = edges(1:end-1)+diff(edges)/2;     %  bin ce
xth = xhis;
figure(1)
clf
hold on
for i=[100,200,500,1000]
    [his,~]=histcounts(r(:,i),edges,'Normalization','Pdf');
    plot(xhis,his,'-');
    t = i*dt;
    yth = 1/(sqrt(4*pi())*(D*t)^(3/2)) * exp(-xth.^2/(4*D*t)) .*xth.*xth;
    plot(xth,yth,'-r')
end
chi(l) = sum((his-yth).^2./yth);


dim = [.6 .45 .3 .3];
str = ['p = ',num2str(p),'; Chi = ',num2str(chi(l))];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlabel('Distance from the source [\mu m]')
ylabel('Particle density per 80 nm shell')
hold off
name=['prob',num2str(10*p)];
print(name,'-dpng');
end

%%
figure(2)
plot(2/k:1/k:1,chi(2:end),'-o');
xlabel('p');
ylabel('Chi2 @ 1e3 * dt');
print('Chi2_1000','-dpng');


