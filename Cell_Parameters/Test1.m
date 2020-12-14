clear all
disp("Puo' prendere la Rolls, ma me la riporti col pieno.")
% Physics constants
kb=1.38e-23;
% Simulation parameters
sr= 1e5; % Hz, sampling rate of the virtual oscilloscope
dur= 7e-3; % seconds, duration of the simulation
tstart=1e-3;
tstop=5e-3;
Vin_amplitude = 10e-3;
tstep=1/sr;
time=(0:tstep:dur)';
Nstep=size(time,1)-1;

% Cell parameters
disp('Computing signal... ')
Ra= 1e7; %Ohm
Rm=1e8; %Ohm
Cm=3e-11; %F
tau=Cm*(1/Rm + 1/Ra)^(-1);


Vin=zeros(Nstep+1,1);
Vm=zeros(Nstep+1,1);
I0=zeros(Nstep+1,1); % Current without noise
I1=I0; % Current with noise 

Vin=step_fun(Vin,time,tstart,tstop,Vin_amplitude);

noise_amp= 1.5e-10;

f10=figure(10);
f10.Visible='off';
clf
pl1=plot(time*1e3,Vin*1e3,'Linewidth',1.5);
prop = {"Input voltage", "Time [ms]", "Voltage [mV]","m","y"};
SetPlot(get(gcf), prop)
printpdf(f10,'Voltage_input');

for i=1:1:Nstep
	dVm=(-Vm(i)+Vin(i)*Rm/(Rm+Ra))/tau;
	I0(i) = Vm(i)/Rm+Cm*dVm;
	Vm(i+1) = Vm(i)+dVm*tstep;
	I1(i)=I0(i)+noise_amp*randn;
end
f20 = figure(20);
f20.Visible='off';
clf
hold on
grid on
pl21=plot(time*1e3,I1*1e12);
pl2=plot(time*1e3,I0*1e12,'r-','Linewidth',1.5);
prop = {"Current", "Time [ms]"," Current [pA]","m","y"};
SetPlot(get(gcf), prop)
hold on
disp('Saving plot... ')
printpdf(f20,'Current')

% f30 = figure(20);
% f30.Visible='off';
% clf
% pl21=plot(time*1e3,I1*1e12); %, time*10^3, I1*10^12)
% prop = {"Current", "Time [ms]"," Current [pA]","m","y"};
% SetPlot(get(gcf), prop)
% hold on
% disp('Saving plot... ')
% printpdf(f30,'Current_noisy')



%Ra, Rm estimation:
peak1=PeakEstim(time,I1,tstart,5); 
Ra1=Vin_amplitude/peak1;
I1_infty_start = ceil((tstart+5*tau)/tstep);	%I1_infty estimation from data
I1_infty_end = ceil(tstop/tstep)-10;
I1_infty=mean(I1(I1_infty_start:I1_infty_end));
Rm1=Vin_amplitude/I1_infty - Ra1;

%Cm estimation via curve fitting
x= time(time>tstart & time<tstart+13*tau);	%Select right time interval
y=I1(time>tstart & time<tstart+13*tau) - I1_infty; %subtract baseline
fitres=fit(x,y,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
tau1=-1/fitres.b;
Cm1=tau1*(1/Ra1 +1/Rm1);
y1=fitres(x)+I1_infty;
pl3=plot(x,y1,'r-',x,y+I1_infty,'b.');
hold off

% Errors
sigpt=noise_amp;
peak1_sig=sigpt/sqrt(5);
Ra1_sig=Ra*peak1_sig/peak1;
I1_infty_sig=sigpt/sqrt(size(I1(I1_infty_start:I1_infty_end),1));
Rm1_sig=sqrt(Rm^2*(I1_infty_sig/I1_infty)^2 + Ra1_sig^2);
bounds=confint(fitres);
bounds=-(1./bounds);
bounds(:,1)=[];
tau1_sig=(bounds(2)-bounds(1))/2;
R=(1/Ra1 + 1/Rm1);
R_sig=sqrt((Ra1_sig/Ra^2)^2 + (Rm1_sig/Rm1^2)^2);
Cm1_sig=Cm1*sqrt((tau1_sig/tau1)^2 + (R_sig/R)^2);







% --------------- Resting potential --------------------
disp('Computing signal with resting potential... ')
Vrest=-0.08; %mV, resting potential
I0r=I0-I0;
I1r=I0r;
Vmr=zeros(Nstep+1,1);
Vmr(1)=Vrest*Ra/(Rm+Ra);
for i=1:1:Nstep+1
	Vc=(Rm*Vin(i)+Ra*Vrest)/(Ra+Rm);
	dVmr=-1/tau*(Vmr(i)-Vc);
	I0r(i) = ((Vmr(i)-Vrest)/Rm + Cm*dVmr);
	I1r(i)=I0r(i)+noise_amp*randn;
	Vmr(i+1) = Vmr(i)+dVmr*tstep;
end

%Ra, Rm estimation:
peak1=PeakEstim(time,I1r,tstart,5);
I1r_infty_start = ceil((tstart+7*tau)/tstep);	%I1_infty estimation from data
I1r_infty_end = ceil(tstop/tstep)-10;
I1r_infty=mean(I1r(I1r_infty_start:I1r_infty_end));
Ra1r=Vin_amplitude/(peak1-I1r_infty);
Rm1r=Vin_amplitude-Vrest/I1r_infty;


tau1r=0;
Cm1r=0;

f40=figure(40); %plot solution for current
f40.Visible='off';
clf
hold on
grid on
pl21=plot(time*1e3,I1r*1e12);
pl2=plot(time*1e3,I0r*1e12,'g-','Linewidth',1.7);
tit="Current with Vr = " + string(-Vrest*1e3) +" mV";
prop = {tit, "Time [ms]", "Current [pA]","m","y"};
SetPlot(get(gcf), prop)
disp('Saving plot... ')
printpdf(f40,'Current_r')


%Cm estimation via curve fitting
xr= time(time>tstart & time<tstart+13*tau);	%Select right time interval
yr=I1r(time>tstart & time<tstart+13*tau) - I1r_infty; %subtract baseline
fitres1=fit(xr,yr,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
tau1r=-1/fitres.b;
Cm1r=tau1r*(1/Ra1r +1/Rm1r);
y1r=fitres(xr)+I1r_infty;
f50=figure(50); %plot solution for current
f50.Visible='off';
clf
hold on
grid on
pl5=plot(xr*1e3,y1r*1e12,'r-',xr*1e3,yr*1e12+I1r_infty,'b.');
prop = {tit + "fit", "Time [ms]", "Current [pA]","m","y"};
SetPlot(get(gcf), prop)
disp('Saving plot... ')
printpdf(f50,'Current_r_fitted')

% Errors
sigpt=noise_amp;
peak1_sig=sigpt/sqrt(5);
Ra1r_sig=Ra1r*peak1_sig/peak1;
I1r_infty_sig=sigpt/sqrt(size(I1r(I1r_infty_start:I1r_infty_end),1));
Rm1r_sig=sqrt(Rm1r^2*(I1r_infty_sig/I1r_infty)^2 + Ra1r_sig^2);
bounds=confint(fitres1);
bounds=-(1./bounds);
bounds(:,1)=[];
tau1r_sig=(bounds(2)-bounds(1))/2;
R=(1/Ra1r + 1/Rm1r);
R_sig=sqrt((Ra1r_sig/Ra1r^2)^2 + (Rm1r_sig/Rm1r^2)^2);
Cm1r_sig=Cm1r*sqrt((tau1r_sig/tau1r)^2 + (R_sig/R)^2);




% Display some results on screen
disp('-------------- Caso base ---------------------------------')
varnames={'I1_infty [pA]','Ra1 [MOhm]','Rm1 [MOhm]','tau1 [ms]','Cm1 [pf]'};
T=table([I1_infty,I1_infty_sig]*1e12, [Ra1,Ra1_sig]*1e-6, [Rm1,Rm1_sig]*1e-6,[tau1,tau1_sig]*1e3,[Cm1,Cm1_sig]*1e12,'VariableNames',varnames);
disp(T);

disp('-------------- Caso resting potential --------------------')
varnames={'I1r_infty [pA]','Ra1r [MOhm]','Rm1r [MOhm]','tau1r [ms]','Cm1r [pf]'};
Tr=table([I1r_infty,I1r_infty_sig]*1e12, [Ra1r,Ra1r_sig]*1e-6, [Rm1r,Rm1_sig]*1e-6,[tau1r,tau1r_sig]*1e3,[Cm1r,Cm1r_sig]*1e12,'VariableNames',varnames);
disp(Tr);


% ------------------ Bessel filtering ------------------
% Initializing structure fields
npoints=6;

signal(1).time=time;
signal(1).I1=I1;
signal(1).If=I1-I1;
signal(1).base=[0,0];
signal(1).fcutoff=0;
signal(1).Q=0;
signal(1).tau=[0,0];
signal(1).Cm=[0,0];
signal(1).Ra=[0,0];
signal(1).Rm=[0,0];
signal(1).fig=figure('Visible',false);

Qstart=tstart; %Transient charge integration limits
Qend=tstart+3*tau;

fcutoff=logspace(0,1,npoints)*1e3;
fprintf('Filtering %d signals: ', npoints)
for i=1:1:npoints
	signal(i).time=time;
	signal(i).I1=I1;
	signal(i).If=I1-I1;
	signal(i).fcutoff=fcutoff(i);
	%filtering
	[b, a] = besself(4,fcutoff(i)*2*pi);
	[bz, az] = impinvar(b,a,sr);
	signal(i).If=filter(bz,az,signal(1).I1);
	%fitting
	time_peak=time(signal(i).If == max(signal(i).If));
	tstart_fit=time_peak+0.1e-3;
	[base,base_sig]=baseline(signal(i).If,tstart, tstop, tstep,tau);
	signal(i).base(1)=base;
	signal(i).base(2)=base_sig;
	[tauf,tauf_sig,x,y,fitf]=fitthis(signal(i).time,signal(i).If,tstart_fit, tau, signal(i).base(1));
	signal(i).tau=[tauf,tauf_sig];
	%Integrating charge
	Qtransient=signal(i).If(signal(i).time>Qstart & signal(i).time<Qend);
	Qtransient=(Qtransient-signal(i).base(1))*tstep;
	signal(i).Q=sum(Qtransient);
	%Computing parameters
	num=signal(i).tau(1)*Vin_amplitude;
	den=(signal(i).Q + signal(i).tau(1)*signal(i).base(1));
	signal(i).Ra(1)=num/den; % Access resistance
	num_sig=num*tauf_sig/tauf;
	den_sig=sqrt(tauf^2*signal(i).base(2)^2 + signal(i).base(1)^2*tauf_sig^2);
	signal(i).Ra(2)=signal(i).Ra(1)*sqrt((num_sig/num)^2+(den_sig/den)^2);
	
	signal(i).Rm(1)=Vin_amplitude/signal(i).base(1) - signal(i).Ra(1); %Membrane resistance
	signal(i).Rm(2)=sqrt((Vin_amplitude*signal(i).base(2)/signal(i).base(1)^2)^2 + signal(i).Ra(2)^2);

	RR=1/signal(i).Ra(1) +1/signal(i).Rm(1); %Membrane capacity
	signal(i).Cm(1)=tauf*RR;
	RR_sig=sqrt((signal(i).Ra(2)/signal(i).Ra(1)^2)^2 + (signal(i).Rm(2)/signal(i).Rm(1)^2)^2);
	signal(i).Cm(2)=signal(i).Cm(1)*sqrt((tauf_sig/tauf)^2 + (RR_sig/RR)^2);



	signal(i).fig=figure('Visible',false);
	tit='Filtered signal: '+ string(round(fcutoff(i)/1e3,1)) +' kHz';
	prop = {tit, "Time [ms]", "Current [pA]","m","y"};
	clf
	plot(time*1e3,signal(i).If*1e12);
	hold on
	plot(x*1e3,y*1e12,'b.',x*1e3,fitf*1e12,'r-');
	SetPlot(get(gcf), prop);
	grid on
	filename='Bess'+string(ceil(fcutoff(i)))+'Hz';
	printpdf(signal(i).fig,filename);
	fprintf('%d ',i)

	%results
	Cutfreq(i,1)=signal(i).fcutoff;
	Racc(i,:)=signal(i).Ra;
	Rmem(i,:)=signal(i).Rm;
	Cmem(i,:)=signal(i).Cm;

end
	disp(' ')
	varnames={'Cutoff freq.','Ra [MOhm]','Rm [MOhm]','Cm [pF]'};
	T=table(Cutfreq,Racc/1e6,Rmem/1e6, Cmem*1e12,'VariableNames',varnames);
	disp(T)


% Signal plus filtered graphs;
% f60=figure('Visible',false);
f70=figure(70);
f70.Visible='on';
f70.Units = 'centimeters';
f70.OuterPosition = [8 8 30 16];
clf
pl1=plot(time*1e3,signal(1).I1*1e12);
pl1.DisplayName='No filtering';
prop = {"Original vs. Filtered Signal", "Time [ms]", "Current [pA]","m","y"};
SetPlot(get(gcf), prop)
legend()
hold on
grid on
for i=1:npoints
	txt= string(round(signal(i).fcutoff/1e3,1)) + ' kHz cutoff';
	plot(time*1e3,signal(i).If*1e12,'DisplayName',txt,'LineWidth',1.2)
end
hold off
printpdf(f70,'Bess_all_freq')













disp("... Racheeeeeeel!")



function [out_vec] = step_fun(in_vec, time, tstart, tstop, amp)
	in_vec=in_vec-in_vec;
	in_vec(time>=tstart)= amp; %Step function: Vin_amplitude*(theta(tstart)-theta(-tstop))
	in_vec(time>=tstop)=0;
	out_vec=in_vec;
end

function [] = SetPlot(fig, prop )
	font="CMU Sans Serif";
	fontbold= "CMU Sans Serif Bold";
	ax=fig.CurrentAxes;
	ax.Title.String=prop(1);
	ax.Title.FontName=fontbold;
	ax.XLabel.String=prop(2);
	ax.XLabel.FontName=font;
	ax.YLabel.String=prop(3);
	ax.YLabel.FontName=font;
	
	switch prop{4}
	case "s"
		TitFS=18;
		LabFS=12;
	case "m"
		TitFS=22;
		LabFS=18;
	case "l"
		TitFS=26;
		LabFS=20;
	otherwise
		TitFS=20;
		LabFS=15;
	end
	ax.Title.FontSize=TitFS;
	ax.XLabel.FontSize=LabFS;
	ax.YLabel.FontSize=LabFS;

	if(size(prop,2)>=5e-3)
		switch prop{5}
		case "y"
			ax.YLim=enlarge(ax.YLim,0.1);
		case "x"
			ax.XLim=enlarge(ax.XLim,0.1);
		case "xy"
			ax.YLim=enlarge(ax.YLim,0.1);
			ax.XLim=enlarge(ax.XLim,0.1);
		otherwise
			return
		end
	end
end

function [outint]=enlarge(inint,hwmuch)
	interv= abs(inint(2)-inint(1));
	piece= 0.5*(hwmuch*interv);
	if inint(1)<1e-12
		outint= [inint(1), inint(2)+piece];
	else
		outint= [inint(1)-piece, inint(2)+piece];
	end
end

function [peak_val]=PeakEstim(time,signal,peak_time,mean_int)
	%This function performs a mean of the values of
	%the vector signal from time_loc to time_loc+mean_int
	%
	timest=time(2)-time(1);
	peak_idx=ceil(peak_time/timest)+1;
	peak_val=mean(signal(peak_idx:peak_idx+mean_int));
end

function []=printpdf(fig_handle,filename)

	set(fig_handle,'Units','Inches');
	pos = get(fig_handle,'Position');
	set(fig_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
	print(fig_handle,filename,'-dpdf','-bestfit');

end

function [I1_infty, I1_infty_sig]=baseline(I1, tstart, tstop, tstep, tau)
	I1_infty_start = ceil((tstart+5*tau)/tstep);	%I1_infty estimation from data
	I1_infty_end = ceil(tstop/tstep)-10;
	I1_infty = mean(I1(I1_infty_start:I1_infty_end));
	how_many=sqrt(I1_infty_end-I1_infty_start);
	I1_infty_sig=std(I1(I1_infty_start:I1_infty_end))/how_many;
end


function [tau1,tau1_sig,x, y, fitfun]=fitthis(xdata, ydata, tstart_fit, tau, baseline)
	x= xdata(xdata>tstart_fit & xdata<tstart_fit+4*tau);	%Select right time interval
	y= ydata(xdata>tstart_fit & xdata<tstart_fit+4*tau) - baseline; %subtract baseline
	fitres=fit(x,y,'exp1','startpoint',[1e-9,-1/0.2545e-3]); %fit specifying initial parameters
	tau1=-1/fitres.b;
	bounds=confint(fitres);
	bounds=-(1./bounds);
	bounds(:,1)=[];
	tau1_sig=(bounds(2)-bounds(1))/2;
	y=y+baseline;
	fitfun=fitres(x)+baseline;
	% figure(666)
	% pl3=plot(x,y1,'r-',x,y,'b.');

end
