%Benjamin Bangasser
%University of Minnesota
%14 June 2010

%Reduplicate the Chan Motor-Clutch model with event stepping
%Use for simulating one dynamic run
% 0 - disengaged clutch
% 1 - engaged clutch

clear
clc
tic

%Set parameters
nm=50; %number of myosin motors
Fm=-2; %motor stall force in pN
vu=-120; %unloaded motor velocity in nm/s
nc=50; %number of molecular clutches
kon=1; %On rate constant in 1/s
koff=0.1; %Off rate constant in 1/s
Fb=-2; %Bond rupture force in pN
Kc=0.8; %Clutch spring constant in pN/nm
Ksub=1; %Substrate spring constant in pN/nm
gain=0; %gain of feedback loop

eventlimit=1e4; %maximum number of events to simulate

%initialize vectors for plots
cstate=zeros(1,nc); %clutch state vector
xc=zeros(1,nc); %clutch position vector
t=zeros(1,eventlimit); %time vector
subpos=zeros(1,eventlimit); %substrate position vector
numcon=zeros(1,eventlimit); %number of engaged clutches vector
numcoff=zeros(1,eventlimit); %number of disengaged clutches vector
vel=zeros(1,eventlimit); %velocity vector
timestep=zeros(1,eventlimit); %timestep vector
Fcdist=zeros(1,eventlimit);
cdisp1=zeros(1,eventlimit);
ndisp1=zeros(1,eventlimit);

% Set inital state (t=0 s, event #0)
i=1;
dt=0.005;
cstate(1:nc)=0; %set all clutches to disengaged
ceng=find(cstate==1); %find indices of engaged clutches
cdisen=find(cstate==0); %find indices of disengaged clutches
vf=vu; %claculate actin filament velocity
xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches (dt=0.005)
xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate position
%vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
%xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches (dt=0.005)
xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
Fc=Kc*(xc-xsub); %calculate force on each clutch
t(i)=0;
subpos(i)=xsub;
numcon(i)=length(ceng);
numcoff(i)=length(cdisen);
vel(i)=-vf;
timestep(i)=dt;
b=1;
texp(b)=0;
cndisp(b)=50;

%Event stepping
for i=2:eventlimit
    
    Fbe=Fb+gain*Ksub*xsub;
    
    %calculate clutch binding times
    if isempty(cdisen)
        tbind=inf;
    else
        tbind=-log(rand(1,length(cdisen)))/kon;
    end

    %calculate clutch unbinding times
    if isempty(ceng)
        tunbind=inf;
    else
        tunbind=-log(rand(1,length(ceng)))./(koff*exp(Fc(ceng)./(Fb+gain*Fc(ceng))));
    end

    %find minimum time and execute that event
    [dtbind indbind]=min(tbind);
    [dtunbind indunbind]=min(tunbind);
    if dtbind<dtunbind %free clutch bind to actin
        cstate(cdisen(indbind))=1;
        dt=dtbind;
    else %engaged clutch disengages from actin
        cstate(ceng(indunbind))=0;
        dt=dtunbind;
    end

    ceng=find(cstate==1); %find indices of engaged clutches
    cdisen=find(cstate==0); %find indices of disengaged clutches
    vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
    %vf=vu;
    xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches
    xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate posisiton
    %vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
    %xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches
    xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
    Fc=Kc*(xc-xsub); %calculate force on each clutch
    
    t(i)=t(i-1)+dt;
    timestep(i)=dt;
    subpos(i)=xsub;
    numcon(i)=length(ceng);
    numcoff(i)=length(cdisen);
    vel(i)=-vf;
    Fcdist(i)=mean(Fc);
    cdisp1(i)=xc(1);
    ndisp1(i)=xsub;
    if t(i)>2*b
        b=b+1;
        texp(b)=t(i);
        cndisp(b)=ndisp1(i)-cdisp1(i)+50;
    end
        
end

subplot(2,3,1)
plot(t,subpos)
xlabel('time (s)')
ylabel('substrate position (nm)')
axis([0 60 min(subpos) 0])

subplot(2,3,2)
plot(t,numcon,'-b')
hold on
plot(t,numcoff,'-r')
xlabel('time (s)')
ylabel('number of clutches engaged/disengaged')
legend('engaged','disengaged')
axis([0 10 0 nc])

subplot(2,3,3)
plot(t,vel)
xlabel('time (s)')
ylabel('retrograde velocity (nm/s)')
axis([0 10 0 120])

cyctime=diff(t(subpos==0));
subplot(2,3,4)
hist(cyctime,50)
xlabel('failure cycle time (s)')
ylabel('frequency')

[hi lo]=hist(cyctime,50);

for b=1:50
    pdf(b)=(hi(b)/length(cyctime))/(lo(2)-lo(1));
    S(b)=sum(cyctime>lo(b))/length(cyctime);
    h(b)=pdf(b)/S(b);
end
lofrac=lo./mean(cyctime);

subplot(2,3,5)
plotyy(lo,S,lo,pdf)
xlabel('time (s)')
legend('pdf(t)','S(t)')

subplot(2,3,6)
plot(lofrac,h)
xlabel('fraction of mean cycle time')
ylabel('hazard rate (1/s)')

fh = figure(1);
set(fh,'color','white')
toc

%hold off
%plot(texp,cndisp)
%axis([0 60 0 300])
