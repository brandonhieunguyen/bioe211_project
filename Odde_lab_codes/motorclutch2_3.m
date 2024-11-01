%Benjamin Bangasser
%University of Minnesota
%14 June 2010

%Reduplicate the Chan Motor-Clutch model with event stepping
%Use for simulating multiple runs on several stiffnesses
%Fixed average velocity calculation
% 0 - disengaged clutch
% 1 - engaged clutch

clear
clc

tic

nm=50; %number of myosin motors
Fm=-2; %motor stall force in pN
vu=-120; %unloaded motor velocity in nm/s
nc=50; %number of molecular clutches
kon=0.3; %On rate constant in 1/s
koff=0.1; %Off rate constant in 1/s
Fb=-2; %Bond rupture force in pN
Kc=0.8; %Clutch spring constant in pN/nm
gain=0; %gain of feedback loop

events=1e5; %number of events to simulate
stiffness=logspace(-2,2,10); %define stiffness vector
retro=zeros(1,length(stiffness)); %initialize retro flow vector
failfreq=zeros(1,length(stiffness)); %initialize failure frequency vector
avnumcon=zeros(1,length(stiffness)); %initialize average number of engaged clutches vector
avcough=zeros(1,length(stiffness)); %initialize average koff vector
avtime=zeros(1,length(stiffness)); %initialize average timestep vector
avtrac=zeros(1,length(stiffness)); %initialize average traction force vector
avsubpos=zeros(1,length(stiffness)); %initialize average substrate position vector
avFc=zeros(1,length(stiffness)); %initialize average clutch force vector
binds=zeros(1,length(stiffness)); %initialize number of binding events
failures=zeros(1,length(stiffness)); %initialize number of failure events
cyct=zeros(1,length(stiffness)); %initialize cycle time
velcv=zeros(1,length(stiffness)); %initialize velocity cv
rebcyc=zeros(1,length(stiffness)); %initialize rebinds per cycle

for jj=1:length(stiffness)
    
    Ksub=stiffness(jj) %Substrate spring constant in pN/nm
    
    %initialize vectors for plots
    cstate=zeros(1,nc); %clutch state vector
    cunbind=zeros(1,nc); %clutch unbind state vector
    crebind=zeros(1,nc); %clutch rebind state vector
    xc=zeros(1,nc); %clutch position vector
    t=zeros(1,events+1); %time vector
    subpos=zeros(1,events+1); %substrate position vector
    numcon=zeros(1,events+1); %number of engaged clutches vector
    numcoff=zeros(1,events+1); %number of disengaged clutches vector
    vel=zeros(1,events+1); %velocity vector
    timestep=zeros(1,events+1); %vector of dt's
    cough=zeros(1,events+1); %koff vector
    Ft=zeros(1,events+1); %traction force vector
    totFc=zeros(1,events+1); %mean engaged clutch tension
    
    % Set inital state
    i=1;
    ceng=find(cstate==1); %find indices of engaged clutches
    cdisen=find(cstate==0); %find indices of disengaged clutches
    vf=vu; %claculate actin filament velocity
    xc(ceng)=xc(ceng)+vf*0.005; %claculate positions of engaged clutches(dt=0.005)
    xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate posisiton
    %vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
    %xc(ceng)=xc(ceng)+vf*0.005; %claculate positions of engaged clutches (dt=0.005)
    xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
    Fc=Kc*(xc-xsub); %calculate force on each clutch
    t(i)=0;
    subpos(i)=xsub;
    numcon(i)=length(ceng);
    numcoff(i)=length(cdisen);
    vel(i)=-vf;
    timestep(i)=0;
    
    %Event stepping
    while i<=events
        i=i+1;
        
        %calculate clutch binding times
        if isempty(cdisen)
            tbind=inf;
        else
            tbind=-log(rand(1,length(cdisen)))/kon;
        end
        
        %calculate clutch unbinding times
        if isempty(ceng)
            tunbind=inf;
            cough(i)=koff;
            totFc(i)=0;
        else
            tunbind=-log(rand(1,length(ceng)))./(koff*exp(Fc(ceng)./(Fb+gain*Fc(ceng))));
            cough(i)=mean(koff*exp(Fc(ceng)./(Fb+gain*Fc(ceng))));
            totFc(i)=mean(Fc(ceng));
        end
        
        %find minimum time and execute that event
        [dtbind indbind]=min(tbind);
        [dtunbind indunbind]=min(tunbind);
        if dtbind<dtunbind %free clutch bind to actin
            cstate(cdisen(indbind))=1;
            dt=dtbind;
            binds(jj)=binds(jj)+1;
            if cunbind(cdisen(indbind))==1 %if clutch has already unbound during the cycle
                crebind(cdisen(indbind))=crebind(cdisen(indbind))+1;
            end
        else %engaged clutch disengages from actin
            cstate(ceng(indunbind))=0;
            dt=dtunbind;
            cunbind(ceng(indunbind))=1;
        end
        
        ceng=find(cstate==1); %find indices of engaged clutches
        cdisen=find(cstate==0); %find indices of disengaged clutches
        Ftrac=Ksub*xsub; %claculate traction force
        vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament velocity
        %vf=vu;
        xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches
        xsub=(Kc*sum(xc(ceng))/(Ksub+length(ceng)*Kc)); %calculate substrate posisiton
        %Ftrac=Ksub*xsub; %claculate traction force - old
        %vf=vu*(1-((Ksub*xsub)/(nm*Fm))); %claculate actin filament - old
        %xc(ceng)=xc(ceng)+vf*dt; %claculate positions of engaged clutches - old
        xc(cdisen)=xsub; %calculate posisiton of disengaged clutches
        Fc=Kc*(xc-xsub); %calculate force on each clutch
        
        if xsub==0 %reset unbind vector at failure event
            cunbind=zeros(1,nc);
        end
        
        t(i)=t(i-1)+dt;
        subpos(i)=xsub;
        numcon(i)=length(ceng);
        numcoff(i)=length(cdisen);
        vel(i)=-vf;
        timestep(i)=dt;
        Ft(i)=Ftrac;
    end
    cyctime=diff(t(subpos==0)); %cycle time
    
    retro(jj)=sum((vel.*timestep)/t(events+1)); %weighted average retrograde flowrate
    failfreq(jj)=(length(numcon(numcon==0))-1)/(t(events+1)); %failures/s
    avnumcon(jj)=sum((numcon.*timestep)/t(events+1)); %average number of engaged clutches
    avcough(jj)=sum((cough.*timestep)/t(events+1)); %average koff
    avtime(jj)=mean(timestep); %average timestep
    avtrac(jj)=sum((Ft.*timestep)/t(events+1)); %average traction force
    avsubpos(jj)=sum((subpos.*timestep)/t(events+1)); %average substrate position
    avFc(jj)=sum((totFc.*timestep)/t(events+1)); %average force
    failures(jj)=(length(numcon(numcon==0))-1);
    cyct(jj)=mean(cyctime); %mean cycle time
    velcv(jj)=((sum(timestep.*((vel-retro(jj)).^2)/t(events+1)))^0.5)/retro(jj); %weighted actin flow coeffieicent of variation
    rebcyc(jj)=sum(crebind)./failures(jj);
end
figure;
subplot(1,2,1)
semilogx(stiffness,retro)
xlabel('Substrate stiffness (pN/nm)')
ylabel('Mean retrograde flow (nm/s)')

subplot(1,2,2)
semilogx(stiffness,-avtrac)
xlabel('Substrate stiffness (pN/nm)')
ylabel('Mean traction force (pN)')

toc