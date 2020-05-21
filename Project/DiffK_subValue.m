clear
clc
tic

%parameters
n_m=75;    %Total number of myosin motors
F_m=-2;    %Single motor stall force, Units: pN
v_u=-120;  %Unloaded motor velocity, Units: nm/s
n_c=75;    %Total number of molecular clutches
k_on=1;    %Pseudo first order on-rate constant, Unit:1/s
k_off=0.1; %Pseudo first order unloaded off-rate constant, Unit:1/s
F_b=-2;    %Characteristic bond rupture force, Unit: pN
K_c=0.8;   %Molecular clutch spring constant, Unit: pN/nm
gain=0;    %Gain of feedback loop
Viscosity=logspace(-5,0,10); %Vector for substrate stiffness
events=1e4; %total events that will be siViscositylated

%Initialization of vectors
avg_retro_1=zeros(1,length(Viscosity));  %Vector to store Retrograde Flow rate
avg_trac_1=zeros(1,length(Viscosity));   %Vector to store mean traction force
avg_retro_2=zeros(1,length(Viscosity));  %Vector to store Retrograde Flow rate
avg_trac_2=zeros(1,length(Viscosity));   %Vector to store mean traction force
avg_retro_3=zeros(1,length(Viscosity));  %Vector to store Retrograde Flow rate
avg_trac_3=zeros(1,length(Viscosity));   %Vector to store mean traction force

%For Ksub= 0.1
for j=1:length(Viscosity)
    K_sub=0.1;                     %Elastic component of the cell
    Mu=Viscosity(j);             %Viscosity Range
    c_state=zeros(1,n_c);        %Clutch state vector
    c_unbind=zeros(1,n_c);       %Clutch unbind state vector
    c_rebind=zeros(1,n_c);       %Clutch rebind state vector
    c_pos=zeros(1,n_c);          %Clutch position vector
    t=zeros(1,events+1);         %Time vector
    sub_pos=zeros(1,events+1);   %Substrate position vector
    n_eng=zeros(1,events+1);     %Number of engaged clutches vector
    n_dis=zeros(1,events+1);     %Number of disengaged clutches vector
    vel=zeros(1,events+1);       %Velocity vector
    timestep=zeros(1,events+1);  %Vector of dt's
    koff_True=zeros(1,events+1); %koff vector from Bell's equation
    F_t=zeros(1,events+1);       %Traction force vector
    tot_Fc=zeros(1,events+1);    %Average engaged clutch tension
    
    i=1; 
    dt=0.005;
    c_eng=find(c_state==1);           %Indices of engaged clutches
    c_disen=find(c_state==0);         %Indices of disengaged clutches
    v_f=v_u;                          %Actin filament velocity
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt; %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %calculate substrate posisiton
    F_sub=0;                          %No engaged clutches
    c_pos(c_disen)=x_sub;             %Position of disengaged clutches
    F_c=K_c*(c_pos-x_sub);            %Force on each clutch
    t(i)=0;
    sub_pos(i)=x_sub;
    n_eng(i)=length(c_eng);
    n_dis(i)=length(c_disen);
    vel(i)=-v_f;
    timestep(i)=0;
    while i<=events
    i=i+1;
    %Time required for clutch binding 
    if isempty(c_disen)
        t_bind=inf;
    else
        t_bind=-log(rand(1,length(c_disen)))/k_on; %Bangasser paper
    end
    
    %Time required for clutch unbinding
    if isempty(c_eng)
        t_unbind=inf;
        koff_true(i)=k_off;
        tot_Fc(i)=0;
    else
        t_unbind=-log(rand(1,length(c_eng)))./(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        koff_true(i)=mean(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        tot_Fc(i)=mean(F_c(c_eng));
    end
    
    [dt_bind, ind_bind]=min(t_bind);       %Minimum time for binding 
    [dt_unbind, ind_unbind]=min(t_unbind); %Minimum time for unbinding 
    
    if dt_bind<dt_unbind %Disengaged clutches engage to actin
        c_state(c_disen(ind_bind))=1;
        dt=dt_bind;
    else %Engaged clutch disengages from actin
        c_state(c_eng(ind_unbind))=0;
        dt=dt_unbind;
    end
    c_eng=find(c_state==1);                %Indices of engaged clutches
    c_disen=find(c_state==0);              %Indices of disengaged clutches
    v_f=v_u*(1-((K_sub*x_sub)/(n_m*F_m)))/(1+(length(c_eng)*Mu*v_u)/(n_m*F_m)); %Actin filament velocity from Chan and Odde 2008 Supp Eqn 3: Substrate position
    F_trac=(K_sub*x_sub) + (length(c_eng)*Mu*v_f);%Traction force
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt;      %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %Substrate posisiton from Chan and Odde 2008 Supp Eqn 5
    %x_sub=(Kc*xi-n_eng*vf*Mu)/(Ks+n_eng*Kc)
    c_pos(c_disen)=x_sub;                  %Position of disengaged clutches= position of substrate
    F_c=K_c*(c_pos-x_sub);                 %Force on each clutch from Chan and Odde 2008 Supp Eqn 2
    
    if x_sub==0 %reset unbind vector at failure event
       c_unbind=zeros(1,n_c);
    end
    t(i)=t(i-1)+dt;
    timestep(i)=dt;
    sub_pos(i)=x_sub;
    n_ceng(i)=length(c_eng);
    n_cdis(i)=length(c_disen);
    vel(i)=-v_f;
    F_t(i)=F_trac;
    cyctime=diff(t(sub_pos==0)); %cycle time
    
    avg_retro_1(j)=sum((vel.*timestep)/t(events+1)); %Average retrograde flow rate
    avg_trac_1(j)=sum((F_t.*timestep)/t(events+1));  %Average traction force
    end
end
%For Ksub= 1
for j=1:length(Viscosity)
    K_sub=1;                     %Elastic component of the cell
    Mu=Viscosity(j);             %Viscosity Range
    c_state=zeros(1,n_c);        %Clutch state vector
    c_unbind=zeros(1,n_c);       %Clutch unbind state vector
    c_rebind=zeros(1,n_c);       %Clutch rebind state vector
    c_pos=zeros(1,n_c);          %Clutch position vector
    t=zeros(1,events+1);         %Time vector
    sub_pos=zeros(1,events+1);   %Substrate position vector
    n_eng=zeros(1,events+1);     %Number of engaged clutches vector
    n_dis=zeros(1,events+1);     %Number of disengaged clutches vector
    vel=zeros(1,events+1);       %Velocity vector
    timestep=zeros(1,events+1);  %Vector of dt's
    koff_True=zeros(1,events+1); %koff vector from Bell's equation
    F_t=zeros(1,events+1);       %Traction force vector
    tot_Fc=zeros(1,events+1);    %Average engaged clutch tension
    
    i=1; 
    dt=0.005;
    c_eng=find(c_state==1);           %Indices of engaged clutches
    c_disen=find(c_state==0);         %Indices of disengaged clutches
    v_f=v_u;                          %Actin filament velocity
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt; %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %calculate substrate posisiton
    F_sub=0;                          %No engaged clutches
    c_pos(c_disen)=x_sub;             %Position of disengaged clutches
    F_c=K_c*(c_pos-x_sub);            %Force on each clutch
    t(i)=0;
    sub_pos(i)=x_sub;
    n_eng(i)=length(c_eng);
    n_dis(i)=length(c_disen);
    vel(i)=-v_f;
    timestep(i)=0;
    while i<=events
    i=i+1;
    %Time required for clutch binding 
    if isempty(c_disen)
        t_bind=inf;
    else
        t_bind=-log(rand(1,length(c_disen)))/k_on; %Bangasser paper
    end
    
    %Time required for clutch unbinding
    if isempty(c_eng)
        t_unbind=inf;
        koff_true(i)=k_off;
        tot_Fc(i)=0;
    else
        t_unbind=-log(rand(1,length(c_eng)))./(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        koff_true(i)=mean(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        tot_Fc(i)=mean(F_c(c_eng));
    end
    
    [dt_bind, ind_bind]=min(t_bind);       %Minimum time for binding 
    [dt_unbind, ind_unbind]=min(t_unbind); %Minimum time for unbinding 
    
    if dt_bind<dt_unbind %Disengaged clutches engage to actin
        c_state(c_disen(ind_bind))=1;
        dt=dt_bind;
    else %Engaged clutch disengages from actin
        c_state(c_eng(ind_unbind))=0;
        dt=dt_unbind;
    end
    c_eng=find(c_state==1);                %Indices of engaged clutches
    c_disen=find(c_state==0);              %Indices of disengaged clutches
    v_f=v_u*(1-((K_sub*x_sub)/(n_m*F_m)))/(1+(length(c_eng)*Mu*v_u)/(n_m*F_m)); %Actin filament velocity from Chan and Odde 2008 Supp Eqn 3: Substrate position
    F_trac=(K_sub*x_sub) + (length(c_eng)*Mu*v_f);%Traction force
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt;      %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %Substrate posisiton from Chan and Odde 2008 Supp Eqn 5
    %x_sub=(Kc*xi-n_eng*vf*Mu)/(Ks+n_eng*Kc)
    c_pos(c_disen)=x_sub;                  %Position of disengaged clutches= position of substrate
    F_c=K_c*(c_pos-x_sub);                 %Force on each clutch from Chan and Odde 2008 Supp Eqn 2
    
    if x_sub==0 %reset unbind vector at failure event
       c_unbind=zeros(1,n_c);
    end
    t(i)=t(i-1)+dt;
    timestep(i)=dt;
    sub_pos(i)=x_sub;
    n_ceng(i)=length(c_eng);
    n_cdis(i)=length(c_disen);
    vel(i)=-v_f;
    F_t(i)=F_trac;
    cyctime=diff(t(sub_pos==0)); %cycle time
    
    avg_retro_2(j)=sum((vel.*timestep)/t(events+1)); %Average retrograde flow rate
    avg_trac_2(j)=sum((F_t.*timestep)/t(events+1));  %Average traction force
    end
end
%For Ksub= 100
for j=1:length(Viscosity)
    K_sub=100;                     %Elastic component of the cell
    Mu=Viscosity(j);             %Viscosity Range
    c_state=zeros(1,n_c);        %Clutch state vector
    c_unbind=zeros(1,n_c);       %Clutch unbind state vector
    c_rebind=zeros(1,n_c);       %Clutch rebind state vector
    c_pos=zeros(1,n_c);          %Clutch position vector
    t=zeros(1,events+1);         %Time vector
    sub_pos=zeros(1,events+1);   %Substrate position vector
    n_eng=zeros(1,events+1);     %Number of engaged clutches vector
    n_dis=zeros(1,events+1);     %Number of disengaged clutches vector
    vel=zeros(1,events+1);       %Velocity vector
    timestep=zeros(1,events+1);  %Vector of dt's
    koff_True=zeros(1,events+1); %koff vector from Bell's equation
    F_t=zeros(1,events+1);       %Traction force vector
    tot_Fc=zeros(1,events+1);    %Average engaged clutch tension
    
    i=1; 
    dt=0.005;
    c_eng=find(c_state==1);           %Indices of engaged clutches
    c_disen=find(c_state==0);         %Indices of disengaged clutches
    v_f=v_u;                          %Actin filament velocity
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt; %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %calculate substrate posisiton
    F_sub=0;                          %No engaged clutches
    c_pos(c_disen)=x_sub;             %Position of disengaged clutches
    F_c=K_c*(c_pos-x_sub);            %Force on each clutch
    t(i)=0;
    sub_pos(i)=x_sub;
    n_eng(i)=length(c_eng);
    n_dis(i)=length(c_disen);
    vel(i)=-v_f;
    timestep(i)=0;
    while i<=events
    i=i+1;
    %Time required for clutch binding 
    if isempty(c_disen)
        t_bind=inf;
    else
        t_bind=-log(rand(1,length(c_disen)))/k_on; %Bangasser paper
    end
    
    %Time required for clutch unbinding
    if isempty(c_eng)
        t_unbind=inf;
        koff_true(i)=k_off;
        tot_Fc(i)=0;
    else
        t_unbind=-log(rand(1,length(c_eng)))./(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        koff_true(i)=mean(k_off*exp(F_c(c_eng)./(F_b+gain*F_c(c_eng))));
        tot_Fc(i)=mean(F_c(c_eng));
    end
    
    [dt_bind, ind_bind]=min(t_bind);       %Minimum time for binding 
    [dt_unbind, ind_unbind]=min(t_unbind); %Minimum time for unbinding 
    
    if dt_bind<dt_unbind %Disengaged clutches engage to actin
        c_state(c_disen(ind_bind))=1;
        dt=dt_bind;
    else %Engaged clutch disengages from actin
        c_state(c_eng(ind_unbind))=0;
        dt=dt_unbind;
    end
    c_eng=find(c_state==1);                %Indices of engaged clutches
    c_disen=find(c_state==0);              %Indices of disengaged clutches
    v_f=v_u*(1-((K_sub*x_sub)/(n_m*F_m)))/(1+(length(c_eng)*Mu*v_u)/(n_m*F_m)); %Actin filament velocity from Chan and Odde 2008 Supp Eqn 3: Substrate position
    F_trac=(K_sub*x_sub) + (length(c_eng)*Mu*v_f);%Traction force
    c_pos(c_eng)=c_pos(c_eng)+v_f*dt;      %Positions of engaged clutches
    x_sub=((K_c*sum(c_pos(c_eng))-(length(c_eng)*v_f*Mu))/(K_sub+length(c_eng)*K_c)); %Substrate posisiton from Chan and Odde 2008 Supp Eqn 5
    %x_sub=(Kc*xi-n_eng*vf*Mu)/(Ks+n_eng*Kc)
    c_pos(c_disen)=x_sub;                  %Position of disengaged clutches= position of substrate
    F_c=K_c*(c_pos-x_sub);                 %Force on each clutch from Chan and Odde 2008 Supp Eqn 2
    
    if x_sub==0 %reset unbind vector at failure event
       c_unbind=zeros(1,n_c);
    end
    t(i)=t(i-1)+dt;
    timestep(i)=dt;
    sub_pos(i)=x_sub;
    n_ceng(i)=length(c_eng);
    n_cdis(i)=length(c_disen);
    vel(i)=-v_f;
    F_t(i)=F_trac;
    cyctime=diff(t(sub_pos==0)); %cycle time
    
    avg_retro_3(j)=sum((vel.*timestep)/t(events+1)); %Average retrograde flow rate
    avg_trac_3(j)=sum((F_t.*timestep)/t(events+1));  %Average traction force
    end
end

figure;
subplot(1,2,1)
semilogx(Viscosity,avg_retro_1)
hold on
semilogx(Viscosity,avg_retro_2)
hold on
semilogx(Viscosity,avg_retro_3)
xlabel('Substrate viscosity (Pa.s)')
ylabel('Mean retrograde flow (nm/s)')
legend('Ksub at 0.1', 'Ksub at 1', 'Ksub at 100')

subplot(1,2,2)
semilogx(Viscosity,-avg_trac_1)
hold on
semilogx(Viscosity,-avg_trac_2)
hold on
semilogx(Viscosity,-avg_trac_3)
xlabel('Substrate viscosity (Pa.s)')
ylabel('Mean traction force (pN)')
legend('Ksub at 0.1', 'Ksub at 1', 'Ksub at 100')
toc