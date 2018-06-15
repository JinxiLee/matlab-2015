%% Options

clear all
index_calculation=1;
load_nk=0;
incident=1;     %1: sub, 2: air
N_Shift=[27]';     27%Physically must be integer, 若不是integer那會算錯
Max_Number_of_Loop=100;
%For single wavelength,
%A:1.989262e-04/1.000000e-04
%B:1.913883e-06/1.000000e-04
%C:1.643246e-10/1.000000e-04
%Hopefully this errorhe
%


MSE_A_OK=0.0001;         %mean(delta_A./A)
MSE_B_OK=0.0001;         %mean(delta_B./B)
MSE_C_OK=0.002;         %mean(delta_C./C)
MSE_total_OK=0.005;

Wavelength_Center_for_Thickness_Weighting=540;
Wavelength_BW_for_Thickness_Weighting=50;

Wavelength_Thickness_Detection=0.54;    %micron

Only_Calc_The_First=0;

Lambda=0;

Sample_Path='D:\Users\TuanShu\130512\Y2\sam1\';
Reference_Path='D:\Users\TuanShu\130512\Y2\glass1\';
Spectroscopy_Path='D:\Users\TuanShu\1208030331_5micron 1.jws.txt';
cd('D:\Users\TuanShu\');
%cd(sprintf('%s\\',Spectroscopy_Path));
Data_Spectroscopy=importdata(Spectroscopy_Path);
array=[0];

Frequency_Downsample_Ratio=5;
Lateral_Downsample_Ratio=1;


index_used_reference_plane=100;

lateral_index_calibration_start=50;
lateral_index_calibration_end=700;

lateral_index_reference_start=50;
lateral_index_reference_end=250;


lateral_index_sample_start=650;
lateral_index_sample_end=650;


Ratio_Upper2Reference=0.835;%0.835; %if the front interface has larger interference eff than referential glass, Ratio_Upper2Reference>1
Thickness=[8.2/1.6]'*1E-6;  %[4.2:0.01:4.6]'*1E-6;    %Initial Value
if load_nk ==1
    n_should=dlmread('n.txt')+1i*dlmread('k.txt');
elseif load_nk ==0
    n_should=1.6;
end
Ratio_Lower2Upper=1;     %if the rear interface has larger interference eff than front interface, Ratio_Lower2Upper>1

%%%
DC_cutoff=1;   %works for all lateral position
Least_Separation=2;
Window_Size_Right_1=4;  %use for the left side of 1st interface and right side of second interface
Window_Size_Left_1=5;  %use for the left side of 1st interface and right side of second interface
Window_Size_Right_2=5;  %use for the left side of 1st interface and right side of second interface
Window_Size_Left_2=4.1;  %use for the left side of 1st interface and right side of second interface
Reference_Window_Size=5;
No_Second_Interface_Threshold=0.1;

%%%
Wavelength_Center=540;

Center_Wavelength_micron=0.54;
Wavelength_Considered_Min=510; %510         %nm
Wavelength_Considered_Max=590;  %580

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=N_f*16;


%%% Global arrays generation

c=3E8;

Max_Frequency=c/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=c/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=c/(Wavelength_Center*1E-9);
Frequency_Considered_Min=c/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=c/(Wavelength_Considered_Min*1E-9);             %Hz


if length(array)>1
    cd(Sample_Path);
    Data=importdata('D0.txt');      % Data_2: the glass data 1
elseif exist(Sample_Path,'file')  == 7
    cd(sprintf('%s\\inter\\',Sample_Path));
    Data=importdata('D0.txt');      % Data_2: the glass data 1
else
    Data=importdata(Sample_Path);      % Data_2: the glass data 1
end

Wavelength=Data(:,1);           %nm
Frequency_Old=c./(Wavelength*1E-9);
Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';
Wavelength_micron=(c./Frequency)*1E6;

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');
Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

Frequency_Center_for_Thickness_Weighting=c./(Wavelength_Center_for_Thickness_Weighting*1E-9);
Frequency_BW_for_Thickness_Weighting=c.*(Wavelength_BW_for_Thickness_Weighting*1E-9)./(Wavelength_Center_for_Thickness_Weighting*1E-9).^2;

% Time-domain

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=c*Time;
Position_micron=Position*1E6;

%%% T 



Spectrum_Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);     
Frequency_Spectroscopy=c./(Wavelength_Spectroscopy*1E-9);

T=interp1(Frequency_Spectroscopy,Spectrum_Spectroscopy_Old,Frequency,'spline'); 

%%% Theory - Sample Model (n1 - n - n2)


% n1 = AN100
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


n_AN100=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_AN100=abs(n_AN100);

n_AN100(isnan(n_AN100))=0;

if incident==1
    t_AN100=(2*(1)./(n_AN100+1));
    n1=n_AN100;
    n2=1;
elseif incident==2
    n1=1;
    n2=n_AN100;
    t_AN100=(2*(n_AN100)./(n_AN100+1));
end


%%% To generate the original spectrum

r_AN100=((1-n_AN100)./(n_AN100+1));

%%% Data Loading
Spectrum_Sample_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;

if length(array)>1
    cd(Sample_Path);
    for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
        Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
        Spectrum_Sample_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
        fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
    end
elseif exist(Sample_Path,'file') == 7
    cd(sprintf('%s\\inter\\',Sample_Path));
    Data_Temp_inter=importdata('D0.txt');      % Data_2: the glass data 1
    
    cd(sprintf('%s\\sam\\',Sample_Path));
    Data_Temp_sam=importdata('D0.txt');      % Data_2: the glass data 1
    
    cd(sprintf('%s\\ref\\',Sample_Path));
    Data_Temp_ref=importdata('D0.txt');      % Data_2: the glass data 1
    
    Data_Temp=Data_Temp_inter-Data_Temp_sam-Data_Temp_ref;
    Spectrum_Sample_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
else

    Data_Temp=importdata(Sample_Path);      % Data_2: the glass data 1
    Spectrum_Sample_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
end

Spectrum_Reference_Frequency(1:length(Wavelength),1:fix(length(array)/Lateral_Downsample_Ratio))=0;
if length(array)>1
    cd(Reference_Path);
    for array_number=1:fix(length(array)/Lateral_Downsample_Ratio)
        Data_Temp=importdata(sprintf('D%i.txt',array(1+(array_number-1)*Lateral_Downsample_Ratio)));      % Data_2: the glass data 1
        Spectrum_Reference_Frequency(:,array_number)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
        fprintf('Loading B-Scan Data...%d/%d\n',array_number,length(array));
    end
elseif exist(Sample_Path,'file') == 7
    cd(sprintf('%s\\inter\\',Reference_Path));
    Data_Temp_inter=importdata('D0.txt');      % Data_2: the glass data 1
    
    cd(sprintf('%s\\sam\\',Reference_Path));
    Data_Temp_sam=importdata('D0.txt');      % Data_2: the glass data 1
    
    cd(sprintf('%s\\ref\\',Reference_Path));
    Data_Temp_ref=importdata('D0.txt');      % Data_2: the glass data 1
    
    Data_Temp=Data_Temp_inter-Data_Temp_sam-Data_Temp_ref;
    Spectrum_Reference_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
else
    Data_Temp=importdata(Reference_Path);
    Spectrum_Reference_Frequency(:,1)=((Data_Temp(:,2)-Data_Temp(end,2)).*((Wavelength*1E-9).^2)/c);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/c);;
end

%plot(Wavelength,Spectrum_Sample_Wavelength);
%xlabel('Wavelength (micron)');
%ylabel('Spectral Power (a.u.)');

%%%

Spectrum_Sample=interp1(Frequency_Old,Spectrum_Sample_Frequency,Frequency); 
Spectrum_Sample(isnan(Spectrum_Sample))=0;
Spectrum_Sample(Frequency<Min_Frequency,:)=0;
Spectrum_Sample(N_f+1:N_t,:)=0;
Signal_Sample=fft(Spectrum_Sample,[],1);  
Signal_Sample(round(size(Signal_Sample,1)/2+1):end,:)=0;
Spectrum_Sample=Spectrum_Sample(1:N_f,:);


Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency); 
Spectrum_Reference(isnan(Spectrum_Reference))=0;
Spectrum_Reference(Frequency<Min_Frequency,:)=0;
Spectrum_Reference(N_f+1:N_t,:)=0;
Signal_Reference=fft(Spectrum_Reference,[],1);  
Signal_Reference(round(size(Signal_Reference,1)/2+1):end,:)=0;
Spectrum_Reference=Spectrum_Reference(1:N_f,:);


clear Spectrum_Sample_Frequency Spectrum_Reference_Frequency

plot(Position_micron,Signal_Sample,Position_micron,Signal_Reference);

plot(Wavelength_micron(Wavelength_micron<1),Spectrum_Sample(Wavelength_micron<1));
xlabel('Wavelength (micron)');
ylabel('Spectral Intesity (a.u.)');
xlim([0.45 0.7]);

plot(Position_micron,Signal_Sample);
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');
xlim([5 35]);

%%% Spectrum generation
DC_cutoff_index=find(Position_micron>DC_cutoff,1,'first');
Least_Separation_index=find(Position_micron>Least_Separation,1,'first');
Window_Size_Right_1_index=find(Position_micron>Window_Size_Right_1,1,'first');
Window_Size_Left_1_index=find(Position_micron>Window_Size_Left_1,1,'first');
Window_Size_Right_2_index=find(Position_micron>Window_Size_Right_2,1,'first');
Window_Size_Left_2_index=find(Position_micron>Window_Size_Left_2,1,'first');
Reference_Window_Size_index=find(Position_micron>Reference_Window_Size,1,'first');
Signal_Sample(1:DC_cutoff_index,:)=0;
Signal_Reference(1:DC_cutoff_index,:)=0;
[Max_Value Max_index]=max(abs(Signal_Sample));
[Max_Value_Ref Max_index_Ref]=max(abs(Signal_Reference));
Signal_Sample_First=Signal_Sample;
Signal_Sample_Second=Signal_Sample;

for p=1:length(Max_index)
    Signal_Sample_Second((Max_index(p)-Least_Separation_index):(Max_index(p)+Least_Separation_index),p)=0;
end

[Max_Value_Second Max_index_Second]=max(abs(Signal_Sample_Second));
If_Second_Exist=Max_Value_Second>Max_Value*No_Second_Interface_Threshold;
If_Second_Lower=Max_index_Second>Max_index;                                    %1 if second is lower, 0 if second is upper

if If_Second_Lower ==1
    [minvalue minindex]=min(abs(Signal_Sample_Second((Max_index+Least_Separation_index+1):Max_index_Second)));
    Separation_index=minindex+Max_index+Least_Separation_index;
elseif If_Second_Lower ==0
    [minvalue minindex]=min(abs(Signal_Sample_Second(Max_index_Second:(Max_index-Least_Separation_index-1))));
    Separation_index=minindex+Max_index_Second-1;
end
    
round((Max_index+Max_index_Second)/2);

First_No_Second_Index=find(If_Second_Exist==0,1,'first');
Last_With_Second_Index=find(If_Second_Exist==1,1,'last');

Safe_Reference_Index=round((First_No_Second_Index+length(If_Second_Exist))/2);
Safe_Sample_Index=round((Last_With_Second_Index+1)/2);

if (Last_With_Second_Index>First_No_Second_Index)
    disp('Second Peak Condition Working Badly!');
else
    fprintf('Film caoted from index: 1 to %d.\n',Last_With_Second_Index);
    if xor(abs(Max_index(Safe_Reference_Index)-Max_index(Safe_Sample_Index))>abs(Max_index(Safe_Reference_Index)-Max_index_Second(Safe_Sample_Index)),mean(If_Second_Lower(1:Safe_Sample_Index))>0.5)
        disp('Film facing down.');
    else
        disp('Film facing up.');
    end
        
end
if Only_Calc_The_First == 1
    Active_Array_Size=1;
else
    Active_Array_Size=length(array(1:Last_With_Second_Index));
end

Signal_Sample_1(1:N_t,1:length(array))=0;
Signal_Sample_2(1:N_t,1:length(array))=0;

for p=1:length(Max_index)
    if If_Second_Exist(p)==0 %Second non-exist
        Signal_Sample_First(1:(Max_index(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_Right_1_index):end,p)=0;
        Signal_Sample_Second(:,p)=0;
        Max_index_Second(p)=Max_index(p); %to avoid wrong calc of thickness (OPD)
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        
    elseif If_Second_Lower(p)==0 %Second exist and second is upper 
        Signal_Sample_First(1:max(Separation_index(p),Max_index(p)-Window_Size_Left_2_index),p)=0;
        %Signal_Sample_First(1:Separation_index(p),p)=0;
        Signal_Sample_First((Max_index(p)+Window_Size_Right_2_index):end,p)=0;
        Signal_Sample_Second(1:(Max_index_Second(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_Second(min(Separation_index(p),Max_index_Second(p)+Window_Size_Right_1_index):end,p)=0;
        %Signal_Sample_Second(Separation_index(p):end,p)=0;
        Signal_Sample_1(:,p)=Signal_Sample_Second(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_First(:,p);

    else %Second exist and second is lower        
        Signal_Sample_First(min(Separation_index(p),Max_index(p)+Window_Size_Right_1_index):end,p)=0;
        %Signal_Sample_First(Separation_index(p):end,p)=0;
        Signal_Sample_First(1:(Max_index(p)-Window_Size_Left_1_index),p)=0;
        Signal_Sample_Second((Max_index_Second(p)+Window_Size_Right_2_index):end,p)=0;
        Signal_Sample_Second(1:max(Separation_index(p),(Max_index_Second(p)-Window_Size_Left_2_index)),p)=0;
        %Signal_Sample_Second(1:Separation_index(p),p)=0;
        
        Signal_Sample_1(:,p)=Signal_Sample_First(:,p);
        Signal_Sample_2(:,p)=Signal_Sample_Second(:,p);
        
    end
    Signal_Reference(1:(Max_index_Ref(p)-Reference_Window_Size_index),p)=0;
    Signal_Reference((Max_index_Ref(p)+Reference_Window_Size_index):end,p)=0;
end

clear Signal_Sample
%clear Signal_Sample_First Signal_Sample_Second


Spectrum_Sample_1=ifft(Signal_Sample_1,[],1);
Spectrum_Sample_1=2*Spectrum_Sample_1(1:N_f,:);                

Spectrum_Sample_2=ifft(Signal_Sample_2,[],1);
Spectrum_Sample_2=2*Spectrum_Sample_2(1:N_f,:);


Spectrum_Reference=ifft(Signal_Reference,[],1);
Spectrum_Reference=2*Spectrum_Reference(1:N_f,:);                



plot(real(Signal_Sample_1));
 

Spectrum_Sample_1=Spectrum_Sample_1(:,1:Active_Array_Size);
Spectrum_Sample_2=Spectrum_Sample_2(:,1:Active_Array_Size);   %!!!!!!!!!! 因為接下來就要算n k 了, 把沒有兩個介面的放進來情況會很怪

Spectrum_Devided=Spectrum_Sample_1./repmat(Spectrum_Reference,1,Active_Array_Size);
Spectrum_Devided_2=Spectrum_Sample_2./repmat(Spectrum_Reference,1,Active_Array_Size);
Spectrum_Devided_3=Spectrum_Sample_2./Spectrum_Sample_1;

subplot(1,2,1);

plot(Position_micron(Position_micron<55),Signal_Sample_1(Position_micron<55),Position_micron(Position_micron<55),Signal_Sample_2(Position_micron<55),Position_micron(Position_micron<55),Signal_Reference(Position_micron<55));
xlabel('Optical Path Difference (micron)');
ylabel('Amplitude (a.u.)');

%plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Devided_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),(abs(Spectrum_Devided_2(Wavelength_micron<1)./abs(Spectrum_Devided(Wavelength_micron<1)))));
%xlabel('Wavelength (micron)');
%ylabel('Interference Power spectral density (a.u.)');
subplot(1,2,2);
plot(Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_1(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Sample_2(Wavelength_micron<1)),Wavelength_micron(Wavelength_micron<1),abs(Spectrum_Reference(Wavelength_micron<1)));
xlabel('Wavelength (micron)');
ylabel('Interference Power spectral density (a.u.)');
fprintf('Peak height ratio:%f\n',Max_Value/Max_Value_Second);

Thickness_Profile=abs((Max_index_Second-Max_index)*(Position_micron(2)-Position_micron(1)));

%plot((1:size(Thickness_Profile,2))*5,Thickness_Profile(1,:));
%xlabel('Lateral Position (micron)');
%ylabel('Thickness (OPD, micron)');

%%


Frequency_Considered=repmat(Frequency(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
Frequency_Considered=downsample(Frequency_Considered,Frequency_Downsample_Ratio);
Aexp=abs(Spectrum_Sample_1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:));  
Aexp=downsample(Aexp,Frequency_Downsample_Ratio);
Bexp=abs(Spectrum_Sample_2(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:));  
Bexp=downsample(Bexp,Frequency_Downsample_Ratio);
Cunwrap=angle(Spectrum_Devided_3(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index,:)); 
Cunwrap=downsample(Cunwrap,Frequency_Downsample_Ratio);
Cunwrap=(unwrap(Cunwrap,[],1));

Wavelength_micron_Considered=Wavelength_micron(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Wavelength_micron_Considered=downsample(Wavelength_micron_Considered,Frequency_Downsample_Ratio);
Center_Wavelength_micron_Index=find(Wavelength_micron_Considered<Center_Wavelength_micron,1,'first');  

if length(Frequency_Considered)>1
    Cunwrap_intercept=interp1(Frequency_Considered,Cunwrap,0,'linear','extrap');
else
    Cunwrap_intercept=-175;
end
%FitTypeLinear=fittype( @(a, b, x) a+b*x, 'independent', {'x'},'dependent', 'z');    %(x-A)^2+(y-B)^2+(z-C)^2=R^2
%Fit_C=fit(Frequency_Considered,Cexp,FitTypeLinear);
%C_intercept=Fit_C.a;               %很奇怪, 如果用這段, 每次算的結果都會不一樣, 可能是clear
                                    %all不會clear fittype之類的, 我乾脆用interp好了

Cunwrap=Cunwrap-2*pi*fix(Cunwrap(Center_Wavelength_micron_Index)/2/pi);%-C_intercept;
%Cexp(size(Cexp,1)+1,:)=2*Cexp(size(Cexp,1),:)-Cexp(size(Cexp,1)-1,:);

%clear Spectrum_Devided Spectrum_Devided_2 Spectrum_Devided_3
Spectrum_Reference_Considered=Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index);
Spectrum_Reference_Considered=downsample(Spectrum_Reference_Considered,Frequency_Downsample_Ratio);

T_Considered=repmat(T(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
T_Considered=downsample(T_Considered,Frequency_Downsample_Ratio);
%Weight_Function=(abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)).*T_Considered)./max((abs(Spectrum_Reference(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index).*T_Considered)));

if incident==1
    n1_Considered=repmat(n1(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
    n1_Considered=downsample(n1_Considered,Frequency_Downsample_Ratio);
    n2_Considered=repmat(n2,size(n1_Considered,1),Active_Array_Size);
elseif incident==2
    n2_Considered=repmat(n2(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
    n2_Considered=downsample(n2_Considered,Frequency_Downsample_Ratio);
    n1_Considered=repmat(n1,size(n2_Considered,1),Active_Array_Size);
end


r_AN100_Considered=repmat(r_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
r_AN100_Considered=downsample(r_AN100_Considered,Frequency_Downsample_Ratio);
t_AN100_Considered=repmat(t_AN100(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index),1,Active_Array_Size);
t_AN100_Considered=downsample(t_AN100_Considered,Frequency_Downsample_Ratio);


Weighting_Function_for_Thickness=gaussmf(Frequency_Considered,[Frequency_BW_for_Thickness_Weighting Frequency_Center_for_Thickness_Weighting]);
Weighting_Function_for_Thickness=Weighting_Function_for_Thickness/sum(Weighting_Function_for_Thickness);

plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,Bexp);
xlabel('Wavelength (micron)');
ylabel('Normalized Inteference Spectrum');

plot(Wavelength_micron_Considered,Cunwrap);
xlabel('Wavelength (micron)');
ylabel('Continuous Phase Spectrum');
cd('D:\Users\TuanShu\');
[max_Signal_Sample_2 maxindex_Signal_Sample_2]=max(abs(Signal_Sample_2));
[max_Signal_Sample_1 maxindex_Signal_Sample_1]=max(abs(Signal_Sample_1));
%dlmwrite('Position_micron.txt',Position_micron-Position_micron(maxindex_Signal_Sample_1),'delimiter','\t','newline','pc');
%dlmwrite('Signal_Sample_1.txt',real(Signal_Sample_1)/max_Signal_Sample_2,'delimiter','\t','newline','pc');
%dlmwrite('Signal_Sample_2.txt',real(Signal_Sample_2)/max_Signal_Sample_2,'delimiter','\t','newline','pc');
max_Aexp=max(Aexp);

Thickness_Profile=repmat(Thickness_Profile(1:Active_Array_Size),length(Frequency_Considered),1);
%%
dn=1E-10;
dt=1E-16;
df=10;       


t1= @ (n) (2*n1_Considered)./(n1_Considered+n); 
t1_r= @ (n) (2*n)./(n1_Considered+n);
t2= @ (n) (2*n)./(n+n2_Considered); 
r1= @ (n) (n1_Considered-n)./(n1_Considered+n);
r1_r= @ (n) (n-n1_Considered)./(n+n1_Considered);
r2= @ (n) (n-n2_Considered)./(n+n2_Considered);

A = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(r1(n))).*repmat(abs(Spectrum_Reference_Considered),1,Active_Array_Size);
B = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) abs(Ratio_Upper2Reference_Now./r_AN100_Considered.*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))).*repmat(abs(Spectrum_Reference_Considered),1,Active_Array_Size);
C = @ (n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now) angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n))))+(4*pi.*Frequency_Considered.*real(n).*Thickness_Now/c);
%.*exp(1i*4*pi.*(Frequency_Considered+df).*n.*Thickness_Now/c)
%-unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c)))),[],1); %雖然沒有絕對的phase, 但是是不是應該只shift pi的整數?
                   
det_array = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9) Q1.*Q5.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 -Q3.*Q5.*Q7 -Q2.*Q4.*Q9 -Q1.*Q6.*Q8;

%For LMA
Q1_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q1 + Lambda*Q1.*(Q5.*Q9 - Q6.*Q8)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q2_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q2 - Lambda*Q2.*(Q4.*Q9 - Q6.*Q7)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q3_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q3 + Lambda*Q3.*(Q4.*Q8 - Q5.*Q7)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q4_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q4 - Lambda*Q1.*(Q2.*Q9 - Q3.*Q8)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q5_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q5 + Lambda*Q2.*(Q1.*Q9 - Q3.*Q7)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q6_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q6 - Lambda*Q3.*(Q1.*Q8 - Q2.*Q7)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q7_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q7 + Lambda*Q1.*(Q2.*Q6 - Q3.*Q5)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q8_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q8 - Lambda*Q2.*(Q1.*Q6 - Q3.*Q4)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);
Q9_P_add_Lambda_PT_inv = @ (Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Lambda) Q9 + Lambda*Q3.*(Q1.*Q5 - Q2.*Q4)./(Q1.*Q5.*Q9 - Q1.*Q6.*Q8 - Q2.*Q4.*Q9 + Q2.*Q6.*Q7 + Q3.*Q4.*Q8 - Q3.*Q5.*Q7);




%T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(1i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(1i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from 分子4*pi > 2*pi
T_abs = @ (n,Thickness_Now) abs(t_AN100_Considered.*t1(n).*t2(n).*exp(1i*2*pi*Frequency_Considered.*n.*Thickness_Now/c)./(1-r1_r(n).*r2(n).*exp(1i*4*pi*Frequency_Considered.*n.*Thickness_Now/c))).^2;          %Note, 120827 remove 2 from 分子4*pi > 2*pi                    

if length(n_should)==1
    n(:,1:Active_Array_Size)=n_should;
else
    n(:,1:Active_Array_Size)=repmat(n_should,1,Active_Array_Size);
end

dn_real_record=0;
dn_imag_record=0;
dthickness_record=0;
%% Start the n k fitting    
Current_Progress=0;

if length(Thickness)==1 %if length(Thickness)>1 > Uniqueness test
    MSE_A_Record(1:length(N_Shift))=0;
    MSE_B_Record(1:length(N_Shift))=0;
    MSE_C_Record(1:length(N_Shift))=0;
else
    MSE_A_Record(1:length(Thickness))=0;
    MSE_B_Record(1:length(Thickness))=0;
    MSE_C_Record(1:length(Thickness))=0;
end
if index_calculation==1   
    Merit_Best=999999999999999999999;
    Ratio_Lower2Upper_Best=1;
    Ratio_Upper2Reference_Best=1;
    for s=1:length(N_Shift)
        Cexp=Cunwrap+2*pi*N_Shift(s);
    for p=1:length(Thickness)
        Thickness_Now(1:length(Frequency_Considered),1:Active_Array_Size)=Thickness(p);%*Thickness_Profile/max(max(Thickness_Profile));

        for q=1:length(Ratio_Lower2Upper)
            for w=1:length(Ratio_Upper2Reference)
                    
                    %Thickness_Now=Thickness(p);
                    Ratio_Lower2Upper_Now=Ratio_Lower2Upper(q);
                    Ratio_Upper2Reference_Now=Ratio_Upper2Reference(w);
                    if length(n_should)==1
                        n(:,1:Active_Array_Size)=n_should;
                    else
                        n(:,1:Active_Array_Size)=repmat(n_should,1,Active_Array_Size);
                    end
                    Thickness_Now(1:length(Frequency_Considered),1:Active_Array_Size)=Thickness(p);%*Thickness_Profile/max(max(Thickness_Profile));
                    Current_Loop=1;
                    qqqq=1;
                    Current_Progress=100*((s-1)+((p-1)+((q-1)+(w-1)/length(Ratio_Upper2Reference))/length(Ratio_Lower2Upper))/length(Thickness))/length(N_Shift);

                    MSE_A=99999999;
                    MSE_B=99999999;
                    MSE_C=99999999;
                    MSE_total=99999999;
                    while (Current_Loop<Max_Number_of_Loop) && (MSE_total>MSE_total_OK)
                        

                        %C_Temp(size(C_Temp,1)+1,:)=2*C_Temp(size(C_Temp,1),:)-C_Temp(size(C_Temp,1)-1,:);
                        
                        
                        dA_dn_real = (A(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_real = (B(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_real = (C(n+dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        %dC_dn_real(size(dC_dn_real,1)+1,:)=2*dC_dn_real(size(dC_dn_real,1),:)-dC_dn_real(size(dC_dn_real,1)-1,:);
                        
                        dA_dn_imag = (A(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dB_dn_imag = (B(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        dC_dn_imag = (C(n+1i*dn,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dn;
                        %dC_dn_imag(size(dC_dn_imag,1)+1,:)=2*dC_dn_imag(size(dC_dn_imag,1),:)-dC_dn_imag(size(dC_dn_imag,1)-1,:);
                        
                        dA_dthickness = (A(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dB_dthickness = (B(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        dC_dthickness = (C(n,Thickness_Now+dt,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now)-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now))/dt;
                        %dC_dthickness(size(dC_dthickness,1)+1,:)=2*dC_dthickness(size(dC_dthickness,1),:)-dC_dthickness(size(dC_dthickness,1)-1,:);
                        
                        Q1=Q1_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q2=Q2_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q3=Q3_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q4=Q4_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q5=Q5_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q6=Q6_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q7=Q7_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q8=Q8_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);
                        Q9=Q9_P_add_Lambda_PT_inv(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness,Lambda);

                        
                        
                        delta_A=Aexp-A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_B=Bexp-B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);
                        delta_C=Cexp-C(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now);                        
                        MSE_A=(sum((delta_A./Aexp).^2)/length(delta_A)).^0.5;
                        MSE_B=(sum((delta_B./Bexp).^2)/length(delta_B)).^0.5;
                        MSE_C=(sum((delta_C./Cexp).^2)/length(delta_C)).^0.5;
                        MSE_total=((sum((delta_A./Aexp).^2)/length(delta_A))+(sum((delta_B./Bexp).^2)/length(delta_B))+(sum((delta_C./Cexp).^2)/length(delta_C))).^0.5;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        dn_real = det_array(delta_A,Q2,Q3,delta_B,Q5,Q6,delta_C,Q8,Q9)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        dn_imag = det_array(Q1,delta_A,Q3,Q4,delta_B,Q6,Q7,delta_C,Q9)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        dthickness = det_array(Q1,Q2,delta_A,Q4,Q5,delta_B,Q7,Q8,delta_C)./det_array(Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9);
                        %dn_real = det_array(delta_A,dA_dn_imag,dA_dthickness,delta_B,dB_dn_imag,dB_dthickness,delta_C,dC_dn_imag,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        %dn_imag = det_array(dA_dn_real,delta_A,dA_dthickness,dB_dn_real,delta_B,dB_dthickness,dC_dn_real,delta_C,dC_dthickness)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        %dthickness = det_array(dA_dn_real,dA_dn_imag,delta_A,dB_dn_real,dB_dn_imag,delta_B,dC_dn_real,dC_dn_imag,delta_C)./det_array(dA_dn_real,dA_dn_imag,dA_dthickness,dB_dn_real,dB_dn_imag,dB_dthickness,dC_dn_real,dC_dn_imag,dC_dthickness);
                        
                        n=n+dn_real+1i*dn_imag;
                        
                        if length(n_should)==1
                            n(isnan(n))=n_should;
                            n(isinf(n))=n_should;
                        else
                            n(isnan(n))=n_should(isnan(n));
                            n(isinf(n))=n_should(isinf(n));
                        end
                        % (The thickness is different for different wavelength) Thickness_Now=Thickness_Now+0.5*dthickness;   %How if I make the thickness a array during the iteration?
                        
                        if length(Thickness)==1 %if length(Thickness)>1 > Uniqueness test
                            Thickness_Now=Thickness_Now+mean(dthickness);   %How if I make the thickness a array during the iteration?
                        end
                        %Thickness_Now=Thickness_Now+0.01*sum(Weighting_Function_for_Thickness.*dthickness);   %How if I make the thickness a array during the iteration?
                        %Thickness_Now=Thickness_Now+0.01*dthickness;
                        %Thickness_Now=Thickness_Now+0.01*dthickness(find(Wavelength_micron_Considered<Wavelength_Thickness_Detection,1,'first'));
                        Current_Loop=Current_Loop+1;
                        fprintf('%d/%d     %f%%\n',Current_Loop,Max_Number_of_Loop,Current_Progress);
                        fprintf('%d/%d\n',MSE_A,MSE_A_OK);
                        fprintf('%d/%d\n',MSE_B,MSE_B_OK);
                        fprintf('%d/%d\n',MSE_C,MSE_C_OK);
                        fprintf('%d/%d\n',MSE_total,MSE_total_OK);
                        
                        if (max(isnan(Thickness_Now)>0))
                            Thickness_Now(isnan(Thickness_Now))=Thickness(p);%*Thickness_Profile(isnan(Thickness_Now))/max(max(Thickness_Profile));
                        end
                        if (max(isinf(Thickness_Now)>0))
                            Thickness_Now(isinf(Thickness_Now))=Thickness(p);%*Thickness_Profile(isinf(Thickness_Now))/max(max(Thickness_Profile));
                        end

                        
                    end
                        %fprintf('%f%%\n',Current_Progress);
                        Merit=sum((T_Considered-T_abs(n,Thickness_Now)).^2,1);
                        If_Merit_min=Merit<Merit_Best;
                        if max(If_Merit_min)
                            Merit_Best=Merit.*If_Merit_min+Merit_Best.*(1-If_Merit_min);
                        %dthickness_Best=dthickness;
                            Ratio_Lower2Upper_Best=Ratio_Lower2Upper_Now.*If_Merit_min+Ratio_Lower2Upper_Best.*(1-If_Merit_min);
                            Ratio_Upper2Reference_Best=Ratio_Upper2Reference_Now.*If_Merit_min+Ratio_Upper2Reference_Best.*(1-If_Merit_min);
                            T_Best=T_abs(n,Thickness_Now);
                            %C_Best=C_Temp;
                            n_Best=n;
                            
                            Thickness_Best=mean(Thickness_Now,1);
                            %Thickness_Best=Thickness_Now(Center_Wavelength_micron_Index,:);
                        end
                        %Thickness_Best_Record=(Thickness_Now);
                            

            end
        end
        
        if length(Thickness)>1 %if length(Thickness)>1 > Uniqueness test
            MSE_A_Record(p)=(sum((delta_A./Aexp).^2)/length(delta_A)).^0.5;
            MSE_B_Record(p)=(sum((delta_B./Bexp).^2)/length(delta_B)).^0.5;
            MSE_C_Record(p)=(sum((delta_C./Cexp).^2)/length(delta_C)).^0.5;
        end
        
        
    end
    
        if length(Thickness)==1 %if length(Thickness)>1 > Uniqueness test
            MSE_A_Record(s)=(sum((delta_A./Aexp).^2)).^0.5;
            MSE_B_Record(s)=(sum((delta_B./Bexp).^2)).^0.5;
            MSE_C_Record(s)=(sum((delta_C./Cexp).^2)).^0.5;
            MSE_A_Record_1(s)=(sum((delta_A).^2)).^0.5;
            MSE_B_Record_1(s)=(sum((delta_B).^2)).^0.5;
            MSE_C_Record_1(s)=(sum((delta_C).^2)).^0.5;
        end
    end
    
end



subplot(3,1,1)
plot(Frequency_Considered,Aexp,Frequency_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best));



subplot(3,1,2)
plot(Frequency_Considered,Bexp,Frequency_Considered,B(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best));


subplot(3,1,3)
plot(Frequency_Considered,Cexp,Frequency_Considered,C(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best));
%%
delta_freq=abs(Frequency_Considered(2)-Frequency_Considered(1));
C1_d =diff(unwrap(angle(1./r1(n).*(Ratio_Lower2Upper_Now.*(t1(n).*t1_r(n).*r2(n))))))/delta_freq;
C2_d =(4*pi.*real(n(1:(length(n)-1))).*Thickness_Now(1:(length(Thickness_Now)-1))/c);
C3_d =(4*pi.*Frequency_Considered(1:(length(Frequency_Considered)-1)).*(diff(real(n))/delta_freq).*Thickness_Now(1:(length(Thickness_Now)-1))/c);


subplot(3,1,1)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C1_d);
        xlabel('Wavelength (micron)');
        ylabel('C''_1');


subplot(3,1,2)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C2_d);
        xlabel('Wavelength (micron)');
        ylabel('C''_2');

subplot(3,1,3)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C3_d);
        xlabel('Wavelength (micron)');
        ylabel('C''_3');
        
subplot(3,1,1)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C1_d./(C1_d+C2_d+C3_d)*100);
        xlabel('Wavelength (micron)');
        ylabel('C''_1 (%)');
        %ylim([-1 1]);

subplot(3,1,2)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C2_d./(C1_d+C2_d+C3_d)*100);
        xlabel('Wavelength (micron)');
        ylabel('C''_2 (%)');
        %ylim([0 100]);
        
subplot(3,1,3)
plot(Wavelength_micron_Considered(1:(length(Frequency_Considered)-1)),C3_d./(C1_d+C2_d+C3_d)*100);
        xlabel('Wavelength (micron)');
        ylabel('C''_3 (%)');
        %ylim([0 100]);
%%

plot(Wavelength_micron_Considered,real(n_Best),Wavelength_micron_Considered,imag(n_Best));

plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered);

%legend('Calculated Transmission','Measured Transmission','Location','Best')
%ylim([0.4 1]);
%set(gca,'YTick',0.6:0.1:1)

%subplot(6,1,4)
%plot(Wavelength_micron_Considered,Aexp,Wavelength_micron_Considered,A(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));
%xlabel('Wavelength (micron)');
%ylabel('A');

%subplot(6,1,5)
%plot(Wavelength_micron_Considered,Bexp,Wavelength_micron_Considered,B(n,Thickness_Now,Ratio_Lower2Upper_Now,Ratio_Upper2Reference_Now));
%xlabel('Wavelength (micron)');
%ylabel('B');

%subplot(6,1,6)
%plot(Wavelength_micron_Considered,Cexp,Wavelength_micron_Considered,C_Temp);
%xlabel('Wavelength (micron)');
%ylabel('C');

%subplot(6,1,4)
%plot(2:(length(dn_real_record)-1),dn_real_record(3:end));
%xlabel('Number');
%ylabel('dn_real_record');

%subplot(6,1,5)
%plot(2:(length(dn_imag_record)-1),dn_imag_record(3:end));
%xlabel('Number');
%ylabel('dn_imag_record');

%subplot(6,1,6)
%plot(2:(length(dthickness_record)-1),dthickness_record(3:end));
%xlabel('Number');
%ylabel('dthickness_record');

%plot((1:length(Thickness_Best))*5/1000,Thickness_Best*1E6,(1:length(Thickness_Best))*5/1000,Thickness*1E6*Thickness_Profile(1,:)/max(max(Thickness_Profile)));
%xlabel('Lateral Position (mm)');
%ylabel('Thickness (micron)');
%legend('Calculated Thickness','Predicted Thickness (Lieanr to OPD)')

%%

%plot(Wavelength_micron_Considered,mean(T_Best,2),Wavelength_micron_Considered,T_Considered(:,1));
%xlabel('Wavelength (micron)');
%ylabel('T');

%dlmwrite('Lateral Position.txt',((1:length(Thickness_Best))*5/1000)','delimiter','\t','newline','pc');
%dlmwrite('Thickness.txt',Thickness_Best'*1E6,'delimiter','\t','newline','pc');
%dlmwrite('OPD.txt',Thickness_Profile(1,:)','delimiter','\t','newline','pc');

%dlmwrite('T_Best.txt',mean(T_Best(:,1:100),2),'delimiter','\t','newline','pc');

%dlmwrite('T_Considered.txt',T_Considered(:,1),'delimiter','\t','newline','pc');

%dlmwrite('Wavelength_micron_Considered.txt',Wavelength_micron_Considered,'delimiter','\t','newline','pc');

%dlmwrite('n.txt',mean(real(n_Best(:,1:100)),2),'delimiter','\t','newline','pc');

%dlmwrite('k.txt',mean(imag(n_Best(:,1:100)),2),'delimiter','\t','newline','pc');

%plot(Wavelength_micron_Considered,A(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Aexp);
%plot(Wavelength_micron_Considered,B(n_Best,Thickness_Best,Ratio_Lower2Upper_Best,Ratio_Upper2Reference_Best),Wavelength_micron_Considered,Bexp);
%plot(Wavelength_micron_Consi
fprintf('MSE_A=%e, MSE_B=%e, MSE_C=%f\n',MSE_A_Record,MSE_B_Record,MSE_C_Record);

MSE_Normalized=((MSE_A_Record.^2+MSE_B_Record.^2+MSE_C_Record.^2)/min(MSE_A_Record.^2+MSE_B_Record.^2+MSE_C_Record.^2)).^0.5;

%MSE_Normalized_by_sum=((MSE_A_Record.^2+MSE_B_Record.^2+MSE_C_Record.^2)/sum(MSE_A_Record.^2+MSE_B_Record.^2+MSE_C_Record.^2)).^0.5;
MSE_Unnormalized=(((MSE_A_Record.^2+MSE_B_Record.^2+MSE_C_Record.^2)).^0.5)/length(Frequency_Considered);


if length(Thickness)==1 %if length(Thickness)>1 > Uniqueness test
    if length(N_Shift)==1
                %% Plot
        subplot(3,1,1)
        plot(Wavelength_micron_Considered,real(n_Best));%,Wavelength_micron_Considered,n1_Considered(:,1));
        xlabel('Wavelength (micron)');
        ylabel('n');
        %ylim([1.55 1.651]);
        xlim([0.51 0.6]);


        %legend('Calculated Film Refractive Index','Refractive Index of Substrate','Location','Best')

        subplot(3,1,2)
        plot(Wavelength_micron_Considered,imag(n_Best));
        xlabel('Wavelength (micron)');
        ylabel('k');
        %ylim([0.008 0.012]);
        xlim([0.51 0.6]);

        %%

        subplot(3,1,3)
        plot(Wavelength_micron_Considered,T_Best,Wavelength_micron_Considered,T_Considered(:,1));
        xlabel('Wavelength (micron)');
        ylabel('T');
        %ylim([0.5 1]);
        xlim([0.51 0.6]);
        dlmwrite('Wavelength_nm_Considered.txt',Wavelength_micron_Considered*1000,'delimiter','\t','newline','pc');
        dlmwrite('n.txt',real(n_Best),'delimiter','\t','newline','pc');
        dlmwrite('k.txt',imag(n_Best),'delimiter','\t','newline','pc');
        dlmwrite('T.txt',T_Best,'delimiter','\t','newline','pc');
        dlmwrite('T_exp.txt',T_Considered,'delimiter','\t','newline','pc');
        Thickness_Best
    else
        mean(Thickness_Now)
        subplot(4,1,1)
        plot(N_Shift,MSE_A_Record);
        xlabel('N');
        ylabel('MSE_A');
    %xlim([4.65 4.8]);

        subplot(4,1,2)
        plot(N_Shift,MSE_B_Record);
        xlabel('N');
        ylabel('MSE_B');
        %xlim([4.65 4.8]);

        subplot(4,1,3)
        plot(N_Shift,MSE_C_Record);
        xlabel('N');
        ylabel('MSE_C');


        Uniqueness_Range=abs(N_Shift(find(MSE_Normalized<1.1,1,'first'))-N_Shift(find(MSE_Normalized<1.1,1,'last')));

        subplot(4,1,4)
        plot(N_Shift,MSE_Normalized);
        xlabel('N');
        ylabel('MSE_T_o_t_a_l');
        ylim([0 100]);
        [invalue minindex]=min(MSE_Normalized);
        N_Shift_min_MSE=N_Shift(minindex);
        %ylim([1 2.5]);
        
        dlmwrite('N_Shift.txt',N_Shift,'delimiter','\t','newline','pc');
        dlmwrite('MSE_Normalized.txt',MSE_Normalized','delimiter','\t','newline','pc');
        dlmwrite('MSE_Unnormalized.txt',MSE_Unnormalized','delimiter','\t','newline','pc');
    end
else
    subplot(4,1,1)
    plot(Thickness*1E6,MSE_A_Record);
    xlabel('Thickness (micron)');
    ylabel('MSE_A');

    subplot(4,1,2)
    plot(Thickness*1E6,MSE_B_Record);
    xlabel('Thickness (micron)');
    ylabel('MSE_B');

    subplot(4,1,3)
    plot(Thickness*1E6,MSE_C_Record);
    xlabel('Thickness (micron)');
    ylabel('MSE_C');
    %xlim([4.6 4.8]);
    Uniqueness_Range=abs(Thickness(find(MSE_Normalized<1.1,1,'first'))-Thickness(find(MSE_Normalized<1.1,1,'last')))*1E6;

    [invalue minindex]=min(MSE_Normalized);
    Thickness_min_MSE=Thickness(minindex)
    
    subplot(4,1,4)
    plot(Thickness*1E6,MSE_Normalized);
    xlabel('Thickness (micron)');
    ylabel('MSE_T_o_t_a_l');
    %xlim([4.6 4.8]);

    dlmwrite('Thickness_array.txt',Thickness*1E6,'delimiter','\t','newline','pc');
    dlmwrite('MSE_Normalized.txt',MSE_Normalized','delimiter','\t','newline','pc');
end

%xlim([4.64 4.8]);


cd('D:\Users\TuanShu\');

mean(abs((t1(n_Best).*t1_r(n_Best).*r2(n_Best)./r1(n_Best).*exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))))

dlmwrite('Position_micron.txt',Position_micron-Position_micron(maxindex_Signal_Sample_1),'delimiter','\t','newline','pc');
dlmwrite('Signal_Sample_1.txt',real(Signal_Sample_1)/max_Signal_Sample_2,'delimiter','\t','newline','pc');
dlmwrite('Signal_Sample_2.txt',real(Signal_Sample_2)/max_Signal_Sample_2,'delimiter','\t','newline','pc');

dlmwrite('Wavelength_nm_Considered.txt',Wavelength_micron_Considered*1000,'delimiter','\t','newline','pc');
dlmwrite('Aexp_OE version.txt',Aexp/max_Aexp,'delimiter','\t','newline','pc');
dlmwrite('Bexp_OE version.txt',Bexp/max_Aexp,'delimiter','\t','newline','pc');
dlmwrite('Cexp_OE version.txt',Cexp,'delimiter','\t','newline','pc');

%%
SD_Thickness=0.05; %micron
T_sum=0;
NNN=1000;
for j=1:NNN
    Thickness_Random=normrnd(Thickness_Best,SD_Thickness.*1E-6);
    T_sum=T_sum+T_abs(n_Best,Thickness_Random);
end
T_mean=T_sum/NNN;
plot(Wavelength_micron_Considered,T_mean,Wavelength_micron_Considered,T_Considered);
dlmwrite('T_ave.txt',T_mean,'delimiter','\t','newline','pc');

%film abs: abs(exp(1i*4*pi.*Frequency_Considered.*n.*Thickness_Now/c))
%目前看起來, 算出來的k都太大了, 穿透光只會剩下4成以下
%因此, 對到的ratio才會這麼小
%因此, 當給正確ratio (~8)的data進去時, 才會變成得要用>1的Ratio_Lower2Upper (也就是誤以為是因為rear
%interface干涉效率太高導致 ratio變大
%事實上, 我猜正確的ratio應該是在6左右, 但目前對出來的數字都是在3左右
%可能是T的公式搞錯了
%Ratio_Lower2Upper=1;     %if the rear interface has larger interference eff than front interface, Ratio_Lower2Upper>1
