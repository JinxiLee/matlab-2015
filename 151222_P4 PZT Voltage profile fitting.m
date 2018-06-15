clear all
cd('D:\MATLAB\141211_PZT waveform file generation');

%%
Data=dlmread('Waveform_f0.004_V10_O1.txt');        
plot(Data(:,1),Data(:,2));

Original_Velocity=1;

Time=Data(:,1);
Time(500000)=[];
Voltage=Data(:,2);
Voltage(end)=[];

plot(Time,Voltage);

%% chnage to ms

Time=Time*1000;

%% Change speed

V_Increasing_Ratio=100/430;

Time_New=Time/V_Increasing_Ratio;

%% Known fitting result

A_known=1.646*10^(-5); 
B_known=2.4167*10^(-3); 
C_known=-1.7256*10^(-1); 

Voltage_Known=A_known*Time+B_known*(Time.^0.5)+C_known;

plot(Time_New,Voltage,Time,Voltage_Known);
%% Fitting

 g = fittype( @(a, b, c, x) a*x+b*(x.^0.5)+c);
 
 Fit_Result = fit(Time(1:floor(length(Time)/2)),Voltage(1:floor(length(Voltage)/2)),g);
 
 Voltage_fit=Fit_Result.a*Time+Fit_Result.b*(Time.^0.5)+Fit_Result.c;
 
 plot(Time,Voltage,Time,Voltage_fit);
Fit_Result.a
Fit_Result.b
Fit_Result.c


%%  filtering the signal