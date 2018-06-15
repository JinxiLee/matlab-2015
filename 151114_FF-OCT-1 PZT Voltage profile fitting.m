clear all
cd('D:\Users\TuanShu\141211_PZT waveform file generation');

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

%% Fitting

 g = fittype( @(a, b, c, x) a*x+b*(x.^0.5)+c);
 
 Fit_Result = fit(Time(1:floor(length(Time)/2)),Voltage(1:floor(length(Voltage)/2)),g);
 
 Voltage_fit=Fit_Result.a*Time+Fit_Result.b*(Time.^0.5)+Fit_Result.c;
 
 plot(Time,Voltage,Time,Voltage_fit);
Fit_Result.a
Fit_Result.b
Fit_Result.c

%%  filtering the signal