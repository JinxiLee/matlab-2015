clear all

cd('D:\AMO\150805');

Data_without_filter=dlmread('without filter.txt');        
Data_with_filter=dlmread('with filter.txt');        
Data_without_filter=Data_without_filter(1:1000,2)'-mean(Data_without_filter(1:200,2));
Data_with_filter=Data_with_filter(1:1000,2)'-mean(Data_with_filter(1:200,2));
[max1 index1]=max(Data_without_filter);
[max2 index2]=max(Data_with_filter);

Data_without_filter=circshift(Data_without_filter,[0 500-index1])/max1;
Data_with_filter=circshift(Data_with_filter,[0 500-index2])/max2;

plot(1:1000,Data_without_filter,1:1000,Data_with_filter);
legend('Without FGL495M','With FGL495M');
ylim([-1.1 1.1]);
%% Raw Data Reading

sampling_rate=128000/4000;  %Hz

scanning_speed=1;   %micron/sec

Pre_ave=1;

for p=1:length(Data)/Pre_ave
    Data_Pre(p)=sum(Data((Pre_ave*(p-1)+1):(Pre_ave*p)));
end
Data_Pre=Data_Pre';
dPosition=scanning_speed/sampling_rate*Pre_ave; %micron


X=(dPosition:dPosition:dPosition*length(Data_Pre))';


plot(X,Data_Pre);
xlabel('Position (micron)');
ylabel('Intensity (a.u.)');


%% Background Substration
Background_start=Data_Pre(1);
Background_end=Data_Pre(end);

Background=interp1([X(1) X(end)]',[Background_start Background_end]',X);

plot(X,Data_Pre,X,Background);

Data_sub=Data_Pre-Background;

plot(X,Data_sub);

%% Filtering
LB=310;
HB=500;

Data_FFT=fft(Data_sub);

plot(real(Data_FFT));

Data_FFT_Filtered=Data_FFT;

Data_FFT_Filtered(1:LB)=0;
Data_FFT_Filtered(HB:end)=0;

Data_New=ifft(Data_FFT_Filtered);

plot(X,abs(Data_New));

%% Averaging

Averging_Factor=1;

Data_New_Averaged=smooth(abs(Data_New),Averging_Factor);

% SNR calc
Noise_Floor=abs(Data_New_Averaged((length(Data_New_Averaged)-999):length(Data_New_Averaged)));

plot(Noise_Floor);
Error=std(Noise_Floor);

SNR=log10(max(abs(Data_New_Averaged))/Error)*20;

