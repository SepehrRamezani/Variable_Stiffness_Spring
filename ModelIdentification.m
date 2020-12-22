
% folder='C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test';
close all
lpvar=5;
ForceCo=[0.1491 0];
PressureCo=[6.6467 52.761];

[file,path]=uigetfile('C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test\*.csv');
data=importdata(append(path,file));
data.data(:,3)=data.data(:,3)-11.1;
DStime=data.data(1,6);
samfreq=1/DStime;
 [bb,aa] = butter(4, lpvar/(samfreq/2),'low');
  FilteredData=filtfilt(bb,aa,data.data(:,3));
indx=(FilteredData<-0.2 & data.data(:,3)>-6.2 );  
%   diff(FilteredData(:,3))
% plot(FilteredData(indx,1))
newData=data.data(indx,:);
lpvar=5;
[bb,aa] = butter(4, lpvar/(samfreq/2),'low');
FilteredData2=filtfilt(bb,aa,newData(:,[2,3,4]));
% hold on
% plot(diff(FilteredData(:,1))./DStime);
plot(newData(:,[2,3,4]));
hold on
plot(FilteredData2);