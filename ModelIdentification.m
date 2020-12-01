
% folder='C:\MyCloud\OneDriveUcf\Real\VariableSpring\Setup_Test';
lpvar=1;

[file,path]=uigetfile('*.csv');
data=importdata(append(path,file));
DStime=data.data(1,6);
samfreq=1/DStime;
 [bb,aa] = butter(4, lpvar/(samfreq/2),'low');
  FilteredData=filter(bb,aa,data.data(:,[2 3]));
plot(FilteredData)