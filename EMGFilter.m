function [EMGfilt] = EMGFilter(EMGdata,hpvar,lpvar,bworder,samfreq)
%The function inputs a column of vector of raw EMG data and filter settings
%and passes it through a high pass filter, a rectifier, and then a low pass
%filter

%High pass filter
if hpvar ~= 0
  [bb,aa] = butter(bworder, hpvar/(samfreq/2),'high');
  datafilt=filter(bb,aa,EMGdata);
end
%Rectifying EMG data
datafilt=abs(datafilt);
%Low pass filter
if lpvar ~= 0
  [bb,aa] = butter(bworder, lpvar/(samfreq/2),'low');
  EMGfilt=filter(bb,aa,datafilt);
end
end