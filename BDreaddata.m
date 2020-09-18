% Opens datafile and reads time series data
% Header of datafile should end with a line of ----
% Heather Graven 2020
function [tp,pp] = BDreaddata(filename)

fid1=fopen(filename,'r');
Header=textscan(fid1,'%s',1,'delimiter','\n'); % Read strings delimited by a carriage return
h=1;
while h==1
    Header=textscan(fid1,'%s',1,'delimiter','\n'); % Read strings delimited by a carriage return until '-' is encountered
    HeaderString=char(Header{1,1}(1));
    if numel(HeaderString)>0
        if strcmp(HeaderString(1),'-')==1 
        h=2;
        end
    end
end
InputText=textscan(fid1,'%n %n','delimiter','\t'); % Read data
fclose(fid1);    
    
tp=InputText{1,1};
pp=InputText{1,2};