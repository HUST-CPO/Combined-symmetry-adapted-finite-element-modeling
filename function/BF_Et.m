function [out1,out2]=BF_Et(num,u,v)
if num==1
	out1=1-v;
	out2=u;
elseif num==2
	out1=v;
	out2=1-u;
elseif num==3
	out1=-v;
	out2=u;
end