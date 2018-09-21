clear all
close all
clc

f=fopen('axi.e','r');
x=fread(f,'single');

nsam = 2048*110;

xr=reshape(x(1:nsam),2048,10,11);

%for i=1:1:100
%  snap(1:10,1:8) = xr(i,1:10,1:8);
%  snap = snap'./max(max(snap));
%  imagesc(snap),caxis([0,0.7])
%  pause(0.4)
%end

for i=1:8
plot(xr(:,1,i)),hold on
xlim([1,100])
pause(0.4)
end
