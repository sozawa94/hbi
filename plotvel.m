% sample matlab code plotting a snapshot of slip velocity in 2D simulation with single processor
coord=importdata('xyz0_0.dat');
x=coord(:,1);
n=size(x);
n=n(1);
fileID = fopen('vel0_0.dat');
A = fread(fileID,'double');
nd=size(A);
nt=int16(nd(1)/n);
vel=reshape(A, [n nt]);
plot(x,vel(:,10)); 
