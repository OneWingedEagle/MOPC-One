function matsolver
 
[filename,filepath]=uigetfile('*.txt', 'Select results file')
file=strcat(filepath,strcat('\',filename));
file
data=csvread(file,0,0);

nrow=int64(data(1,1));
ncol=int64(data(1,2));
ndata=size(data(:,1),1);
M=zeros(nrow,ncol);
ndata
for k=2:ndata-nrow
  r=int64(data(k,1))+1;
  c=int64(data(k,2))+1;
 % if(data(k,3)>-1e8)
 % data(k,3)
 % end
  M(r,c)=data(k,3)+1i*data(k,4);
end 

b=zeros(nrow,1);
ix=0;
for k=ndata-nrow+2:ndata
  ix=ix+1;
  b(ix)=data(k,1)+1i*data(k,2);
end


Mt=transpose(M);

C=Mt*M;
bt=Mt*b;


x=M\b;

xc=zeros(ncol,2);
for k=1:ncol
  xc(k,1)=real(x(k));
  xc(k,2)=imag(x(k));
end


fid = fopen('x.txt','wt');  % Note the 'wt' for writing in text mode
for k=1:ncol
  

 fprintf(fid,'%f %f\n',real(x(k)),imag(x(k)));


end
    fclose(fid);


end