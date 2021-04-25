function charter
  
  
 N=1000;
 
 dx=0.01/N;
 
 dd=zeros(N,1);
 
 xx=zeros(N,1);
 
 for k=1:N
   x=1.61+k*dx;
   xx(k)=x;
   %  dd(k)=2*sqrt(1-1./x)+(1-1./x)-x;

    dd(k)=x^4-2*x^3-x^2+2*x+1;
    if(abs(dd(k))<1e-6)
    x
    dd(k)
    end
 %  dd(k)=sqrt(x-1./x)+sqrt(1-1./x)-x;
end
 y=1.6185*sqrt(2);
 y
 
  figure(5)
 plot(xx,dd,'-r','LineWidth',2);
 

%filename= 'D:\Works and Studies\Photonic\Lumerical\models\faradary_isolator\results-7-180s-108.txt';
[filename,filepath]=uigetfile('*.txt', 'Select results file')
file=strcat(filepath,strcat('\',filename));
data0=csvread(file,3,0);

[filename1,filepath1]=uigetfile('*.txt', 'Select results file')
file1=strcat(filepath1,strcat('\',filename1));
data1=csvread(file1,3,0);
     

data=data0;

ndata=size(data(:,1),1);

reverse=0;
if(reverse==1)
for j=1:ndata
data(j,1)=1./data0(ndata-j+1,1);
data(j,2)=data0(ndata-j+1,2);
data(j,3)=data0(ndata-j+1,3);
end

fid = fopen('results_lam.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'results vs lambda\n');  
  
fprintf(fid,'******\n');
fprintf(fid,'[Fn *  Rotation * Transmitance ]\n');  
  for p=1:ndata
  

  fprintf(fid,'%f, %f, %f\n',data(p,1),data(p,2),data(p,3));


  end
    fclose(fid);
end

Fn=data(:,1);


wn1=Fn(1);
wn2=Fn(ndata);

colT='-r';
colR='-r';


 figure(1)

 plot(Fn(:,1),data(:,3),colR,'LineWidth',2);
  hold on
 plot(data1(:,1),data1(:,3),'-.k','LineWidth',2);
 axis([wn1,wn2,0,1]);
 
 figure(2)
 plot(Fn(:,1),data(:,2),colT,'LineWidth',2);
   hold on
 plot(data1(:,1),data1(:,2),'-.k','LineWidth',2);

 tmax=1.1*max(abs(data(:,2)));
   axis([wn1,wn2,0,90]);

  hold on
end