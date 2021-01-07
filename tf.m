function tf

clear all

colors = {'-or', '-ob', '-oc', '-ok','-*r', '-*b', '-*c', '-*k','-xr', '-xb', '-xc', '-xk'};
%C = {'k','b','r','g','y'}; % Cell array of colros.



[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)
%fid = fopen('input.txt','rt');

line=getNewDataLine(fid);
numbs = str2num(line);
geometry=numbs(1);
rec=0;
if(geometry==1 && length(numbs)>1)
rec=numbs(2);
end
rec
line=getNewDataLine(fid);
numbs = str2num(line);
d1=numbs(1);
d2=0;
if(length(numbs)>1)
d2=numbs(2);
end
d3=0;
if(length(numbs)>2)
d3=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
fi=numbs(1);

line=getNewDataLine(fid);
numbs = str2num(line);
theta=numbs(1);
global inc_mode;
inc_mode=0;
if(length(numbs)>1)
inc_mode=numbs(2);
end

line=getNewDataLine(fid);
numbs = str2num(line);
ax=numbs(1);
ay=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
eps1=numbs(1);
eps3=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
ndata=length(numbs)

epsbx=numbs(1);
if(ndata>1)
epsby=numbs(2);
epsbz=numbs(3);
gamab=numbs(4);
end

line=getNewDataLine(fid);
numbs = str2num(line);
epsax=numbs(1);
if(ndata>1)
epsay=numbs(2);
epsaz=numbs(3);
gamaa=numbs(4);
end

if(ndata>1)

epsb=[epsbx 0 -1i*gamab;0 epsby 0;1i*gamab 0 epsbz ];

epsa=[epsax 0 -1i*gamaa;0 epsay 0;1i*gamaa 0 epsaz ];
else
epsb=[epsbx];

epsa=[epsax];
end

line=getNewDataLine(fid);
numbs = str2num(line);
Na=numbs(1);
global ndef;
global defstart;
ndef=0;
defstart=0;
if(length(numbs)>1)
ndef=numbs(2);
defstart=floor(Na/2);
end
if(length(numbs)>2)
defstart=numbs(3);
end
line=getNewDataLine(fid);
numbs = str2num(line);
transmit=numbs(1);
rotation=numbs(2);


if(transmit>0)
colT=colors{transmit};
else
colT=colors{1};
end

if(rotation>0)
colR=colors{rotation+1};
else
colR=colors{2};
end



line=getNewDataLine(fid);
numbs = str2num(line);

wn1=numbs(1);
wn2=numbs(2);

global ndiv;
ndiv=numbs(3);

line=getNewDataLine(fid);
numbs = str2num(line);
nGx=numbs(1)
nGy=numbs(2)

line=getNewDataLine(fid);
numbs = str2num(line);
plotFT=numbs(1);
plotWave=numbs(2);



%=================

t1=cputime;

dwn=(wn2-wn1)/ndiv;

cf=(plotWave==0);

Fr_hom=zeros(1*cf*ndiv+1,1);

 d=d1+d2+d3;
 
 TE=0;
for m=1:1

if(m==2)
 TE=1;
 colT='-ok';
end
for p=1:1*cf*ndiv+1
    
    %p
    
     Fn(p)=wn1+dwn*(p-1);
    
    wvlen=Fn(p);
   %  wvlen= d*Fn(p)
   %  wvlen=d/Fn(p);
    % wvlen
   % Fn(p)=wn1+dwn*(p-1);
   % k1=sqrt(eps1)*2*pi*Fn(p)/ay;
    % Fr_hom(p)=-180*gamab/sqrt(epsbx)*Na*Fn(p);
    

   [Ts Rs,Fr]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,d1,d2,d3,Na,wvlen,p,theta,fi,TE);
           
      if(real(Ts)>1) 
      %  Ts=1;
      end
      
    Tr(p)=real(Fr);
    Tt(p)=real(Ts);
    
    uu=Ts+Rs;

end


t2=cputime;

comptation_time=t2-t1;
comptation_time;
        
        Tr';
        Tt';
	
result=zeros(1*cf*ndiv+1,3);


fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'%d\t%d\n',nGx,nGy);
fprintf(fid,'[wn *  Rotation * Transmitance ]\n');  


  for p=1:1*cf*ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= Tr(p);
  result(p,3)= Tt(p);

  fprintf(fid,'%f\t%f\t%f\n',result(p,1),result(p,2),result(p,3));
  fprintf(fidx,'%f\t%f\n',result(p,1),Fr_hom(p));


  end
    fclose(fid);
     fclose(fidx);
     
  disp('Results:');
  disp('[wn *  Rotation * Transmitance ]');	
  disp(result);
  
if(rotation &&length(Tr)>1)
                figure(1)
             plot(Fn,Tr,colR);
             
 
             axis([wn1,wn2,-90,90]);
             hold on
             
        %   plot(Fn,Fr_hom,'+k');
        %     hold on

end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,colT);
                 axis([wn1,wn2,min(Tt),max(Tt)]);
                %axis([wn1,wn2,0,1]);
                 hold on
            
end

end

end



function [Ts Rs Fr]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,...
d1,d2,d3,Na,wvlen,p,theta,fi, TE)

nk=size(epsb,2);

global MM;

global Kapa;
global bE;
global ndiv;



 lam=wvlen;%*1e-6;
%%% eps
ev=1.602e-19;

wp=8.24;
wc=0.048;
c=3.0e8;
w=2*pi*c/lam;

em=1-wp^2/(w*(w+1i*wc));
em;

nAg=sqrt(em);

 Ts=cos(2*pi*p/ndiv)^2;
 Rs=0;
  Fr=0;%90*sin(2*pi*p/ndiv);
 
 e1=epsa(1,1);
 e2=epsb(1,1);
  n0=sqrt(1./eps1);
  
 A = 4.92719645;
 B = 7.27691471;
 C =11.5786091E+2;
 D = 42.7173925;
 E =100 ;
 
 lam2=lam*lam;
 
 nsq=A+B*lam2/(lam2-C)+D*lam2/(lam2-E);
 
 n1=sqrt(nsq);

 n2=1.46;
 
 n1=sqrt(e1);
 n2=sqrt(e2);

 %n1=2.17;
 %n2=1.49;

 %nAg=.01+3i;
 nAg=.5;%.9+.02i;
 n3=nAg;
 p0=cosd(theta)/n0;
 method=3;
 
if(method==0)

 beta1=2*pi*n1*d1*cosd(theta)/lam;
 if(TE==0)
  p1=cosd(theta)/n1;
else
  p1=cosd(theta)*n1;
  end
 M1=zeros(2,2);
 M2=zeros(2,2);
  M3=zeros(2,2);

M1(1,1)=cos(beta1);
M1(1,2)=-i*sin(beta1)/p1;
M1(2,1)=-i*sin(beta1)*p1;
M1(2,2)=M1(1,1);

 beta2=2*pi*n2*d2*cosd(theta)/lam;
 
 if(TE==0)
  p2=cosd(theta)/n2;
else
  p2=cosd(theta)*n2;
end

M2(1,1)=cos(beta2);
M2(1,2)=-i*sin(beta2)/p2;
M2(2,1)=-i*sin(beta2)*p2;
M2(2,2)=M2(1,1);

  
 beta3=2*pi*n3*d3*cosd(theta)/lam;
 if(TE==0)
  p3=cosd(theta)/n3;
else
  p3=cosd(theta)*n3;
  end
M3(1,1)=cos(beta3);
M3(1,2)=-i*sin(beta3)/p3;
M3(2,1)=-i*sin(beta3)*p3;
M3(2,2)=M3(1,1);

  ML=M1*M2;


   MML=ML;

  for n=2:Na
    MML=MML*ML;

  end
  
   nf=n1;

  
  MM=MML;

  m11=MM(1,1);
  m12=MM(1,2);
  m21=MM(2,1);
  m22=MM(2,2);
  
  t1=2*p0/(m11+m12*p0+m21+m22*p0);
  
  t2=t1*conj(t1);
  
elseif(method==1)

 M1=zeros(2,2);
 M2=zeros(2,2);
  M3=zeros(2,2);

 %%%%%%%%%%% layer 1
%  Xi1=1-sind(theta)^2/n1^2;
  Xi1=sqrt(1-sind(theta)^2/n1^2);
  phi1=2*pi*n1*d1/lam*Xi1;
  
   P1=zeros(2,2);
   P1(1,1)=exp(1i*phi1);
   P1(2,2)=exp(-1i*phi1);
  
  D1=zeros(2,2);
   D1(1,1)=1;
     D1(1,2)=-n1*Xi1;
   D1(2,1)=1;
    D1(2,2)=n1*Xi1;
    
    M1=inv(D1)*P1*D1;


%%%%%%%%%%% layer 2
  Xi2=sqrt(1-sind(theta)^2/n2^2);
  phi2=2*pi*n2*d2/lam*Xi2;
  
   P2=zeros(2,2);
   P2(1,1)=exp(1i*phi2);
   P2(2,2)=exp(-1i*phi2);
  
  D2=zeros(2,2);
   D2(1,1)=1;
     D2(1,2)=-n2*Xi2;
   D2(2,1)=1;
    D2(2,2)=n2*Xi2;
    
    M2=inv(D2)*P2*D2;

  ML=M1*M2;
  MU=M2*M1;

   MML=ML;
   MMU=MU;

  for n=2:Na/2
    MML=MML*ML;
    MMU=MMU*MU;
  end
  
  %MM=MML*MMU;
 nf=n1;
 df=d3;
   betaf=2*pi*nf*df*cosd(theta)/lam;
 if(TE==0)
  pf=cosd(theta)/nf;
else
  pf=cosd(theta)*nf;
end

  Mf=zeros(2,2);
Mf(1,1)=cos(betaf);
Mf(1,2)=-i*sin(betaf)/pf;
Mf(2,1)=-i*sin(betaf)*pf;
Mf(2,2)=Mf(1,1);

MM=MML*Mf*MMU;
 % MM=MMU;
  

   %%%%%%%%%%%
    
    D0=zeros(2,2);
    
       D0(1,1)=1;
     D0(1,2)=1;
   D0(2,1)=cosd(theta);
    D0(2,2)=-cosd(theta);
  
  ss=1;
   if(ss==1)
   MM=inv(D0)*MM*D0;
   end
  
    m11=MM(1,1);
  m12=MM(1,2);
  m21=MM(2,1);
  m22=MM(2,2);

  
  if(ss==0)
  t1=2*p0/(m11+m12*p0+m21+m22*p0);
  
  t2=t1*conj(t1);
  else
t2=1./(m11*conj(m11));
end  

elseif(method==2)

 M1=zeros(2,2);
 M2=zeros(2,2);
 %%%%%%%%%%% layer 1
kz1=2*pi*n1*cosd(theta)/lam;
  phi1=kz1*d1;
  
   P1=zeros(2,2);
   P1(1,1)=exp(1i*phi1);
   P1(2,2)=exp(-1i*phi1);
  
  D1=zeros(2,2);
   D1(1,1)=1;
     D1(1,2)=1;
   D1(2,1)=-kz1;
    D1(2,2)=kz1;
    
    M1=D1*P1*inv(D1);


%%%%%%%%%%% layer 2
  kz2=2*pi*n2*cosd(theta)/lam;
  phi2=kz2*d2;
   P2=zeros(2,2);
   P2(1,1)=exp(1i*phi2);
   P2(2,2)=exp(-1i*phi2);
  
  D2=zeros(2,2);
   D2(1,1)=1;
     D2(1,2)=1;
   D2(2,1)=-kz2;
    D2(2,2)=kz2;
    
  M2=D2*P2*inv(D2);

  ML=M1*M2;
%MU=M2*M1;
  MU=M1*M2;
   MML=ML;
   MMU=MU;

  for n=2:Na/2
    MML=MML*ML;
    MMU=MMU*MU;
  end
 % n3=n1;
   %%%%%%%%%%% defect layer
   kz3=2*pi*n3*cosd(theta)/lam;
  phi3=kz3*d3;

   P3=zeros(2,2);
   P3(1,1)=exp(1i*phi3);
   P3(2,2)=exp(-1i*phi3);
  
  D3=zeros(2,2);
   D3(1,1)=1;
     D3(1,2)=1;
   D3(2,1)=-kz3;
    D3(2,2)=kz3;
    
    M3=D3*P3*inv(D3);
  
 % MM=MML*M3*MMU;
   MM=MML*MMU;

   %%%%%%%%%%%
  kz0=2*pi*n0*cosd(theta)/lam;  
    D0=zeros(2,2);
    
   D0(1,1)=1;
   D0(1,2)=1;
   D0(2,1)=-kz0;
   D0(2,2)=kz0;

%   M=D0*MM*inv(D0);
    M=inv(D0)*MM*D0;

    m11=M(1,1);
  
  t2=1./(m11*conj(m11));

elseif(method==3)

 M1=zeros(4,4);
 M2=zeros(4,4);
 %%%%%%%%%%% layer 1

 dn1=0.1*n1;
 n1p=n1+dn1;
  n1m=n1-dn1;

 beta1p=2*pi*n1p*d1*cosd(theta)/lam;
  p1p=cosd(theta)/n1p;
  
   beta1m=2*pi*n1m*d1*cosd(theta)/lam;
  p1m=cosd(theta)/n1m;
  
M1(1,1)=cos(beta1p);
M1(1,2)=-i*sin(beta1p)/p1p;
M1(2,1)=-i*sin(beta1p)*p1p;
M1(2,2)=M1(1,1);
M1(3,3)=cos(beta1m);
M1(3,4)=-i*sin(beta1m)/p1m;
M1(4,3)=-i*sin(beta1m)*p1m;
M1(4,4)=M1(3,3);


%%%%%%%%%%% layer 2

 dn2=0.0*n2;
 n2p=n2+dn2;
  n2m=n2-dn2;
  
 beta2p=2*pi*n2p*d2*cosd(theta)/lam;
  p2p=cosd(theta)/n2p;
  
   beta2m=2*pi*n2m*d2*cosd(theta)/lam;
  p2m=cosd(theta)/n2m;
  
M2(1,1)=cos(beta2p);
M2(1,2)=-i*sin(beta2p)/p2p;
M2(2,1)=-i*sin(beta2p)*p2p;
M2(2,2)=M2(1,1);
M2(3,3)=cos(beta2m);
M2(3,4)=-i*sin(beta2m)/p2m;
M2(4,3)=-i*sin(beta2m)*p2m;
M2(4,4)=M2(3,3);


  D0=zeros(4,4);  
  D0(1,1)=1;
  D0(1,2)=1;
  D0(2,1)=cosd(theta);
  D0(2,2)=-cosd(theta);
  D0(3,3)=1;
  D0(3,4)=1;
  D0(4,3)=cosd(theta);
  D0(4,4)=-cosd(theta);  
    
  ML=M1*M2;
  MU=M2*M1;
 % MU=M1*M2;
   MML=ML;
   MMU=MU;

  for n=2:Na/2
    MML=MML*ML;
    MMU=MMU*MU;
  end
  
   %%%%%%%%%%% layer in
 M3=zeros(4,4);
   n3=n1;
  d3=0*d1; 
 dn3=0.0*n3;
 n3p=n3+dn3;
  n3m=n3-dn3;

 beta3p=2*pi*n3p*d3*cosd(theta)/lam;
  p3p=cosd(theta)/n3p;
  
   beta3m=2*pi*n3m*d3*cosd(theta)/lam;
  p3m=cosd(theta)/n3m;
  
M3(1,1)=cos(beta3p);
M3(1,2)=-i*sin(beta3p)/p3p;
M3(2,1)=-i*sin(beta3p)*p3p;
M3(2,2)=M3(1,1);
M3(3,3)=cos(beta3m);
M3(3,4)=-i*sin(beta3m)/p3m;
M3(4,3)=-i*sin(beta3m)*p3m;
M3(4,4)=M3(3,3);

  
  MM=MML*M3*MMU;
  
    D0i=inv(D0);

   MM=D0i*MM*D0;

  
   m11=MM(1,1);
  m12=MM(1,2);
  m21=MM(2,1);
  m22=MM(2,2);
  m33=MM(3,3);
  t2=1./(m11*conj(m11));
t2=1./(m33*conj(m33));

Fr=0.5*atan(m11/m33)/pi*180;

end
  
  Ts=abs(t2);

end


function line=getNewDataLine(fid)

TF=1;
k=0;
while (k<100 && TF==1)
line=fgets(fid);
TF = strncmpi(line,'/',1);
k=k+1;

end


end


function [h]=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The function plotting 3-dimensional arrow
%
% h=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The inputs are:
%       x,y,z =  vectors of the starting point and the ending point of the
%           arrow, e.g.:  x=[x_start, x_end]; y=[y_start, y_end];z=[z_start,z_end];
%       head_frac = fraction of the arrow length where the head should  start
%       radii = radius of the arrow
%       radii2 = radius of the arrow head (defult = radii*2)
%       colr =   color of the arrow, can be string of the color name, or RGB vector  (default='blue')
%
% The output is the handle of the surfaceplot graphics object.
% The settings of the plot can changed using: set(h, 'PropertyName', PropertyValue)
%
% example #1:
%        arrow3d([0 0],[0 0],[0 6],.5,3,4,[1 0 .5]);
% example #2:
%        arrow3d([2 0],[5 0],[0 -6],.2,3,5,'r');
% example #3:
%        h = arrow3d([1 0],[0 1],[-2 3],.8,3);
%        set(h,'facecolor',[1 0 0])
%
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% July 2010 (C)

if nargin==5
    radii2=radii*2;
    colr='blue';
elseif nargin==6
    colr='blue';
end
if size(x,1)==2
    x=x';
    y=y';
    z=z';
end

x(3)=x(2);
x(2)=x(1)+head_frac*(x(3)-x(1));
y(3)=y(2);
y(2)=y(1)+head_frac*(y(3)-y(1));
z(3)=z(2);
z(2)=z(1)+head_frac*(z(3)-z(1));
r=[x(1:2)',y(1:2)',z(1:2)'];

N=10;
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

X1=[];Y1=[];Z1=[];
j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii*cos(theta).*(P1-Pc) + radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(2:3,j)=R1(:,1);
    Y1(2:3,j)=R1(:,2);
    Z1(2:3,j)=R1(:,3);
    j=j+1;
end

r=[x(2:3)',y(2:3)',z(2:3)'];

dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;

normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;

j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii2*cos(theta).*(P1-Pc) + radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(4:5,j)=R1(:,1);
    Y1(4:5,j)=R1(:,2);
    Z1(4:5,j)=R1(:,3);
    j=j+1;
end

X1(1,:)=X1(1,:)*0 + x(1);
Y1(1,:)=Y1(1,:)*0 + y(1);
Z1(1,:)=Z1(1,:)*0 + z(1);
X1(5,:)=X1(5,:)*0 + x(3);
Y1(5,:)=Y1(5,:)*0 + y(3);
Z1(5,:)=Z1(5,:)*0 + z(3);

h=surf(X1,Y1,Z1,'facecolor',colr,'edgecolor','none');

light('Color','r')


end

