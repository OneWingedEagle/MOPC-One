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


%fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'%d\t%d\n',nGx,nGy);
fprintf(fid,'[wn *  Rotation * Transmitance ]\n');  


  for p=1:1*cf*ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= Tr(p);
  result(p,3)= Tt(p);

  fprintf(fid,'%f, %f, %f\n',result(p,1),result(p,2),result(p,3));
 % fprintf(fidx,'%f\t%f\n',result(p,1),Fr_hom(p));


  end
    fclose(fid);
    % fclose(fidx);
     
  disp('Results:');
  disp('[wn *  Rotation * Transmitance ]');	
  disp(result);
  
if(rotation &&length(Tr)>1)
                figure(1)
             plot(Fn,Tr,colR);
             
 
             axis([wn1,wn2,0,40]);
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



 lam=wvlen;
  d=d1;
%%% eps

%c=3.0e8;
omega=2*pi/lam;


Ts=1;
 Rs=0;
  Fr=0;%90*sin(2*pi*p/ndiv);
  
   e1=epsb(1,1);
   gam=epsb(1,3)/(-1i);

  D0=DynamicMatrix(1,0,omega);
    D0;
    omega
    MM1=TransferMatrix(e1,gam,omega,d);
    
   MM=MM1;
   
  for n=2:Na
    MM=MM1*MM;
  end

  TT=inv(D0)*MM*D0;
 %  TT=D0*MM*inv(D0);
 
  angle=TT(3,3)/TT(1,1);
  
  Fr=arg(angle)*90./pi();
  
  t2=1./(TT(1,1)*conj(TT(1,1)));

 Ts=abs(t2);
  
end


function [T]=TransferMatrix(eps,gamma,omega,d)
  
     Q=gamma/eps;

  Dn=DynamicMatrix(eps,Q,omega);
  
  Pn=PropagationMatrix(eps,Q,omega,d);

  D0=DynamicMatrix(1,0,omega);

   T=Dn*Pn*inv(Dn);
    
end

function [Dn]=DynamicMatrix(e1, Q,omega)

  pv1=(1./sqrt(2))*transpose([1 1i 0]);
  pv2=(1./sqrt(2))*transpose([-1 -1i 0]);
  pv3=(1./sqrt(2))*transpose([1 -1i 0]);  
  pv4=(1./sqrt(2))*transpose([-1 1i 0]);
  

  kz=omega*[0 0 1]';
  
  kv1=sqrt(e1*(1+Q))*kz;
  kv2=-sqrt(e1*(1+Q))*kz;
  kv3=sqrt(e1*(1-Q))*kz; 
  kv4=-sqrt(e1*(1-Q))*kz; 

  
  qv1=1/omega*cross(kv1,pv1);
  qv2=1/omega*cross(kv2,pv2);    
  qv3=1/omega*cross(kv3,pv3);
  qv4=1/omega*cross(kv4,pv4);
  
  Dn=[pv1(1) pv2(1) pv3(1) pv4(1);...
     qv1(2) qv2(2) qv3(2) qv4(2);...
    pv1(2) pv2(2) pv3(2) pv4(2);...
    qv1(1) qv2(1) qv3(1) qv4(1)];
end

function [Pn]=PropagationMatrix(e1, Q,omega,d)


  kz=omega*[0 0 1]';
  
  kv1=sqrt(e1*(1+Q))*kz;
  kv2=-sqrt(e1*(1+Q))*kz;
  kv3=sqrt(e1*(1-Q))*kz; 
  kv4=-sqrt(e1*(1-Q))*kz; 
  
   phi1=kv1(3)*d;
   phi2=kv2(3)*d;;
   phi3=kv3(3)*d;;  
   phi4=kv4(3)*d;;
   
   Pn=zeros(4,1);
   
   Pn(1,1)=exp(-1i*phi1);
    Pn(2,2)=exp(-1i*phi2);
     Pn(3,3)=exp(-1i*phi3);
      Pn(4,4)=exp(-1i*phi4);

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


