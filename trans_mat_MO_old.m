function trans_mat_MO_old

clear all

colors = {'-or', '-ob', '-oc', '-ok','-*r', '-*b', '-*c', '-*k','-xr', '-xb', '-xc', '-xk'};

[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)
%%  ==== Indicent angle
line=getNewDataLine(fid);
numbs = str2num(line);
theta= numbs(1);

%%=========  stack befor defect
dd1=[0 0 0];
eps1=[1 1 1];
gama1=[0 0 0];
line=getNewDataLine(fid);
numbs = str2num(line);
dd1(1)=numbs(1);

if(length(numbs)>1)
dd1(2)=numbs(2);
end
if(length(numbs)>2)
dd1(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
eps1(1)=numbs(1);

if(length(numbs)>1)
eps1(2)=numbs(2);
end
if(length(numbs)>2)
eps1(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
gama1(1)=numbs(1);

if(length(numbs)>1)
gama1(2)=numbs(2);
end
if(length(numbs)>2)
gama1(3)=numbs(3);
end
%%=============

%%=========   defect stack
ddef=[0 0 0];
epsdef=[1 1 1];
gamadef=[0 0 0];
line=getNewDataLine(fid);
numbs = str2num(line);
ddef(1)=numbs(1);

if(length(numbs)>1)
ddef(2)=numbs(2);
end
if(length(numbs)>2)
ddef(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
epsdef(1)=numbs(1);

if(length(numbs)>1)
epsdef(2)=numbs(2);
end
if(length(numbs)>2)
epsdef(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
gamadef(1)=numbs(1);

if(length(numbs)>1)
gamadef(2)=numbs(2);
end
if(length(numbs)>2)
gamadef(3)=numbs(3);
end
%%=============

%%=========  stack befor defect
dd2=[0 0 0];
eps2=[1 1 1];
gama2=[0 0 0];
line=getNewDataLine(fid);
numbs = str2num(line);
dd2(1)=numbs(1);

if(length(numbs)>1)
dd2(2)=numbs(2);
end
if(length(numbs)>2)
dd2(3)=numbs(3);
end

line=getNewDataLine(fid);
numbs = str2num(line);
eps2(1)=numbs(1);

if(length(numbs)>1)
eps2(2)=numbs(2);
end
if(length(numbs)>2)
eps2(3)=numbs(3);
end
line=getNewDataLine(fid);
numbs = str2num(line);
gama2(1)=numbs(1);

if(length(numbs)>1)
gama2(2)=numbs(2);
end
if(length(numbs)>2)
gama2(3)=numbs(3);
end
%%==========================


line=getNewDataLine(fid);
numbs = str2num(line);
N1=numbs(1);
Ndef=0;
N2=0;
if(length(numbs)>1)
Ndef=numbs(2);
end
if(length(numbs)>2)
N2=numbs(3);
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

%=================

t1=cputime;

dwn=(wn2-wn1)/ndiv;

Fr_hom=zeros(ndiv+1,1);

for p=1:ndiv+1
    
    %p
    
     Fn(p)=wn1+dwn*(p-1);
    
    wvlen=Fn(p);


   [Ts Rs,Fr]=TransferMatrixMultiLayer(eps1,epsdef,eps2,gama1,gamadef,...
   gama1,dd1,ddef,dd2,N1,Ndef,N2,wvlen,p,theta);
           
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
	
result=zeros(ndiv+1,3);


%fidx = fopen('analyt.txt','wt');  % Note the 'wt' for writing in text mode
fid = fopen('results.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'[nGx *  nGy]\n');  
  
fprintf(fid,'**************b\n');  
fprintf(fid,'[wn *  Rotation * Transmitance ]\n');  


  for p=1:ndiv+1
  
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



function [Ts Rs Fr]=TransferMatrixMultiLayer(eps1,epsdef,eps2,...,
  gama1,gamadef,gama2,dd1,ddef,dd2,N1,Ndef,N2,wvlen,p,theta)


 lam=wvlen;

omega=2*pi/lam;

Ts=1;
 Rs=0;
  Fr=0;
 
  %$$$$============= stacks before defect
  T1=TransferMatrix(eps1(1),gama1(1),omega,dd1(1));
  T2=TransferMatrix(eps1(2),gama1(2),omega,dd1(2));
  T3=TransferMatrix(eps1(3),gama1(3),omega,dd1(3));
  %%%%%%%%%%%%%%%
   T123=T1*T2*T3;
   TT1=T123;
     
  for n=2:N1
    TT1=T123*TT1;
  end
  
    %$$$$============= defect layer
  if(Ndef>0)
  T1=TransferMatrix(epsdef(1),gamadef(1),omega,ddef(1));
  T2=TransferMatrix(epsdef(2),gamadef(2),omega,ddef(2));
  T3=TransferMatrix(epsdef(3),gamadef(3),omega,ddef(3));
   T123=T1*T2*T3;
   Tdef=T123;
  for n=2:Ndef
    Tdef=T123*Tdef;
  end
  
end


  %$$$$============= stacks after defect
  if(N2>0)
  T1=TransferMatrix(eps2(1),gama2(1),omega,dd2(1));
  T2=TransferMatrix(eps2(2),gama2(2),omega,dd2(2));
  T3=TransferMatrix(eps2(3),gama2(3),omega,dd2(3));
  %%%%%%%%%%%%%%%
   T123=T1*T2*T3;
   TT2=T123;
  for n=2:N2
    TT2=T123*TT2;
  end
end

  T=TT1;
  if(Ndef>0)
  T=Tdef*T;
  end;
  if(N2>0)
  T=TT2*T;
  end;
  
  
  D0=DynamicMatrix(1,0,omega);
  
  TT=inv(D0)*T*D0;
 
  ang=TT(3,3)/TT(1,1);
  
  Fr=angle(ang)*90./pi();
  
  t2=1./(TT(1,1)*conj(TT(1,1)));

 Ts=abs(t2);

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


function [T]=TransferMatrix(eps,gamma,omega,d)
  
     Q=gamma/eps;

  Dn=DynamicMatrix(eps,Q,omega);
  
  Pn=PropagationMatrix(eps,Q,omega,d);

   T=Dn*Pn*inv(Dn);
    
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
