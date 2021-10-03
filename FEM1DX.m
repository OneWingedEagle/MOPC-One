

function FEM1D
  tic
clear


EPS0=8.854E-12;
XMU0=1.2566306E-6;
T=0;

global Nelm;
global Nnodes;

cryst_const=1;
x0=0;
x1=20;
Nelm=100;
dx=(x1-x0)/Nelm;

eps1=5.5+0.0*1j;
eps2=5.5+0*1j;



Nnodes=Nelm+1;


global indices;
indices=zeros(Nnodes+2,1);

ix=0;


for i=1:Nnodes
if(i==1 || i==Nnodes)
indices(i)=0;
 continue;
 end;
ix=ix+1;
indices(i)=ix;
end

ix=ix+1;
indices(Nnodes+1)=ix;
ix=ix+1;
indices(Nnodes+2)=ix;

global xcoord;
xcoord=linspace(x0,x1,Nnodes);
global elms;
elms=zeros(Nelm,2);
global epss;
epss=zeros(Nelm,1);

for i=1:Nelm
  elms(i,1)=i;
  elms(i,2)=i+1;
end

global neq;
global K;
global b;
%c=1/sqrt(EPS0*XMU0);
%lam=0.5e-6;

wn1=0.3;
wn2=.9;
wn1=1.801;
wn2=1.83;
divw=30;
T=zeros(divw+1,1);
wvec=zeros(divw+1,1);
dw=(wn2-wn1)/divw;
global k1;

for iw=1:divw+1;
wn=wn1+(iw-1)*dw;
wvec(iw)=wn;

k1=2*pi*wn/cryst_const;
f=k1;
lam=1./f;

w2c2=k1^2;


for ie=1:Nelm
  
  if(ie<Nelm/4 || ie>Nelm*3/4)
  epss(ie)=eps1;
  else
 epss(ie)=eps2;
 end

end


neq=0;
for i=1:Nnodes
if(indices(i)>0)
 neq=neq+1;
 end;
end

neq=neq+2;



SetK_Mat(w2c2);

b=zeros(neq,1);

SetRHS(k1,w2c2);

K;
sol=K\b;

T(iw)=abs(sol(neq));

end

  figure(2)
    
  plot(wvec,T,'-or');
 hold on

Ez=complex(Nnodes,1);
%Ezm=zeros(Nnodes,1);
Ezn=zeros(Nnodes,1);
for n=1:Nnodes
 
  index= indices(n);

     if(index>0) 
      Ez(n)=sol(index);
  end
end


for n=Nnodes+1:Nnodes+2
 
  index= indices(n);

  if(n==Nnodes+1)
  Ez(1)=sol(index);
  end
  if(n==Nnodes+2)
  Ez(Nnodes)=sol(index);
  end

end

for n=1:Nnodes
 
   x=xcoord(n);
   E0=1*exp(1i*k1*x);
   Ez(n)=Ez(n)+E0;
       
  Ezn(n)= abs(Ez(n));
end

 colorRe='r';
 colorIm='b';
 
    figure(1)
    
    plot(xcoord,real(Ez),colorRe);
     hold on
     plot(xcoord,imag(Ez),colorIm);
      hold on
     plot(xcoord,Ezn,'k');
    hold on

toc

end

function SetK_Mat(w2c2)
global K;
global neq;
global Nelm;
global Nnodes;
global elms;
global indices;
global epss;
K=zeros(neq,neq);
S=zeros(neq,neq);
for ie=1:Nelm
    gg=grdNi_gradNj(ie);
    
    NiNj=Ni_dot_Nj(ie);
   eps=epss(ie);
 % w2c2eps=epss(ie)*w2c2;
   for j=1:2
    n1=elms(ie,j); 
    indx1=indices(n1);
     if(indx1>0)
      for k=1:2
       n2=elms(ie,k);
       indx2=indices(n2);
       if(indx2>0)
         K(indx1,indx2)=  K(indx1,indx2)+gg(j,k);
         S(indx1,indx2)=  S(indx1,indx2)+eps*NiNj(j,k);
       end
      end
     end
    end
  end
   

for ie=1:Nelm
    
    NiNjMix=Ni_dot_NjMix(ie);
   eps=epss(ie);
 % w2c2eps=epss(ie)*w2c2;
   for j=1:2
    n1=Nnodes+j; 
    indx1=indices(n1);
     if(indx1>0)
      for k=1:2
       n2=elms(ie,k);
       indx2=indices(n2);
       if(indx2>0)
         S(indx1,indx2)=  S(indx1,indx2)+eps*NiNjMix(j,k);
       end
      end
     end
    end
  end

  for ie=1:1
    gg=grdNi_gradNjTR(ie);
  
    NiNj=Ni_dot_NjTR(ie);
   eps=epss(ie);
 % w2c2eps=epss(ie)*w2c2;
   for j=1:2
    n1=Nnodes+j; 
    indx1=indices(n1);
     if(indx1>0)
      for k=1:2
        n2=Nnodes+k; 
       indx2=indices(n2);
       if(indx2>0)
         K(indx1,indx2)=  K(indx1,indx2)+gg(j,k);
         S(indx1,indx2)=  S(indx1,indx2)+eps*NiNj(j,k);
       end
      end
     end
    end
  end

K=K+w2c2*S;
end


function SetRHS(k,w2c2)
global b;
global elms;
global xcoord;
global Nelm;
global indices;
global epss;

k2=k*k;

for ie=1:Nelm
  
    n1=elms(ie,1); 
    n2=elms(ie,2);
    
    eps=epss(ie);
    w2c2eps=eps*w2c2;

    x1=xcoord(n1);
    x2=xcoord(n2);
    
    d=x2-x1;
    
    br1=(cos(k*(d+x1))-cos(k*x1))/(d*k2)-sin(k*x1)/k;
    br2=-(cos(k*(d+x1))-cos(k*x1))/(d*k2)+sin(k*(x1+d))/k;
    
    bm1=(sin(k*(d+x1))-sin(k*x1))/(d*k2)+cos(k*x1)/k;
    bm2=-(sin(k*(d+x1))-sin(k*x1))/(d*k2)-cos(k*(x1+d))/k;
    
     index1=indices(n1);
     index2=indices(n2);
     
     if(index1>0)
     b(index1)= b(index1)-br1*(w2c2eps-k2)-1i*bm1*(w2c2eps-k2);
     end;
     if(index2>0)
     b(index2)= b(index2)-br2*(w2c2eps-k2)-1i*bm2*(w2c2eps-k2);
     end;     

end

  
end

function [gg]=grdNi_gradNj(ie)
gg=zeros(2,2);
global elms;
global xcoord;

  n1=elms(ie,1);
  n2=elms(ie,2);
  dx=xcoord(n2)-xcoord(n1);
  g=1./dx;
  
  gg(1,1)=g*g*dx;
  gg(1,2)=-g*g*dx;
  gg(2,1)=-g*g*dx;
  gg(2,2)=g*g*dx;
  
end

function [NiNj]=Ni_dot_Nj(ie)
NiNj=zeros(2,2);
global elms;
global xcoord;

  n1=elms(ie,1);
  n2=elms(ie,2);
  dx=xcoord(n2)-xcoord(n1);

  NiNj(1,1)=dx/3;
  NiNj(1,2)=dx/6;
  NiNj(2,1)=dx/6;
  NiNj(2,2)=dx/3;

  
end


function [gg]=grdNi_gradNjTR(ie)
gg=zeros(2,2);
global elms;
global xcoord;
global Nelm;
  n1=elms(1,1);
  n2=elms(Nelm,2);
  dx=xcoord(n2)-xcoord(n1);
  g=1./dx;
  
  gg(1,1)=g*g*dx;
  gg(1,2)=-g*g*dx;
  gg(2,1)=-g*g*dx;
  gg(2,2)=g*g*dx;
  
end

function [NiNj]=Ni_dot_NjTR(ie)
NiNj=zeros(2,2);
global elms;
global xcoord;
global Nelm;
  n1=elms(1,1);
  n2=elms(Nelm,2);
  dx=xcoord(n2)-xcoord(n1);

  NiNj(1,1)=dx/3;
  NiNj(1,2)=dx/6;
  NiNj(2,1)=dx/6;
  NiNj(2,2)=dx/3;

  
end



function [NiNj]=Ni_dot_NjMix(ie)
NiNj=zeros(2,2);
global elms;
global xcoord;
global Nelm;
global k1;

  n11=elms(1,1);
  n22=elms(Nelm,2);
  L=xcoord(n22)-xcoord(n11);
  
  n1=elms(ie,1);
  n2=elms(ie,2);
   
  x1=xcoord(n1);
  x2=xcoord(n2);
  d=x2-x1;
  div=10;
  dx=d/div;
  sum11=0;
  sum12=0;
    sum21=0;
  sum22=0;
  for k=1:div
    
    x=x1+(k-1)*dx;
   
  yy1=(1-x/L)*exp(-1i*k1*x);;
  yy2=x/L*exp(1i*k1*x);;
  
  y1=(-x/d+n1+1);
  y2=(x/d-n2+2);
  sum11=sum11+yy1*y1*dx;
  sum12=sum12+yy1*y2*dx;
  
  sum21=sum21+yy2*y1*dx;
  sum22=sum22+yy2*y2*dx;
  
  end

    NiNj(1,1)=sum11;
  NiNj(1,2)=sum12;
  NiNj(2,1)=sum21;
  NiNj(2,2)=sum22;

  
end

