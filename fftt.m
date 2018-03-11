function [f,espectro]=fftt(x,t)

[n,c]=size(x);
if n<c
	x=x.'; %x es ahora una columna
end
[n2,c]=size(t);
if n2<c
	t=t.';     %t es ahora una columna
end

[n]=size(x);
[n2]=size(t);

if n~=n2
    sprintf('Vectores no tienen la misma dimensión')
    return;
end

N=length(x);

if mod(N,2)~=0
    x=x(1:N-1);
    t=t(1:N-1);
    N=length(x);
end


dt=t(2:2:N)-t(1:2:N-1);
if norm(dt-dt(1))>1e-5
     sprintf('Tiempo no esta equiespaciado')
    return;
end   

T=dt(1);
F=1/T;

f = F*(0:N/2)/N;

espectro=fft(x,N); %FFT bruta
espectro=espectro/(N/2); %correccion de amplitud
espectro(1)=espectro(1)/2;%correccion de amplitud componente estatica
espectro=espectro(1:N/2+1);


[n,c]=size(espectro);
if n<c
	espectro=espectro.'; %espectro es ahora una columna
end
[n2,c]=size(f);
if n2<c
	f=f.';     %f es ahora una columna
end