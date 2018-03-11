function [A,V,Pos,t,sensores,G,Mag,tmag]=cadera(str)
%Todos los datos deben estar en la misma carpeta. Ir a "HOME"->"Set Path"->Elegir carpeta con los datos
%Los archivos de los datos deben estar en .csv

%prompt = '¿Nombre del Archivo?  '; 
%str = input(prompt,'s'); % convierte el nombre del archivo ingresado en un string
M= csvread(str,1,0);%lee el archivo con los datos. NO considera títulos (primera fila)
[F,C]=size(M); %guarda dimensiones de M
ctf=2; %contador. considera que el primer tiempo es nulo
ctf2=2;
C1=1; %factor de conversión para pasar de mg a m/s^2
for i=1:F
    if M(i,3)~=0
        ctf=ctf+1; %cuenta el largo de las columnas de datos (tiempo distinto de cero)
    else
    end
end

for i=1:F
    if M(i,15)~=0
        ctf2=ctf2+1; %cuenta el largo de las columnas de datos (tiempo distinto de cero)
    else
    end
end
ctf=ctf-1;
ctf2=ctf2-1;
A=zeros(3,ctf,C/20); %Almacenan datos de aceleración para cada sensor. hay C/8 sensores. *la numeración no necesariamente es consecutiva
t=zeros(1,ctf);%vector de tiempo para graficar
tmag=zeros(1,ctf2);
t=M(1:ctf,3);
tmag=M(1:ctf2,15);
% Aceleraciones
for i=1:ctf
    y=1;
    for j=4:20:C
        A(1,i,y)=M(i,j)*C1; %convierte de mg a m/s^2
        A(2,i,y)=M(i,j+2)*C1;
        A(3,i,y)=M(i,j+4)*C1;
        y=y+1;
    end
end
% for i=1:3
%     for k=1:y-1
%     A(i,:,k)=A(i,:,k)-sum(A(i,:,k))/ctf;
%     end
% end
% Velocidades
for i=1:y-1
    V(1,:,i) = cumtrapz(t,A(1,:,i)); %integra aceleración
    V(2,:,i) = cumtrapz(t,A(2,:,i));
    V(3,:,i) = cumtrapz(t,A(3,:,i));    
end
% Posiciones
for i=1:y-1
    Pos(1,:,i) = cumtrapz(t,V(1,:,i)); %integra velocidad
    Pos(2,:,i) = cumtrapz(t,V(2,:,i));
    Pos(3,:,i) = cumtrapz(t,V(3,:,i));    
end
sensores=y-1;
%Giroscopio
G=zeros(3,ctf,C/20); %Almacenan datos de los giroscopios en grados/s para cada sensor. hay C/8 sensores. *la numeración no necesariamente es consecutiva
for i=1:ctf
    y1=1;
    for j=10:20:C
        G(1,i,y1)=M(i,j); 
        G(2,i,y1)=M(i,j+2);
        G(3,i,y1)=M(i,j+4);
        y1=y1+1;
    end
end

%Magnetómetro
Mag=zeros(3,ctf2,C/20); %Almacenan datos de los giroscopios en grados/s para cada sensor. hay C/8 sensores. *la numeración no necesariamente es consecutiva
for i=1:ctf2
    y1=1;
    for j=16:20:C
        Mag(1,i,y1)=M(i,j); 
        Mag(2,i,y1)=M(i,j+2);
        Mag(3,i,y1)=M(i,j+4);
        y1=y1+1;
    end
end
end