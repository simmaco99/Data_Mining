function Analisi_Immagine()
close
% Programma che chiede all'utente un immagine e compie varie operazioni
% 1) mostra i valori singolari e comprimi l'immagine usando solo alcuni
%   valori singolari scelti dall'utente
% 2) comprime l'immagine mantenendo l'errore sotto una soglia e la
% compressione sopra un'altra
% 3) aggiunge il rumore all'immagine
% 4) applica la tecnica del filtraggio lineare basata sulla discretizzazione dell'equazione del calore, Sobel o Prewitt
% 5) segmenta l'immagine

dir = '/Users/simmaco/MATLAB/Data_Mining/Programmi/Immagini Test';
A=leggi_immagine(dir);
select= menu('Seleziona il tipo di operazione', 'Compressione','Aggiunta di rumore',"Filtraggio","Segmentazione");
switch select
    case 1
    A = double(A);
    select2 = menu("Scegli il tipo di compressione", "Scegliendo il numero dei valori singolari da usare","Scegliendo il numero massimo di valori singolari e il passo");
    switch  select2
        case 1
            S= svd(A(:,:,1));
            plot(S,'*') ;
            drawnow
            title("Valori singolari (1 canale)"); xlabel("numero valori singolari"); ylabel("valori singolari")
            k= input("Scegliere il numero di valori singolari da utilizzare (anche vettori)");
            [~,C,err]= Comprimi(A,k);
            Compressione = C';
            Errore= err;
            k = k';
            
        case 3
            C = input("Inserisci la percentuale di compressione: ");
            toll = input("Inserisci un'errore massimo: ");
            [~,k,C,err] =Comprimi2(A,toll,C);
            Compressione =C';
            Errore = err;
            k = k';
        case 2
            S= svd(A(:,:,1));
            plot(S,'*') ;
            drawnow
            title("Valori singolari (1 canale)"); xlabel("numero valori singolari"); ylabel("valori singolari")
            k= input("Scegliere il numero di valori singolari da utilizzare (anche vettori): ");
            passo =input("Scegli ogni quanti valori singolari disegnare l'immagine: ");
            k = passo:passo:k;
            k = [1 k];
            [~,C,err]= Comprimi(A,k);
            Compressione = C';
            Errore= err;
            k = k';
    end
   
    table(k, Compressione,Errore)
    figure();
    plot(k,Errore(:,1),'ro-')
    hold on
    plot(k,Errore(:,2),'go-')
    plot(k,Errore(:,3),'bo-')
    xlabel('numero valori singolari')
    ylabel('Errore in norma di frobenius')
    legend('rosso','verde','blu')
    case 2
        B=Rumore(A); %non funziona add e moltiplicativo
        drawnow
        sc=input("Vuoi filtrare l'immagine? (s/n) ","s");
        if sc =='s'
            M = input("Ogni quante iterazioni vuoi vedere l'immagine ");
            filtraggio_lineare(B,M); 
        end
    case 3
        M = input("Ogni quante iterazioni vuoi vedere l'immagine ");
        filtraggio_lineare(A,M);
    case 4
        %F = [ -1 0 1; -2 0 2; -1 0 1];
        Segmentazione(A);
end
end

function A =leggi_immagine(dir)
% A =leggi_immagine() permette di leggere una qualsiasi immagine presente
% nella directory dir (nel nome dei file delle immagini deve essere
% presente la stringa BN (se l'immagine e' in bianco e nero)
%
% INPUT
% dir   percorso dalla radice della directory delle immagini

loc = pwd;
cd(dir)
a = string(split(ls()));
select = menu("Seleziona l'immagine",a(1:end-1));
name = a(select);
A = imread(name);
if contains(name,'BN')
    A = rgb2gray(A);
end

cd(loc)
end


function  [Ak,C,err]=Comprimi(A,k)
close all;
% [Ak,C,err]=Comprimi(A,k) comprime l'immagine A usando la SVD
% troncata a un fissato numero di valori singolari
% INPUT
% A     matrice dell'immagine
% k     vettore contenente il numero dei valori singolari da usare
% OUTPUT
% Ak    ultima matrice compressa (usando k(end) valori singolari)
% C     fattori di compressione a seconda del numero di valori singolari
% err   norma dell'erore in norma di Frobenius a seconda del numero di valori singolari

s=size(A);
s= s(1:2);
b=size(A,3);
A =double(A);

U=zeros(s(1),s(1),b);
V = zeros(s(2),s(2),b);
S=zeros(size(A));
sigma = zeros(min(s),b);
for i =1:b
    [U(:,:,i),S(:,:,i),V(:,:,i)]=svd(A(:,:,i));
    sigma(:,i) =diag(S(:,:,i));
end

l = length(k);
C = (sum(s) +1)/prod(s)*k;
err = zeros(l,b);
Ak = zeros(s(1),s(2),b);
nc =2;
nr = floor(l/nc)+1;
subplot(nr,nc,1)
imshow(uint8(A))
title('Originale')

for i = 1:l
    for can = 1:b
        Ak(:,:,can)=U(:,1:k(i),can) * S(1:k(i),1:k(i),can)* V(:,1:k(i),can)';
    end
    subplot(nr,nc,i+1)
    str = sprintf("k = %d",k(i));
    imshow(uint8(Ak))
    title(str)
    for can =1:b
    err(i,can) = norm(A(:,:,can)-Ak(:,:,can),'fro')/norm(A(:,:,can),'fro');
    end
    
    
end
end

function [Ak,k,C,err]=Comprimi2(A,toll,C)
% [Ak,k,C,err]=Comprimi2(A,toll) comprime l'immagine A usando la SVD
% fino a quando non viene richiesta la compressione minima oppure l'errore
% non e' minore di quello richiesto
% INPUT
% A     matrice dell'immagine
% toll  errore massimo accettabile
% C     fattore di compressione desiderato
% OUTPUT
% Ak    matrice compressa
% k     numero di valori singolari utilizzati
% C     fattori di compressioni a seconda del numero di valori singolari
% err   norma di Frobenius dell'errore a seconda del numero di valori singolari

s =size(A(:,:,1));
b = size(A,3);
l = min(s);
err=zeros(l,b);
A =double(A);

Ak = zeros(s(1),s(2), b);
U=zeros(s(1),s(1),b);
V=zeros(s(2),s(2),b);
S=zeros(s(1),s(2),b);
sigma =zeros(min(s),b);

for i = 1: b
    [U(:,:,i),S(:,:,i),V(:,:,i)]=svd(A(:,:,i));
    sigma(:,i) =diag(S(:,:,i));
end

j=1;

while (j<l) && (j< C * prod(s)/(sum(s) + 1)) && (err(j) <=toll)
    j=j+1;
    for i = 1: b
        Ak(:,:,i) = Ak(:,:,i)+ sigma(j,i) .* U(:,j,i) *V(:,j,i)';
    end
    for can =1:b
        err(j-1,can) = norm(A(:,:,can)-Ak(:,:,can),'fro')./norm(A(:,:,can),'fro');
    end
end
k = 1:j;
C = (sum(s) +1)/prod(s)*k;
err = err(1:j-1,:);


end

function C=Rumore(A)
tipo = menu('Scegli il tipo di rumore', 'Poisson','Gaussiano','Sale e pepe');

ss= size(A);
bn =size(A,3);
s= ss(1:2);

switch tipo
    case 1
         C = imnoise(A,'poisson');   
    case 2
        C = imnoise(A,'gaussian');
    case 3
            A=double(A);
         C= A/255;
         perc = input('Percentuale di rumore: ');
         np = floor(prod(s)*perc); %numero di pixel da cambiare
         C = reshape(C,[prod(s),1,bn]);
         tipo= randi(2,np,1)-1;  
         dx = zeros(np,1,bn);
         for i =1:bn
             dx(:,:,i) = tipo;
         end
         C(randi(prod(s),np,1),1:bn) = dx;
         C =reshape(C,ss);
         C = C*255;
         C = uint8(C);
end

subplot(1,2,1) 
imshow(uint8(A));
title('Immagine originale');
subplot(1,2,2)
imshow(C);
title('Immagine con rumore');
end

function C=filtraggio_lineare(A,M,F)
%  C=filtraggio_lineare(A,M,F) applica il filtro lineare F ad un'immagine
%  e mostra a video il confronto tra la matrice originale e quella filtrata ogni M iterazioni
%  e chiede se continuare.
%  C=filtraggio_lineare(A,M) applica uno tra i seguenti filtri lineari
%       Discretizzazione dello schema del calore a 5 punti
%       Filtro di Sobel 
%       Filtro di Prewitt
% ad un immagine scelta dall'utente e mostra a video il confronto tra la
% matrice originale e quella filtrata ogni M iterazioni e chiede se
% continuare. 
% 
A=double(A);
subplot(1,2,1)
imshow(uint8(A));
title('Immagine originale');
k=0;
s = size(A);
%SI SPECIFICA IL FILTRO
if nargin ==2 %il filtro non viene passato
    
    filtro = menu('Scegliere il filtro','Equazione del calore','Sobel','Prewitt');

    switch filtro
        case 1; F= [ 0 .25 0 ; .25 0 .25; 0 .25 0 ];
        case 2; F = [ -1 0 1; -2 0 2; -1 0 1];
        case 3; F = [ -1 0 1; -1 0 1; -1 0 1];
    end
    if filtro~=1
        dir=menu('Specificare la direzione di applicazione','orizzontale','verticale');
        if dir==2
                F=F';
        end
    end
end
%PARTE DI FILTRAGGIO
while true
    C= zeros(size(A));
    k = k+1;
    for i =2:s(1)-1
        for j=2:s(2)-1
               
            C(i,j,:) =sum(sum(F.*A(i-1:i+1,j-1:j+1,:)));
        end
    end
    A=C;
    if mod(k,M)==0 
            subplot(1,2,2);imshow(uint8(C));str=sprintf('%d applicazioni del filtro',k); title(str); 
            drawnow
            cnt= menu('Continuare', 'si','no');
            if cnt ==2; break; end
    end
   
end
end

function C=Segmentazione(A)
F = [ -1 0 1; -1 0 1; -1 0 1];
A =double(A);
C= zeros(size(A));
s = size(A);
for i =2:s(1)-1
    for j=2:s(2)-1
        C(i,j,:) =(sum(sum(F.*A(i-1:i+1,j-1:j+1,:))).^2 +sum(sum(F'.*A(i-1:i+1,j-1:j+1,:))).^2).^(1/2);
    end
end
subplot(2,1,1); imshow(uint8(A)); title("Immagine originale")
subplot(2,1,2);imshow(uint8(C)); title("Immagine segmentata"); 


end
