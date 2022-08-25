function [A,x]=information_retrieval(file,keywords,query)
% INPUT
% file      string array con i nomi dei file da usare come documenti
% keywords  string array con le varie keywords
% query     string array con le parole delle query
% OUTPUT
% A matrice di 0/1 dove la i-esima colonna corrisponde al i-esimo file
% x vettore di preferenza 
%
% L'algoritmo funziona in due modi in base alla lunghezza di file:
% se length(file)==1 allora i ''documenti" tra i quali ricerchiamo sono
% tante quante le linee del documento denominato file
% se length(file)>1  allora i file contiene il nome di tutti i documenti
%
% ESEMPIO
% keywords = string({'marzo','aprile','giugno','presidente','consiglio','europa','argento','italia','concerto')


keywords= upper(keywords);
query = upper(query);
A = MakeMat(file,keywords);
b = zeros(length(keywords),1);
for i =1:length(query)
    b = b+ (keywords== query(i))';
end
%soluzione ai minimi quadrati
x=A\b;
disp([x;norm(b-A*x)/norm(b)])

[U,S,V] = svd(A);
plot(A'*U(:,1),A'*U(:,2),'o')
sigma = diag(S);
for k = 1: length(sigma)
    x(1,:) = U*
end

function A = MakeMat(file,b)
NMAX = 10000; %numero massimo di documenti apribili
if length(file)==1
    A= zeros(length(b),NMAX);
    % file e b cell array con nome dei file e nome dei keywords
    fid = fopen(file,'r');
    linea=upper(fgetl(fid));
    col=1;
    while ischar(linea)
        for i =1:length(b)
            key=b(i);
            if contains(linea,key)
                A(i,col) = 1;
            end
        end
        linea =upper(fgetl(fid));
        col=col+1;
    end
    fclose(fid);
    B= zeros(length(b),col);
    B(1:size(A,1),1:size(A,2))= A;
    A=B;
else
    A= zeros(length(b), length(file));
    for col = 1: length(file)
        fid = fopen(file(col),'r');
        linea = upper(fgetl(fid));
        i = 1: length(b); % contiene indice parole  ancora da cercare
        while ~isempty(i) && ischar(linea)
            for j = 1:length(i)
                key = b(i(j));
                if contains(linea,key)
                    A(i(j),col) = 1;
                    i(j)=[];
                end
            end
            linea = upper(fgetl(fid));
        end
        fclose(fid);
    end
end
end
