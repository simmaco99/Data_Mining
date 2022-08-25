close all
clearvars -except U1 U2;
load fisheriris
bool=1;
while bool
k = input("Inserisci il numero di cluster: ");
it = input("Inserisci il numero di volte che vuoi iterare l'algoritmo di Llyoid: ");
select = menu("Scegliere la distanza da utilizzare","euclidean","seuclidean","cityblock","chebychev","squaredeuclidean","minkowski");

switch select
    case 1
        dist = "euclidean";
    case 2
        dist = "seuclidean"; 
    case 3 
        dist ="cityblock";
    case 4 
        dist = "chebychev";
    case 5 
        dist= "squaredeuclidean";
    case 6
        dist= "minkowski";
end

d_fun =@(x,y) pdist2(x,y,dist);

x = meas(:,1:2);

C = KMEANSS(x,d_fun,k,it);
plot_cluster(x,C);
drawnow
s = input("Vuoi eseguire una nuova clusterizzazione? (s/n)","s");
if s=='n'
    bool = 0;
end
end

clearvars -except U1 U2;

function C = KMEANSS(x,d_fun,k,it)
[C,cost]=llyoid(x,d_fun,k);
for i = 2:it
    [Cn,costn]=llyoid(x,d_fun,k);
    if costn<cost
        C=Cn;
        cost=costn;
    end
end
end

function [C,cost]=llyoid(x,d_fun,k)
% Inizializzazione: si scelgono casualmente i k-centroidi 
N = size(x,1);
s = randperm(N);
s=s(1:k);
z = x(s,:); 
D = d_fun(x,z);
[~,C]=min(D,[],2);

Cn=0;
while ~isequal(C,Cn)
    z = calcola_centroidi(x,C);
    D = d_fun(x,z);
    Cn = C;
    [~,C]=min(D,[],2);
    
end
cost = sum(min(D,[],2));
end

function z = calcola_centroidi(x,C)
    k = max(C);
    z=zeros(k,size(x,2));
    for i = 1: k
        z(i,:) =mean(x(C==i,:));
    end
end

function plot_cluster(x,C)
col = 'krgbcmy';
k = max(C);
for i = 1: k 
    plot(x(C==i,1),x(C==i,2), [col(i) 'o'],'MarkerFaceColor' ,col(i), 'MarkerSize',6);
    hold on
end
legend("cluster"+string((1:k)'))
end

