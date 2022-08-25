clc;
close;
clearvars -except U1 U2;
select = menu("Seleziona il tipo di dati","Grasso corporeo","Prezzo delle case");

switch select
	case 1
	    A = table2array(readtable("/Users/simmaco/MATLAB/Data_Mining/Programmi/Bodyfat.csv"));
        col=3;
    case 2
       	A = table2array(readtable("/Users/simmaco/MATLAB/Data_Mining/Programmi/houses.csv"));
        col = 4;
end
	
corr_sign =0.9; %valori significativi con correlazione in modulo maggiore di corr_sign
perc = 0.80; %percentuale di significativia' spiegata dalla componente principale
Z=zscore(A);
bool=1;
while bool
[n,p]=size(Z);
R=corrcoef(Z);

C_pos=(R -eye(size(R)))>=corr_sign;
[ip,jp]=find(C_pos); %indici riga e colonna maggiori di corr_sign
ip=ip(1:length(ip)/2);
jp = jp(1:length(jp)/2);

C_neg=-(R-eye(size(R)))>=corr_sign;
[in,jn]=find(C_neg); %indici riga e colonna maggiori di corr_sign
in=in(1:length(in)/2);
jn= jn(1:length(jn)/2);

if ~isempty(ip)
    fprintf("Sono correlati positivamente\n")
    fprintf("%d\t%d\n",[ip jp]')
end
if  ~isempty(in)
    fprintf("Sono correlati negativamente\n")
    fprintf("%d\t%d\n",[in jn]')
end

[V,lambda]=eig(R,'vector');

% oridino gli autovalori e cambio ordine degli autovalori
V = V(:,end:-1:1);
lambda= lambda(end:-1:1);

somma_lambda = cumsum(lambda);
perc_lamda = somma_lambda/somma_lambda(end);
dim_PC = sum(perc_lamda<perc)+1; %numero di componenti principali
fprintf("Usando %d componenti principli, spieghiamo il %2.f %\n %%",dim_PC,perc_lamda(dim_PC)*100)
y = Z * V(:,1:dim_PC); %le colonne sono le componenti principlali

for i = 1 : p
    close all
    x = Z(:,i);
    y= y(:,1);
    centro =[mean(x) mean(y)];
    dist = centro - [ x y];
    dist = sqrt(sum(dist.^2,2));

    er = (dist>=mean(dist)+3*var(dist));
    er=find(er);
    numeri = 1:n;
    numeri(er)=[];
     plot(Z(numeri,i),y(numeri,1),'.')
     hold on
     plot(Z(er,i), y(er,1),'r*')
    for j=1:n
        text(Z(j,i), y(j,1),num2str(j))
    end
    C= corrcoef(Z(:,i),y(:,1));
    C =C(1,2);
    title([num2str(i) ' variabile: coefficiente di correlazione: ' num2str(C,2)])
    for k = 1: numel(er)
        drawnow
        str= "\n"+num2str(er(k)) + " eliminare? (y/n) " ;
        bol=input(str,'s');
        if bol=="y"
            Z(er(k),:)=[];
            y(er(k),:)=[];
            n=n-1;
        end
    end
end
righe = p/col;
for i = 1:p
    subplot(righe,col,i)
    x = Z(:,i);
    y= y(:,1);
    plot(Z(:,i),y(:,1),'.')
    C= corrcoef(Z(:,i),y(:,1));
    C =C(1,2);
    title([num2str(i) ' var: coeff. correlazione: ' num2str(C,2)])
end
drawnow
in = input("\nVuoi ripetere l'analisi? (s/n)","s");
if in =='n'
    bool=0;
end
end

clearvars -except U1 U2;
