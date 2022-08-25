clear all
close
col = 'krgbcmy';

dim = 20;  % dimensione data set
k_max = 20; 
passo =5;
select = menu("Scegli cosa classificare","Numeri","Lettere Maiuscole","Lettere minuscole", "Lettere","Numeri e lettere");

data = menu("Scegliere il file di dati","English Handwriting","English HandWriting preprocessati","Computer");

switch data
    case 1 
        file   =  pwd()+"Caratteri/Img/Sample%03d/img%03d-%03d.png";
        tot_imm=52;
        num_char=62;
    case 2
        file = pwd()+"Caratteri/Cut/Sample%03d/img%03d-%03d.png";
        tot_imm=52;
        num_char=10;
    case 3
        file = pwd()+"Caratteri/Fnt/Sample%03d/img%03d-%05d.png";
        tot_imm=1016;
        num_char = 62;
    case 4
end


if select~=1 && data==4
    error("Tale dataset contiene solo numeri")
end

num=string(0:9)';
L = string(char(double('A')+ (0:25)'));
l = string(char(double('a')+ (0:25)'));

switch select
    case 1 
       Caratteri = num;
        in = 1;
        fine = 10;
    case 2 
        Caratteri = L;
        in = 11;
        fine = 36;
    case 3
       Caratteri = l;
       in =37;
       fine =62;
    case 4
       Caratteri = [L; l];
       in = 11;
       fine = 62;
    case 5
        Caratteri = [num ;L ; l];
        in =1;
        fine = 62;
end
if data ==2  || data ==3
    A =double(imread(sprintf(file,1,1,1)));
else
    A= double(rgb2gray(imread(sprintf(file,1,1,1))));
end
[dim1,dim2]=size(A);
np =numel(A);

data_name = "Dati"+num2str(data)+".mat";
if ~contains(ls,data_name) && data~=4
tr = randperm(tot_imm,dim)';
train = zeros(np,dim,num_char);
U = zeros(np,k_max,fine-in+1);
for j = 1: fine-in+1
%     disp(j)
    c = in+j-1;
    for i = 1: dim
        if data==2 || data ==3
            A= double(imread(sprintf(file,c,c,tr(i))));
        else
            A= double(rgb2gray(imread(sprintf(file,c,c,tr(i)))));
        end
        train(:,i,j) = reshape(A,[np 1]);
    end
    [U(:,:,j),~,~]=svds(train(:,:,j),k_max);
end
end
caricamento = 0;

while true
   if bool~=1
   b=menu("Vuoi un altro carattere","Si","No");
   if b==2
       break
   end
   end
     clc
      close all
    bool =2;
   n = tr(1);
   while ~isempty(find(tr==n, 1))
       c = randi([in fine]); %scegli il carattere
       n = randi(tot_imm);
   end
   if data ==2 || data ==3
        A= double(imread(sprintf(file,c,c,n)));
   else
       A= double(rgb2gray(imread(sprintf(file,c,c,n))));
   end
   z = reshape(A,[np,1]);
   residuo = calcolo_err(z,U);
   
   subplot(2,1,1)
   imshow(uint8(reshape(z,dim1,dim2)))
   subplot(2,1,2)
   beauty = 1:passo:k_max;
   dd = length(beauty);
   for i = 1:dd 
       plot(residuo(:,beauty(i)),col(i),"LineWidth",2);
       hold on
   end
   set(gca,"XTick",1:length(Caratteri))
   set(gca,"XLim",[1 length(Caratteri)])
   set(gca,"XTickLabel",Caratteri)
   legend(string((beauty)'),'Location','best');

  
   T = table(Caratteri);
   T2 = array2table(residuo(:,beauty));
   T2.Properties.VariableNames="Res"+string(beauty);
   disp([T T2])
   
   [residuo,j] = sort(residuo(:,k_max));
   fprintf("I possibili caratteri sono\n")
   residuo = residuo(1:2);
   j = j(1:2);
   carattere = Caratteri(j);
   disp(table(carattere,residuo))
   fprintf("Il carattere era %c",Caratteri(c-in+1));

end



function res=calcolo_err(z,U)
    % dato un vettore calcola il vettore dei residui 
    sz = size(U);
    nchar = sz(3);
    nk = sz(2);
    res = zeros(nchar,nk);
   
    for i = 1: nchar
         dx =  U(:,:,i)'*z;
        for s = 1: nk
        res(i,s) = norm(z - U(:,1:s,i)*(dx(1:s)));
        end
    end
end



        

