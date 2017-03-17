//my_GA.m
function  [gy]= my_GA( popsize,pc,pm,max_iter )
%A simple GA to solve unconstrained optimization problem
% by Dong Longrui March,17 ,2017

%encoding routine
%length of a chromsome is 21+18=39.
l=39;
pop=round(rand(popsize,l)); %initial population is randomly generated
t=1;
%decoding routine
fitvalue = my_fitness( pop );
[ bestone,bestfit,s(t) ] = my_bestone( pop,fitvalue );
xb_1(t)=decodechrom( bestone(1,1:21) );
xb_2(t)=decodechrom( bestone(1,22:39) );
yb(t)=bestfit;

while(t<max_iter)
    %Crossover
    newpop1=my_crosssover( pop,pc );
    %Mutation
    newpop2=my_mutation( newpop1,pm );
    %Selection
    newpop3=my_selection(newpop2,fitvalue);
    pop=newpop3;
    t=t+1;
    fitvalue = my_fitness( pop );
    [ bestone,bestfit,s(t) ] = my_bestone( pop,fitvalue );
    xb_1(t)=decodechrom( bestone(1,1:21) );
    xb_2(t)=decodechrom( bestone(1,22:39));
    yb(t)=bestfit;
    [g_y(t),I]=max(yb);
    
end
gx_1=-3.0+xb_1(I)*15.1/(2^21-1);
gx_2=4.1+xb_2(I)*1.7/(2^18-1);
gy=g_y(t);

fprintf('Optimal value: gy = %0.5f\n',gy);
fprintf('Optimal solution: gx_1=%0.5f\n',gx_1);
fprintf('gx_2=%0.5f\n',gx_2);
fprintf('The generation of optimal solution=%d\n',I);

figure(1)
n=1:max_iter;
plot(n,yb,'g');
hold on;plot(n,g_y,'r');
% plot(I,0,'.');
plot(I,gy,'*');
plot([I I], [0 gy],'--k');text(I,gy,['(',num2str(I),',',num2str(gy),')']);
set(gcf,'color','w');
xlabel('Generation');ylabel('z');
legend('Green:Optimal solution each generation','Red:The best solution up till now');
title('Evolution Process');
% figure(2)
% xa1=-3.0:0.001:12.1;xa2=4.1:0.001:5.8;
% [x1,x2]=meshgrid(xa1,xa2);z=21.5+x1.*sin(x1.*4*pi)+x2.*sin(x2.*20*pi);
% mesh(x1,x2,z);
% hold on,plot3(gx_1,gx_2,gy,'p');

end

//decodechrom.m
function [ depop ] = decodechrom( pop )
%   二进制染色体解码
[~,py]=size(pop);
for i=1:py
    depop(:,i)=2.^(py-1).*pop(:,i);
    py=py-1;
end
depop=sum(depop,2);

end

//exam_problem.m
function [ objvalue ] = exam_problem( x_1,x_2 )
%   A example of problem

%objvalue=21.5+pop1.*sin.(pop1.*4*pi)+pop2.*sin.(pop2.*20*pi);
x_1=-3.0+x_1*15.1/(2^21-1);
x_2=4.1+x_2*1.7/(2^18-1);
objvalue=21.5+x_1*sin(x_1*4*pi)+x_2*sin(x_2*20*pi);

end

//my_fitness.m
function fitvalue = my_fitness( pop )
%   calculate fitnesses of the population 
[px,~]=size(pop);
for i=1:px
    fitvalue(i)=0;
end
pop1=pop(:,1:21);  %x1 chromsome
pop2=pop(:,22:39); %x2 chromsome
depop1=decodechrom(pop1);  
depop2=decodechrom(pop2);  %对染色体的两部分分别解码
for i=1:px
    fitvalue(i)=exam_problem( depop1(i),depop2(i) );
    %计算种群个体的适应度值
end 


end

//my_bestone.m
function [ bestone,bestfit,s ] = my_bestone( pop,fitvalue )
%   hold the best solution in the pocket by now
s=1;
[n,~]=size(pop);
bestone=pop(1,:);
bestfit=fitvalue(1);
for i=2:n
    if bestfit<=fitvalue(i)
        bestfit=fitvalue(i);
        bestone=pop(i,:);
        s=i;
    end
    
end

end

my_crossover.m
function [ newpop1 ] = my_crosssover( pop,pc )
%   One-Cut Point Crossover
[px,py]=size(pop);
for k=1:(px/2)
    if pc>=round(rand(1,1))
        i=0;
        j=0;
        while(i==j)
            i=1+(px-1)*round(rand(1,1));
            j=1+(px-1)*round(rand(1,1));
            
        end
        p=2+(py-2)*round(rand(1,1)); %p:the cut position
        
        pop(i,:)=cat(2,pop(i,1:p-1),pop(j,p:py));
        pop(j,:)=cat(2,pop(j,1:p-1),pop(i,p:py));
    end
    
    
    
    
end


newpop1=pop;








end

//my_mutation
function [ newpop2 ] = my_mutation( pop,pm )
%   mutate in 基因位
[px,py]=size(pop);
for k=1:px
    for j=1:py %j:cut position
        if pm>=rand(1,1)
            pop(k,j)=1-pop(k,j); % mutate at jth gene
        end
    end
    
end
newpop2=pop;

end

//my_selection
function [ newpop3 ] = my_selection( pop,fitvalue )
%   Roulette wheel selection
[px,~]=size(pop);
totalfit=sum(fitvalue);  %calculate the total fitness
for i=1:px
    p(i)=fitvalue(i)/totalfit;   %calculate the selection probability
end
for k=1:px
    q(k)=0;
    for j=1:k
        q(k)=q(k)+p(j);     %calculate the cumulative probability    
    end   
end

% for k=1:popsize
%     r=rand(1,1);
%     if r<=q(1)
%         newpop3(k)=pop(1);
%     end
%     
%     
%     
% end
rs=rand(px,1);
r=sort(rs);       %产生popsize个随机数并对其从小到大排序
j=1;
i=1;
while i<=px      %select the same size population from old pop
    if (r(i)<q(j))
        newpop3(i,:)=pop(j,:);
        i=i+1;  %累计概率大于随机数，选择该个体
    else
        j=j+1;  %否则，继续寻找适应度强的个体
    end 
end

    


end

//my_runGA.m
%main fun
clc;clear;
close all;

popsize = 50;
pc = 0.25;
pm = 0.01;

max_iter=1000;
Runtime=1;%运行1次

for i=1:Runtime
    [gy(i) ]= my_GA( popsize,pc,pm,max_iter );
end








