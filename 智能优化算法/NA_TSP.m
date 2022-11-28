
function GA_TCP()
clear all;clc;close all;
%% 1.参数初始化
    N=50;                     %城市的个数
    M=200;                    %种群的个数
    C=500;                    %迭代次数
    m=2;                      %适应值归一化淘汰加速指数
    Pc=0.5;                   %交叉概率
    Pmutation=0.2;            %变异概率
    pos=10*randn(N,2);        %生成城市的坐标
    D=zeros(N,N);             %生成城市之间距离矩阵

%计算两点间的距离
    for i=1:N
        for j=i+1:N
            dis=(pos(i,1)-pos(j,1)).^2+(pos(i,2)-pos(j,2)).^2;
            D(i,j)=dis^(0.5);                                     
            D(j,i)=D(i,j);
        end
    end
    
%% 2.1生成初始群体
    popm=zeros(M,N);
    for i=1:M
        popm(i,:)=randperm(N);             %将N序号随机打乱
    end
    %%2.2随机选择一个种群画出相应的路径
    R=popm(1,:);
    
    figure(1);
    scatter(pos(:,1),pos(:,2),'ro');
    plot_route(pos,R);                     %画出种群各城市之间的连线
    title('18个城市坐标的最初路线图');
    str=['初代，','总距离:',num2str(dis)];
    text(1000,5000,str);
    %% 3.初始化种群及其适应函数
    fitness=zeros(M,1);
    len=zeros(M,1);%%M为种群的个数
    for i=1:M
        len(i,1)=myLength(D,popm(i,:));
    end
    maxlen=max(len);
    minlen=min(len);
    fitness=fit(len,m,maxlen,minlen);
    rr=find(len==minlen);%%在数组len当中找到等于minlen的数组下标记录进rr
    R=popm(rr(1,1),:);%%将相应种群当中的编码记录进R
%     for i=1:N
%         fprintf('%d ',R(i));
%     end     %%输出对应的编码
%     fprintf('\n');
    fitness=fitness/sum(fitness);%%适应度函数进行归一化处理
%%这里C为迭代的次数，这里的数组用于记录每次迭代得到的种群中的最小距离
    distance_min=zeros(C+1,1);  %%用于存储各次迭代的最小的种群的距离
    while C>=0
%         fprintf('迭代第%d次\n',500-C);
        %%选择操作%%
        nn=0;
        popm_sel=[];
        
        for i=1:size(popm,1)
            len_1(i,1)=myLength(D,popm(i,:));
        
            jc=rand*0.8;
            for j=1:size(popm,1)
                if fitness(j,1)>=jc
                    nn=nn+1;
                    popm_sel(nn,:)=popm(j,:);
                    break;
                end
            end
        end
%         %%每次选择都保存最优的种群%%选择函数存在一定的问题，没有达到筛选的效果
%         kk=;
%         popm_sel=[];
%         popm_sel=kk;

        [len_m len_index]=min(len_1);
        popm_sel=[popm_sel;popm(len_index,:)];%%考虑到随机性，避免最小距离的情况没有被保存，找出最小距离对应的基因添加到popm_sel

        %popm_sel(1,:)=popm(len_index,:);

        %%交叉操作
        nnper=randperm(nn);

        for i=1:nn*Pc
            A=popm_sel(nnper(i),:);
            B=popm_sel(nnper(i+1),:);
            [A,B]=cross(A,B);
            popm_sel(nnper(i),:)=A;
            popm_sel(nnper(i+1),:)=B;
        end
        %%变异操作
        for i=1:nn
            pick=rand;
            while pick==0
                pick=rand;
            end
            if pick<=Pmutation
                popm_sel(i,:)=Mutation(popm_sel(i,:));
            end
        end
        %%求适应度函数
        NN=size(popm_sel,1);
        len=zeros(NN,1);
        for i=1:NN
            len(i,1)=myLength(D,popm_sel(i,:));
        end
        maxlen=max(len);
        minlen=min(len);
        distance_min(501-C,1)=minlen;
        fitness=fit(len,m,maxlen,minlen);
        rr=find(len==minlen);%%找到最优基因的下标
        %fprintf('minlen=%d\n',minlen);
        R=popm_sel(rr(1,1),:);%%找到最优基因
%         for i=1:N
%             fprintf('%d ',R(i));
%         end
%         fprintf('\n');
        popm=[];
        
        popm=popm_sel;
        C=C-1;
        %pause(1);

        %size(popm)
    end

    figure(2)
    plot_route(pos,R);
    title('18个城市坐标的最终优化路线图');
    str=['初代，','总距离:',num2str(dis)];
    text(1000,5000,str);
    figure(3)
    plot(distance_min);
end
%% 城市点间连线
function plot_route(a,R)
    scatter(a(:,1),a(:,2),'rx');
    hold on;
    plot([a(R(1),1),a(R(length(R)),1)],[a(R(1),2),a(R(length(R)),2)]);
    hold on;
    for i=2:length(R)
        x0=a(R(i-1),1);
        y0=a(R(i-1),2);
        x1=a(R(i),1);
        y1=a(R(i),2);
        xx=[x0,x1];
        yy=[y0,y1];
        plot(xx,yy);
        hold on;
    end
end

%% 个体距离计算函数  mylength.m
function len=myLength(D,p)  %输入的D为各个点之间的距离矩阵，p为游览顺序的一个一维数组
    [N,NN]=size(D);
    len=D(p(1,N),p(1,1)); %第一段距离作为初始值
    for i=1:(N-1)
        len=len+D(p(1,i),p(1,i+1));%%按照顺序把相应的距离相加,相当于对于单个旅行顺序求解距离
    end
end

%% 适应度函数fit.m
function fitness=fit(len,m,maxlen,minlen)%%输入的值分别种群当中不同个体对应的距离，适应度加速淘汰指数，最大距离值和最小距离值
    fitness=len;
    for i =1:length(len)
        fitness(i,1)=(1-(len(i,1)-minlen)/(maxlen-minlen+0.0001)).^m;%%为了使得优化速度更快，这里就把目标函数故意求平方
    end
end

%% 交叉操作函数  cross.m
function [A,B]=cross(A,B)
    L=length(A);
    if L<10
        W=L;
    elseif ((L/10)-floor(L/10))>=rand&&L>10
        W=ceil(L/10)+8;
    else
        W=floor(L/10)+8;
    end
    %%这里的W为需要交叉的位数，如果小于十则全部交叉，等于10交叉9，大于十则分为不同概率交叉，使得交叉存在较高的概率
    p=unidrnd(L-W+1);%%这里是随机产生一个交叉的标志位置
    %fprintf('p=%d ',p);
    for i=1:W
        x=find(A==B(1,p+i-1));
        y=find(B==A(1,p+i-1));
        [A(1,p+i-1),B(1,p+i-1)]=exchange(A(1,p+i-1),B(1,p+i-1));
        [A(1,x),B(1,y)]=exchange(A(1,x),B(1,y));
    end
    %这里的for循环完成交叉的过程


end
%% 对调函数 exchange.m
function [x,y]=exchange(x,y)
    temp=x;
    x=y;
    y=temp;

end
%% 变异函数 Mutation.m
function a=Mutation(A)
    index1=0;
    index2=0;
    nnper=randperm(size(A,2));
    index1=nnper(1);
    index2=nnper(2);%%产生两个随即下标，然后选择前两个作为下标进行互换
    %fprintf('index1=%d ',index1);
    %fprintf('index2=%d ',index2);

    temp=0;
    temp=A(index1);
    A(index1)=A(index2);
    A(index2)=temp;
    a=A;
end
