
function GA_TCP()
clear all;clc;close all;
%% 1.������ʼ��
    N=50;                     %���еĸ���
    M=200;                    %��Ⱥ�ĸ���
    C=500;                    %��������
    m=2;                      %��Ӧֵ��һ����̭����ָ��
    Pc=0.5;                   %�������
    Pmutation=0.2;            %�������
    pos=10*randn(N,2);        %���ɳ��е�����
    D=zeros(N,N);             %���ɳ���֮��������

%���������ľ���
    for i=1:N
        for j=i+1:N
            dis=(pos(i,1)-pos(j,1)).^2+(pos(i,2)-pos(j,2)).^2;
            D(i,j)=dis^(0.5);                                     
            D(j,i)=D(i,j);
        end
    end
    
%% 2.1���ɳ�ʼȺ��
    popm=zeros(M,N);
    for i=1:M
        popm(i,:)=randperm(N);             %��N����������
    end
    %%2.2���ѡ��һ����Ⱥ������Ӧ��·��
    R=popm(1,:);
    
    figure(1);
    scatter(pos(:,1),pos(:,2),'ro');
    plot_route(pos,R);                     %������Ⱥ������֮�������
    title('18��������������·��ͼ');
    str=['������','�ܾ���:',num2str(dis)];
    text(1000,5000,str);
    %% 3.��ʼ����Ⱥ������Ӧ����
    fitness=zeros(M,1);
    len=zeros(M,1);%%MΪ��Ⱥ�ĸ���
    for i=1:M
        len(i,1)=myLength(D,popm(i,:));
    end
    maxlen=max(len);
    minlen=min(len);
    fitness=fit(len,m,maxlen,minlen);
    rr=find(len==minlen);%%������len�����ҵ�����minlen�������±��¼��rr
    R=popm(rr(1,1),:);%%����Ӧ��Ⱥ���еı����¼��R
%     for i=1:N
%         fprintf('%d ',R(i));
%     end     %%�����Ӧ�ı���
%     fprintf('\n');
    fitness=fitness/sum(fitness);%%��Ӧ�Ⱥ������й�һ������
%%����CΪ�����Ĵ�����������������ڼ�¼ÿ�ε����õ�����Ⱥ�е���С����
    distance_min=zeros(C+1,1);  %%���ڴ洢���ε�������С����Ⱥ�ľ���
    while C>=0
%         fprintf('������%d��\n',500-C);
        %%ѡ�����%%
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
%         %%ÿ��ѡ�񶼱������ŵ���Ⱥ%%ѡ��������һ�������⣬û�дﵽɸѡ��Ч��
%         kk=;
%         popm_sel=[];
%         popm_sel=kk;

        [len_m len_index]=min(len_1);
        popm_sel=[popm_sel;popm(len_index,:)];%%���ǵ�����ԣ�������С��������û�б����棬�ҳ���С�����Ӧ�Ļ�����ӵ�popm_sel

        %popm_sel(1,:)=popm(len_index,:);

        %%�������
        nnper=randperm(nn);

        for i=1:nn*Pc
            A=popm_sel(nnper(i),:);
            B=popm_sel(nnper(i+1),:);
            [A,B]=cross(A,B);
            popm_sel(nnper(i),:)=A;
            popm_sel(nnper(i+1),:)=B;
        end
        %%�������
        for i=1:nn
            pick=rand;
            while pick==0
                pick=rand;
            end
            if pick<=Pmutation
                popm_sel(i,:)=Mutation(popm_sel(i,:));
            end
        end
        %%����Ӧ�Ⱥ���
        NN=size(popm_sel,1);
        len=zeros(NN,1);
        for i=1:NN
            len(i,1)=myLength(D,popm_sel(i,:));
        end
        maxlen=max(len);
        minlen=min(len);
        distance_min(501-C,1)=minlen;
        fitness=fit(len,m,maxlen,minlen);
        rr=find(len==minlen);%%�ҵ����Ż�����±�
        %fprintf('minlen=%d\n',minlen);
        R=popm_sel(rr(1,1),:);%%�ҵ����Ż���
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
    title('18����������������Ż�·��ͼ');
    str=['������','�ܾ���:',num2str(dis)];
    text(1000,5000,str);
    figure(3)
    plot(distance_min);
end
%% ���е������
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

%% ���������㺯��  mylength.m
function len=myLength(D,p)  %�����DΪ������֮��ľ������pΪ����˳���һ��һά����
    [N,NN]=size(D);
    len=D(p(1,N),p(1,1)); %��һ�ξ�����Ϊ��ʼֵ
    for i=1:(N-1)
        len=len+D(p(1,i),p(1,i+1));%%����˳�����Ӧ�ľ������,�൱�ڶ��ڵ�������˳��������
    end
end

%% ��Ӧ�Ⱥ���fit.m
function fitness=fit(len,m,maxlen,minlen)%%�����ֵ�ֱ���Ⱥ���в�ͬ�����Ӧ�ľ��룬��Ӧ�ȼ�����ָ̭����������ֵ����С����ֵ
    fitness=len;
    for i =1:length(len)
        fitness(i,1)=(1-(len(i,1)-minlen)/(maxlen-minlen+0.0001)).^m;%%Ϊ��ʹ���Ż��ٶȸ��죬����Ͱ�Ŀ�꺯��������ƽ��
    end
end

%% �����������  cross.m
function [A,B]=cross(A,B)
    L=length(A);
    if L<10
        W=L;
    elseif ((L/10)-floor(L/10))>=rand&&L>10
        W=ceil(L/10)+8;
    else
        W=floor(L/10)+8;
    end
    %%�����WΪ��Ҫ�����λ�������С��ʮ��ȫ�����棬����10����9������ʮ���Ϊ��ͬ���ʽ��棬ʹ�ý�����ڽϸߵĸ���
    p=unidrnd(L-W+1);%%�������������һ������ı�־λ��
    %fprintf('p=%d ',p);
    for i=1:W
        x=find(A==B(1,p+i-1));
        y=find(B==A(1,p+i-1));
        [A(1,p+i-1),B(1,p+i-1)]=exchange(A(1,p+i-1),B(1,p+i-1));
        [A(1,x),B(1,y)]=exchange(A(1,x),B(1,y));
    end
    %�����forѭ����ɽ���Ĺ���


end
%% �Ե����� exchange.m
function [x,y]=exchange(x,y)
    temp=x;
    x=y;
    y=temp;

end
%% ���캯�� Mutation.m
function a=Mutation(A)
    index1=0;
    index2=0;
    nnper=randperm(size(A,2));
    index1=nnper(1);
    index2=nnper(2);%%���������漴�±꣬Ȼ��ѡ��ǰ������Ϊ�±���л���
    %fprintf('index1=%d ',index1);
    %fprintf('index2=%d ',index2);

    temp=0;
    temp=A(index1);
    A(index1)=A(index2);
    A(index2)=temp;
    a=A;
end
