%�������� ��� ����������� p, d ��� r ��� ���������� ��� �� excel
load matlab

% ���������� ��� ������� ��� �����������
A = zeros(43+16, 43*16);
C = zeros(1, 43+16);

% ����������� �xi = 1, ��� ���� i = 1:43
for i=1:43
	A(i, ((i-1)*16+1):((i-1)*16+16)) = 1;
end
C(1:43) = 1;

% ����������� ��� ��� ��������� ��� ������ �� ����� ���� ���������
for i=44:43+16
	A(i, i:16:688) = 1;
end
C(44:43+16) = [4 2 2 2 2 2 2 2 2 4 2 3 4 4 3 3]; % ��� �� ������ excel
C(44:43+16) = C(44:43+16) + 1;

% ���� �� ����������� ����� ��������
eqin = ones(1, 43);
eqin(44:43+16) = -1;

% ���������� ������� ������� ���������
dosmenes = [4 2 2 2 2 2 2 2 2 4 2 3 4 4 3 3];
kathigites = 1:16;
ektimisi = zeros(43, 16);
for i=1:43
	temp = ceil(rand(1)*length(kathigites));
    ektimisi(i, kathigites(temp)) = 1;
    dosmenes(kathigites(temp)) = dosmenes(kathigites(temp)) - 1;
    if dosmenes(kathigites(temp)) == 0
        x=find(kathigites==kathigites(temp));
        kathigites=[kathigites(1:x-1) kathigites(x+1:length(kathigites))];
    end
end

maxtries = 10; % �������� ������� ������� ��� ����������

% ����� ��� ���������� simulated annealing
f = 0;
x = [];
for i=1:maxtries
    [x_new f_new] = sm_anneal(p,ektimisi, A, C, eqin, 43, 16);
    if(f_new > f)
        x = x_new;
        f = f_new;
    end
end

% �������� ������ ���������
disp('� ������� ������� ���������: ');
x

x2 = bintprog(p,A,C);
% �������� ������ ���������
disp('� ������� ������� ���������: ');
x2