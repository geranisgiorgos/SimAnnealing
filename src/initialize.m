%εισαγωγή των διανυσμάτων p, d και r που εισήχθησαν από το excel
load matlab

% δημιουργία των πινάκων των περιορισμών
A = zeros(43+16, 43*16);
C = zeros(1, 43+16);

% περιορισμός Σxi = 1, για κάθε i = 1:43
for i=1:43
	A(i, ((i-1)*16+1):((i-1)*16+16)) = 1;
end
C(1:43) = 1;

% περιορισμός για τις πτυχιακές που μπορεί να δώσει κάθε καθηγητής
for i=44:43+16
	A(i, i:16:688) = 1;
end
C(44:43+16) = [4 2 2 2 2 2 2 2 2 4 2 3 4 4 3 3]; % από το αρχείο excel
C(44:43+16) = C(44:43+16) + 1;

% όλοι οι περιορισμοί είναι ισοτικοί
eqin = ones(1, 43);
eqin(44:43+16) = -1;

% δημιουργία τυχαίας αρχικής εκτίμησης
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

maxtries = 10; % μέγιστος αριθμός κλήσεων του αλγορίθμου

% κλήση του αλγόριθμου simulated annealing
f = 0;
x = [];
for i=1:maxtries
    [x_new f_new] = sm_anneal(p,ektimisi, A, C, eqin, 43, 16);
    if(f_new > f)
        x = x_new;
        f = f_new;
    end
end

% Εκτύπωση πίνακα αναθέσεων
disp('Ο τελικός πίνακας αναθέσεων: ');
x

x2 = bintprog(p,A,C);
% Εκτύπωση πίνακα αναθέσεων
disp('Ο τελικός πίνακας αναθέσεων: ');
x2