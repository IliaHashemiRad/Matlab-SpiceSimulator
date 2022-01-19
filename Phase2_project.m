clear; clc;

syms s t

fileName = 'E:\NetList.txt'; %file address
fileID = fopen(fileName);
tline = fgetl(fileID);

R_List = [];
R_NodeList = {};
R_ValueList = [];

C_List = [];
C_NodeList = {};
C_ValueList = [];

L_List = [];
L_NodeList = {};
L_ValueList = [];

V_List = [];
V_NodeList = {};
V_ValueList = {};

I_List = [];
I_NodeList = {};
I_ValueList = {};

VCVS_List = [];
VCVS_NodeList = {};
VCVS_CtrlNodeList = {};
VCVS_GValue = {};

VCCS_List = [];
VCCS_NodeList = {};
VCCS_CtrlNodeList = {};
VCCS_GValue = {};

CCCS_List = [];
CCCS_NodeList = {};
CCCS_CtrlNodeList = {};
CCCS_GValue = {};

CCVS_List = [];
CCVS_NodeList = {};
CCVS_CtrlNodeList = {};
CCVS_GValue = {};

ML1_List = [];
ML1_NodeList = {};
ML1_ValueList = [];
ML2_List = [];
ML2_NodeList = {};
ML2_ValueList = [];
ML_MValueList = [];

Nodes = [];

while ischar(tline)
    CurrentLine = split(tline, ',');
    Nodes(end + 1) = str2double(CurrentLine(3));
    Nodes(end + 1) = str2double(CurrentLine(4));
    if tline(1) == 'R'
        R_List{end + 1} = CurrentLine(2);
        R_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        R_ValueList(end + 1) = str2double(CurrentLine(5));
    end
    if tline(1) == 'C'
        C_List{end + 1} = CurrentLine(2);
        C_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        C_ValueList(end + 1) = str2double(CurrentLine(5));
    end
    if tline(1) == 'L'
        L_List{end + 1} = CurrentLine(2);
        L_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        L_ValueList(end + 1) = str2double(CurrentLine(5));
    end
    if tline(1) == 'V'
        V_List{end + 1} = CurrentLine(2);
        V_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        V_ValueList{end + 1} = CurrentLine(5);
    end
    if tline(1) == 'I'
        I_List{end + 1} = CurrentLine(2);
        I_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        I_ValueList{end + 1} = CurrentLine(5);
    end
    if tline(1) == 'Z'
        VCVS_List{end + 1} = CurrentLine(2);
        VCVS_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        VCVS_CtrlNodeList{end + 1} = [str2double(CurrentLine(5)), str2double(CurrentLine(6))];
        VCVS_GValue{end + 1} = CurrentLine(7);
    end
    if tline(1) == 'H'
        VCCS_List{end + 1} = CurrentLine(2);
        VCCS_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        VCCS_CtrlNodeList{end + 1} = [str2double(CurrentLine(5)), str2double(CurrentLine(6))];
        VCCS_GValue{end + 1} = CurrentLine(7);
    end 
    if tline(1) == 'Y'
        CCVS_List{end + 1} = CurrentLine(2);
        CCVS_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        CCVS_CtrlNodeList{end + 1} = [str2double(CurrentLine(5)), str2double(CurrentLine(6))];
        CCVS_GValue{end + 1} = CurrentLine(7);
    end
    if tline(1) == 'T'
        CCCS_List{end + 1} = CurrentLine(2);
        CCCS_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        CCCS_CtrlNodeList{end + 1} = [str2double(CurrentLine(5)), str2double(CurrentLine(6))];
        CCCS_GValue{end + 1} = CurrentLine(7);
    end
    if strcmp(CurrentLine(1), 'ML')
        ML1_List{end + 1} = CurrentLine(2);
        ML1_NodeList{end + 1} = [str2double(CurrentLine(3)), str2double(CurrentLine(4))];
        ML1_ValueList(end + 1) = str2double(CurrentLine(5));
        ML2_List{end + 1} = CurrentLine(6);
        ML2_NodeList{end + 1} = [str2double(CurrentLine(7)), str2double(CurrentLine(8))];
        ML2_ValueList(end + 1) = str2double(CurrentLine(9));
        ML_MValueList(end + 1) = str2double(CurrentLine(10));
    end
    
    tline = fgetl(fileID);
end
fclose(fileID);

V_ValueList = laplace(str2sym(string(V_ValueList)));
I_ValueList = laplace(str2sym(string(I_ValueList)));
VCVS_GValue = str2sym(string(VCVS_GValue));
VCCS_GValue = str2sym(string(VCCS_GValue));
CCCS_GValue = str2sym(string(CCCS_GValue));
CCVS_GValue = str2sym(string(CCVS_GValue));

% creating G matrix
NumNodes = length(unique(Nodes)) - 1;
G = sym(zeros(NumNodes));

for k = 1:NumNodes
    for j = 1:k
        if( k ~= j)
            index = cellfun(@(x) isequal(x, [k j]) || isequal(x, [j k]), R_NodeList);
            R = R_ValueList(index);
            g = -1./ R;
            G(k, j) = G(k, j) + sum(g,'all');

            % C
            index = cellfun(@(x) isequal(x, [k j]) || isequal(x, [j k]), C_NodeList);
            C = C_ValueList(index);
            Yc = -C*s;
            G(k, j) = G(k, j) + sum(Yc,'all');

            % L
            index = find(cellfun(@(x) isequal(x, [k j]) || isequal(x, [j k]), L_NodeList));
            L = L_ValueList(index);
            Yl = -1./(L*s);
            G(k, j) = G(k, j) + sum(Yl,'all');

            
        else
            index = cellfun(@(x) isequal(x(1), k) || isequal(x(2), k), R_NodeList);
            R = R_ValueList(index);
            g = 1./ R;
            G(k, j) = G(k, j) + sum(g,'all');

            % C
            index = cellfun(@(x) isequal(x(1), k) || isequal(x(2), k), C_NodeList);
            C = C_ValueList(index);
            Yc = C*s;
            G(k, j) = G(k, j) + sum(Yc,'all');

            % L
            index = cellfun(@(x) isequal(x(1), k) || isequal(x(2), k), L_NodeList);
            L = L_ValueList(index);
            Yl = 1./(L*s);
            G(k, j) = G(k, j) + sum(Yl,'all');

        end
        % for VCCS
        index1 = cellfun(@(x) isequal(x(1), k), VCCS_NodeList);
        index2 = cellfun(@(x) isequal(x(2), k), VCCS_NodeList);

        indexCtrl1 = cellfun(@(x) isequal(x(1), j), VCCS_CtrlNodeList);
        indexCtrl2 = cellfun(@(x) isequal(x(2), j), VCCS_CtrlNodeList);

        VCCS = VCCS_GValue .* (index1 .* indexCtrl1 + index2 .* indexCtrl2);
        G(k, j) = G(k, j) + sum(VCCS, 'all');

        VCCS = VCCS_GValue .* (index1 .* indexCtrl2 + index2 .* indexCtrl1);
        G(k, j) = G(k, j) - sum(VCCS, 'all');
        % done

        G(j, k) = G(k, j);
    end
end

% creating B matrix
N_V = length(V_List);
B1 = sym(zeros(NumNodes, N_V));
if(N_V > 0)
    for i = 1:N_V
        if V_NodeList{i}(1) ~= 0 && V_NodeList{i}(2) ~= 0
            B1(V_NodeList{i}(1), i) = 1;
            B1(V_NodeList{i}(2), i) = -1;
        elseif V_NodeList{i}(1) == 0
            B1(V_NodeList{i}(2), i) = -1;
        else
            B1(V_NodeList{i}(1), i) = 1;
        end
    end
end
%for VCVS
N_VCVS = length(VCVS_List);
B2 = sym(zeros(NumNodes, N_VCVS));
if(N_VCVS > 0)
    for i = 1:N_VCVS
        if VCVS_NodeList{i}(1) ~= 0 && VCVS_NodeList{i}(2) ~= 0
            B2(VCVS_NodeList{i}(1), i) = 1;
            B2(VCVS_NodeList{i}(2), i) = -1;
        elseif VCVS_NodeList{i}(1) == 0
            B2(VCVS_NodeList{i}(2), i) = -1;
        else
            B2(VCVS_NodeList{i}(1), i) = 1;
        end
    end
end

%for CCCS
N_CCCS = length(CCCS_List);
B3 = sym(zeros(NumNodes, N_CCCS));
if(N_CCCS > 0)
    for i = 1:N_CCCS
        if CCCS_CtrlNodeList{i}(1) ~= 0 && CCCS_CtrlNodeList{i}(2) ~= 0
            B3(CCCS_CtrlNodeList{i}(1), i) = 1;
            B3(CCCS_CtrlNodeList{i}(2), i) = -1;
        elseif CCCS_CtrlNodeList{i}(1) == 0
            B3(CCCS_CtrlNodeList{i}(2), i) = -1;
        else
            B3(CCCS_CtrlNodeList{i}(1), i) = 1;
        end

       
    end
end

B = [B1, B2, B3];

% creating C matrix
C = B.' ;
if N_VCVS > 0
    for k = 1:N_VCVS
        if VCVS_CtrlNodeList{k}(1) ~= 0 && VCVS_CtrlNodeList{k}(2) ~= 0
            C(k + N_V, VCVS_CtrlNodeList{k}(1)) = C(k + N_V, VCVS_CtrlNodeList{k}(1)) - VCVS_GValue(k);
            C(k + N_V, VCVS_CtrlNodeList{k}(2)) = C(k + N_V, VCVS_CtrlNodeList{k}(2)) + VCVS_GValue(k);
        elseif VCVS_CtrlNodeList{k}(1) == 0
            C(k + N_V, VCVS_CtrlNodeList{k}(2)) = C(k + N_V, VCVS_CtrlNodeList{k}(2)) + VCVS_GValue(k);
        else
            C(k + N_V, VCVS_CtrlNodeList{k}(1)) = C(k + N_V, VCVS_CtrlNodeList{k}(1)) - VCVS_GValue(k);
        end
    end
end

% updating B matrix for CCCS
if(N_CCCS > 0)
    for i = 1:N_CCCS
        if CCCS_NodeList{i}(1) ~= 0 && CCCS_NodeList{i}(2) ~= 0
            B(CCCS_NodeList{i}(1), N_V + N_VCVS + i) = B(CCCS_NodeList{i}(1), N_V + N_VCVS + i) + CCCS_GValue(i);
            B(CCCS_NodeList{i}(2), N_V + N_VCVS + i) = B(CCCS_NodeList{i}(2), N_V + N_VCVS + i) - CCCS_GValue(i);
        elseif CCCS_NodeList{i}(1) == 0
            B(CCCS_NodeList{i}(2), N_V + N_VCVS + i) = B(CCCS_NodeList{i}(2), N_V + N_VCVS + i) - CCCS_GValue(i);
        else
            B(CCCS_NodeList{i}(1), N_V + N_VCVS + i) = B(CCCS_NodeList{i}(1), N_V + N_VCVS + i) + CCCS_GValue(i);
        end
    end
end

% creating D matrix
D = sym(zeros(N_V + N_VCVS + N_CCCS));

% Creating A matrix
A = [G, B;
     C, D];

% creating i matrix
i = sym(zeros(NumNodes, 1));

if ~isempty(I_List)
    for k = 1:NumNodes
        index = cellfun(@(x) isequal(x(1), k), I_NodeList);
        I = I_ValueList .* index;
        i(k) = i(k) - sum(I, 'all');

        index = cellfun(@(x) isequal(x(2), k), I_NodeList);
        I = I_ValueList .* index;
        i(k) = i(k) + sum(I, 'all');
    end
end

% creating e matrix
e = V_ValueList.' ;
vcvs_new = sym(zeros(length(VCVS_List), 1));
cccs_new = sym(zeros(length(CCCS_List), 1));
e = [e;
     vcvs_new;
     cccs_new];

% creating z matrix
z = [i;
     e];

% for CCVS
N_CCVS = length(CCVS_List);
for n = 1:N_CCVS
    i = CCVS_CtrlNodeList{n}(1);
    j = CCVS_CtrlNodeList{n}(2);
    k = CCVS_NodeList{n}(1);
    l = CCVS_NodeList{n}(2);

    A = [A, sym(zeros(size(A, 1), 2));
         sym(zeros(2, size(A, 2) + 2))];

    if i ~= 0 && j ~= 0
        A(i, end - 1) = 1;
        A(j, end - 1) = -1;
        A(end - 1, i) = 1;
        A(end - 1, j) = -1;
    elseif i ==0
        A(j, end - 1) = -1;
        A(end - 1, j) = -1;
    else
        A(i, end - 1) = 1;
        A(end - 1, i) = 1;
    end

    if k ~= 0 && l ~= 0
        A(k, end) = 1;
        A(l, end) = -1;
        A(end, k) = 1;
        A(end, l) = -1;
    elseif k ==0
        A(l, end) = -1;
        A(end, l) = -1;
    else
       A(k, end) = 1;
       A(k, end) = 1;
    end

    A(end, end - 1) = -CCVS_GValue(n);

    z = [z;
         sym(zeros(2, 1))];

end

% for ML
N_ML = length(ML_MValueList);
for n = 1:N_ML
    i = ML1_NodeList{n}(1);
    j = ML1_NodeList{n}(2);
    k = ML2_NodeList{n}(1);
    l = ML2_NodeList{n}(2);
    
    A = [A, sym(zeros(size(A, 1), 2));
         sym(zeros(2, size(A, 2) + 2))];
    
    if i ~= 0 && j ~= 0
        A(i, end - 1) = 1;
        A(j, end - 1) = -1;
        A(end - 1, i) = 1;
        A(end - 1, j) = -1;
    elseif i ==0
        A(j, end - 1) = -1;
        A(end - 1, j) = -1;
    else
        A(i, end - 1) = 1;
        A(end - 1, i) = 1;
    end

    if k ~= 0 && l ~= 0
        A(k, end) = 1;
        A(l, end) = -1;
        A(end, k) = 1;
        A(end, l) = -1;
    elseif k ==0
        A(l, end) = -1;
        A(end, l) = -1;
    else
       A(k, end) = 1;
       A(k, end) = 1;
    end

    A(end - 1, end - 1) = -s*ML1_ValueList(n);
    A(end - 1, end) = -s*ML_MValueList(n);
    A(end, end - 1) = -s*ML_MValueList(n);
    A(end, end) = -s*ML2_ValueList(n);

    z = [z;
         sym(zeros(2, 1))];

end

% solving the final equation
X = linsolve(A, z);
X = simplify(simplifyFraction(X));
X = ilaplace(X);

fileName2 = 'Output5.txt';
fileID2 = fopen(fileName2,'W');

% Printing Outputs
for k = 1:length(R_List)
    if R_NodeList{k}(1) ~= 0 && R_NodeList{k}(2) ~= 0
        V = X(R_NodeList{k}(1)) - X(R_NodeList{k}(2));
    elseif R_NodeList{k}(2) == 0
        V = X(R_NodeList{k}(1));
    else
        V = -X(R_NodeList{k}(2));
    end
    I = V ./ R_ValueList(k);
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(R_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(R_List{k}));
end

for k = 1:length(C_List)
    if C_NodeList{k}(1) ~= 0 && C_NodeList{k}(2) ~= 0
        V = X(C_NodeList{k}(1)) - X(C_NodeList{k}(2));
    elseif C_NodeList{k}(2) == 0
        V = X(C_NodeList{k}(1));
    else
        V = -X(C_NodeList{k}(2));
    end
    I = ilaplace(laplace(V) .* s*C_ValueList(k));
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(C_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(C_List{k}));
end

for k = 1:length(L_List)
    if L_NodeList{k}(1) ~= 0 && L_NodeList{k}(2) ~= 0
        V = X(L_NodeList{k}(1)) - X(L_NodeList{k}(2));
    elseif L_NodeList{k}(2) == 0
        V = X(L_NodeList{k}(1));
    else
        V = -X(L_NodeList{k}(2));
    end
    I = ilaplace(laplace(V) ./ (s*L_ValueList(k)));
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(L_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(L_List{k}));
end

for k = 1:length(V_List)
    V = ilaplace(V_ValueList(k));
    I = X(k + NumNodes);
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(V_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(V_List{k}));
end

for k = 1:length(I_List)
    if I_NodeList{k}(1) ~= 0 && I_NodeList{k}(2) ~= 0
        V = X(I_NodeList{k}(1)) - X(I_NodeList{k}(2));
    elseif I_NodeList{k}(2) == 0
        V = X(I_NodeList{k}(1));
    else
        V = -X(I_NodeList{k}(2));
    end
    I = ilaplace(I_ValueList(k));
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(I_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(I_List{k}));
end

for k = 1:length(VCVS_List)
    if VCVS_NodeList{k}(1) ~= 0 && VCVS_NodeList{k}(2) ~= 0
        V = X(VCVS_NodeList{k}(1)) - X(VCVS_NodeList{k}(2));
    elseif VCVS_NodeList{k}(2) == 0
        V = X(VCVS_NodeList{k}(1));
    else
        V = -X(VCVS_NodeList{k}(2));
    end
    I = X(k + NumNodes + N_V);
    P = V .* I;
    
    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(VCVS_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(VCVS_List{k}));
end

for k = 1:length(VCCS_List)
    if VCCS_NodeList{k}(1) ~= 0 && VCCS_NodeList{k}(2) ~= 0
        V = X(VCCS_NodeList{k}(1)) - X(VCCS_NodeList{k}(2));
    elseif VCCS_NodeList{k}(2) == 0
        V = X(VCCS_NodeList{k}(1));
    else
        V = -X(VCCS_NodeList{k}(2));
    end
    
    if VCCS_CtrlNodeList{k}(1) ~= 0 && VCCS_CtrlNodeList{k}(2) ~= 0
        VCtrl = X(VCCS_CtrlNodeList{k}(1)) - X(VCCS_CtrlNodeList{k}(2));
    elseif VCCS_CtrlNodeList{k}(2) == 0
        VCtrl = X(VCCS_CtrlNodeList{k}(1));
    else
        VCtrl = -X(VCCS_CtrlNodeList{k}(2));
    end
    I = VCtrl .* VCCS_GValue(k);
    P = V .* I;
    
    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(VCCS_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(VCCS_List{k}));
end

for k = 1:N_CCCS
    if CCCS_NodeList{k}(1) ~= 0 && CCCS_NodeList{k}(2) ~= 0
        V = X(CCCS_NodeList{k}(1)) - X(CCCS_NodeList{k}(2));
    elseif CCCS_NodeList{k}(2) == 0
        V = X(CCCS_NodeList{k}(1));
    else
        V = -X(CCCS_NodeList{k}(2));
    end
    I = X(NumNodes + N_V + N_VCVS + k) .* CCCS_GValue(k);
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(CCCS_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(CCCS_List{k}));
end

for k = 1:N_CCVS
    if CCVS_NodeList{k}(1) ~= 0 && CCVS_NodeList{k}(2) ~= 0
        V = X(CCVS_NodeList{k}(1)) - X(CCVS_NodeList{k}(2));
    elseif CCVS_NodeList{k}(2) == 0
        V = X(CCVS_NodeList{k}(1));
    else
        V = -X(CCVS_NodeList{k}(2));
    end
    I = X(NumNodes + N_V + N_VCVS + N_CCCS + 2*k);
    P = V .* I;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(CCVS_List{k}), string(V), string(I), string(P));
    PlotPrinter(V, I, P, string(CCVS_List{k}));
end

for k = 1:length(ML_MValueList)
    if ML1_NodeList{k}(1) ~= 0 && ML1_NodeList{k}(2) ~= 0
        V1 = X(ML1_NodeList{k}(1)) - X(ML1_NodeList{k}(2));
    elseif ML1_NodeList{k}(2) == 0
        V1 = X(ML1_NodeList{k}(1));
    else
        V1 = -X(ML1_NodeList{k}(2));
    end
    I1 = X(NumNodes + N_V + N_VCVS + N_CCCS + 2*N_CCVS + 2*k-1);
    P1 = V1 .* I1;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(ML1_List{k}), string(V1), string(I1), string(P1));
    PlotPrinter(V1, I1, P1, string(ML1_List{k}));

    if ML2_NodeList{k}(1) ~= 0 && ML2_NodeList{k}(2) ~= 0
        V2 = X(ML2_NodeList{k}(1)) - X(ML2_NodeList{k}(2));
    elseif ML2_NodeList{k}(2) == 0
        V2 = X(ML2_NodeList{k}(1));
    else
        V2 = -X(ML2_NodeList{k}(2));
    end
    I2 = X(NumNodes + N_V + N_VCVS + N_CCCS + 2*N_CCVS + 2*k);
    P2 = V2 .* I2;

    fprintf(fileID2, '<%s><%s><%s><%s>\n', string(ML2_List{k}), string(V2), string(I2), string(P2));
    PlotPrinter(V2, I2, P2, string(ML2_List{k}));
end

fclose(fileID2);

%% functions
function PlotPrinter(V, I, P, Name)
    
    fig = figure('Visible','off');
    
    subplot(1,3,1)
    fplot(V,[0 50],'Color','b',LineWidth=1.5);
    title(sprintf('\n $V(t) \\ Plot \\ of \\ %s$ \n', Name) ,'Interpreter','latex',FontSize=10);
    xlabel('$t(s)$','Interpreter','latex',FontSize=8);
    ylabel('$V(v)$','Interpreter','latex',FontSize=8);

    subplot(1,3,2)
    fplot(I,[0 50],'Color','g',LineWidth=1.5);
    title(sprintf('\n $I(t) \\ Plot \\ of \\ %s$ \n', Name) ,'Interpreter','latex',FontSize=10);
    xlabel('$t(s)$','Interpreter','latex',FontSize=8);
    ylabel('$I(A)$','Interpreter','latex',FontSize=8);
    
    subplot(1,3,3)
    fplot(P,[0 50],'Color','r',LineWidth=1.5);
    title(sprintf('\n $P(t) \\ Plot \\ of \\ %s$ \n', Name) ,'Interpreter','latex',FontSize=10);
    xlabel('$t(s)$','Interpreter','latex',FontSize=8);
    ylabel('$P(w)$','Interpreter','latex',FontSize=8);

    exportgraphics(fig,'Out.pdf','Append',true);

end












