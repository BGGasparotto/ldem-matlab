function ini = ANSYS_Copia(IFILE,CONT1,CONT2,CONT3,CONT4,CONT5,CONT6)
% ANSYS_Copia - Exporta a malha LDEM para Mechanical APDL.
%
% TYPE 1: temporario, durante a criacao dos elementos
% TYPE 2: aco (amarelo)
% TYPE 3: matriz epoxi (verde)
% TYPE 4: fibras nas tres camadas (vermelho)
% TYPE 5: barras de regioes enfraquecidas (cinza)
% TYPE 6: fissura no aco (azul)
%
% A classificacao por TYPE e visual. Os valores constitutivos de EP, E1 e
% E2 do LDEM NAO sao transferidos automaticamente para os materiais APDL.

% Manter os globais usados pelo programa principal.
global INT INB U V W U0 V0 W0 BDX BDY BDZ VLX VLY VLZ
global ACX ACY ACZ FRX FRY FRZ FX FY FZ MAS MSS CN ENR
global EDG E1 E2 LI DEQB F YNGB GFRB LCRB EPRB EP SS EM
global EN KR KR1 DSS DM FS CAF DT GFR DEQ LCR TEST EST
global DF DF1 DF2 DMP1 DMP2 TDAF CON M N L M1 N1 L1 MNL
global NNO NNS NNT NBA NBS NBT INCO NCON STEP NSTEP NDAMF
global RUN NRUN LCO TIME TTOT PSN ROH YNG EPR EPMAX EPMIN
global MEAN COEV ENEX ENCN ENEL ENGD ENDP ENIN DEN DERIV
global XB YB ZB UBAR VBAR WBAR REACAOX REACAOY REACAOZ
global STRAIN STRESS DM DS

% ========================================================================
% 1. PREPARACAO DO ARQUIVO APDL
% ========================================================================
assert(IFILE > 0, 'Identificador de arquivo APDL invalido.');
assert(LCO > 0 && NBT > 0, 'LCO e NBT devem ser positivos.');
assert(numel(CN) >= 2*NBT, 'Vetor CN menor do que 2*NBT.');
assert(numel(XB) >= NBT && numel(YB) >= NBT && numel(ZB) >= NBT, ...
       'Coordenadas de centros incompletas.');

fprintf(IFILE, '/PREP7\n');
fprintf(IFILE, '/TITLE, TEMPO = %g\n', TIME);

% Mantidas as mesmas definicoes do arquivo original.
for tipo = 1:6
    fprintf(IFILE, 'ET,%d,4\n', tipo);     % LINK4
end
fprintf(IFILE, 'MP,EX,1,2E10\n');
fprintf(IFILE, 'MP,PRXY,1,0.25\n');
fprintf(IFILE, 'MP,DENS,1,2500\n');
fprintf(IFILE, 'R,1,0.13,0.0014,0.36\n\n');

% ========================================================================
% 2. NOS E CONECTIVIDADES
% ========================================================================
fprintf(IFILE, '! COORDENADAS NODAIS\n');
for I = 1:NNT
    fprintf(IFILE, 'N,%d,%g,%g,%g\n', I, ...
            single(U(3,I)), single(V(3,I)), single(W(3,I)));
end
fprintf(IFILE, '\n! CONECTIVIDADES DOS ELEMENTOS\n');

% Contagem do estado original das barras (independe do TYPE geometrico).
CONT1 = 0; CONT2 = 0; CONT3 = 0;
CONT4 = 0; CONT5 = 0; CONT6 = 0;

% A ordem de criacao e preservada: elemento I corresponde a barra I.
fprintf(IFILE, 'TYPE,1\n');
for I = 1:NBT
    if DM(I)
        CONT1 = CONT1 + 1;
    elseif FS(I)
        CONT2 = CONT2 + 1;
    else
        CONT5 = CONT5 + 1;
    end
    fprintf(IFILE, 'E,%d,%d,!EL%d\n', CN(2*I-1), CN(2*I), I);
end
fprintf(IFILE, '\n');

% ========================================================================
% 3. VETORES DE CENTROS E COTAS DOS NOS: SEMPRE COLUNAS DE NBT POSICOES
% ========================================================================
xb = XB(1:NBT); xb = xb(:);
yb = YB(1:NBT); yb = yb(:);
zb = ZB(1:NBT); zb = zb(:);

NOi = CN(1:2:2*NBT); NOi = NOi(:);
NOf = CN(2:2:2*NBT); NOf = NOf(:);
zi = W0(NOi); zi = zi(:);
zf = W0(NOf); zf = zf(:);

tol = 1e-9;                   % tolerancia em metros
cruzaZ = abs(zi-zf) > tol;   % barras com componente na direcao Z

% ========================================================================
% 4. GEOMETRIA DOS TRES DEGRAUS - MANTIDAS AS MASCARAS DO MAIN1.M
% ========================================================================
% Regiao de Z aproximadamente 2 a 3*LCO
Logic1  = (zb <= 2.99*LCO) & (zb >= 1.99*LCO);
Logic2  = (yb < 46.99*LCO);
Logic3  = (yb > (47.01+82)*LCO);
Logic4  = (xb < 10*LCO);
Logic5  = (xb > 90*LCO);
Logic17 = (abs(zb-3.0*LCO) < tol) & cruzaZ;

Mask1  = Logic1 & Logic2;
Mask2  = Logic1 & Logic3;
Mask3  = Logic1 & Logic4;
Mask4  = Logic1 & Logic5;
Mask14 = Logic17 & Logic2;
Mask15 = Logic17 & Logic3;
Mask16 = Logic17 & Logic4;
Mask17 = Logic17 & Logic5;

MaskReducao1 = Mask1 | Mask2 | Mask3 | Mask4 | ...
               Mask14 | Mask15 | Mask16 | Mask17;

% Regiao de Z aproximadamente 1 a 2*LCO
Logic6  = (zb < 1.99*LCO) & (zb > 0.99*LCO);
Logic7  = (yb < 56.99*LCO);
Logic8  = (yb > (56.99+62)*LCO);
Logic9  = (xb < 20*LCO);
Logic10 = (xb > 80*LCO);
Logic18 = (abs(zb-2.0*LCO) < tol) & cruzaZ;

Mask5  = Logic6 & Logic7;
Mask6  = Logic6 & Logic8;
Mask7  = Logic6 & Logic9;
Mask8  = Logic6 & Logic10;
Mask18 = Logic18 & Logic7;
Mask19 = Logic18 & Logic8;
Mask20 = Logic18 & Logic9;
Mask21 = Logic18 & Logic10;

MaskReducao2 = MaskReducao1 | Mask5 | Mask6 | Mask7 | Mask8 | ...
               Mask18 | Mask19 | Mask20 | Mask21;

% Regiao de Z aproximadamente 0 a 1*LCO
Logic11 = (zb <= 0.99*LCO);
Logic12 = (yb < 66.99*LCO);
Logic13 = (yb > (66.99+42)*LCO);
Logic14 = (xb < 30*LCO);
Logic15 = (xb > 70*LCO);
Logic19 = (abs(zb-1.0*LCO) < tol) & cruzaZ;

Mask9  = Logic11 & Logic12;
Mask10 = Logic11 & Logic13;
Mask11 = Logic11 & Logic14;
Mask12 = Logic11 & Logic15;
Mask22 = Logic19 & Logic12;
Mask23 = Logic19 & Logic13;
Mask24 = Logic19 & Logic14;
Mask25 = Logic19 & Logic15;

MaskReducao = MaskReducao2 | Mask9 | Mask10 | Mask11 | Mask12 | ...
              Mask22 | Mask23 | Mask24 | Mask25;

% ========================================================================
% 5. FISSURA NO ACO: 39 mm DE COMPRIMENTO, 45 GRAUS, CENTRO (50,88) mm
% ========================================================================
% As coordenadas do programa estao em metros; LCO = 0.001 m (1 mm).
a = 39*LCO;          % comprimento total em metros
theta = 45;          % angulo em graus
xg = 50*LCO;         % centro em X
yg = 88*LCO;         % centro em Y

xmin = xg - (a/2)*cosd(theta);
xmax = xg + (a/2)*cosd(theta);
ymin = yg - (a/2)*sind(theta);
ymax = yg + (a/2)*sind(theta); %#ok<NASGU> % para conferencia

% Reta central da fissura: y = ymin + (x-xmin)*tan(theta).
% Faixa vertical de +/- 1*LCO ao redor dessa reta (largura 2 mm em Y).
yReta = ymin + (xb-xmin)*tand(theta);
Logic20 = (zb >= 3*LCO);                 % somente aco
Logic21 = (xb > xmin) & (xb < xmax);    % limites longitudinais
Logic22 = (yb > (yReta-LCO)) & (yb < (yReta+LCO));
MaskFissuraGeom = Logic20 & Logic21 & Logic22;

% ========================================================================
% 6. FIBRAS NO REPARO, EM Z/LCO = 0.5, 1.5 E 2.5
% ========================================================================
tolRel = 1e-8;
x = xb/LCO; y = yb/LCO; z = zb/LCO;

MaskPlano05 = abs(z-0.5) < tolRel;
MaskPlano15 = abs(z-1.5) < tolRel;
MaskPlano25 = abs(z-2.5) < tolRel;
MaskPlanosFibra = MaskPlano05 | MaskPlano15 | MaskPlano25;

xInteiro = abs(x-round(x)) < tolRel;
yInteiro = abs(y-round(y)) < tolRel;
xMeio = abs(x-(floor(x)+0.5)) < tolRel;
yMeio = abs(y-(floor(y)+0.5)) < tolRel;

% Fibras em X: x semi-inteiro, y inteiro;
% fibras em Y: x inteiro, y semi-inteiro.
MaskFibraX = MaskPlanosFibra & xMeio & yInteiro;
MaskFibraY = MaskPlanosFibra & xInteiro & yMeio;

% ========================================================================
% 7. CLASSIFICACAO EXCLUSIVA, POR PRIORIDADE
% ========================================================================
MaskEliminadas = MaskReducao;
MaskFibras = (MaskFibraX | MaskFibraY) & ~MaskEliminadas;
MaskFissura = MaskFissuraGeom & ~MaskEliminadas & ~MaskFibras;
MaskAco = (zb >= 3*LCO) & ~MaskEliminadas & ~MaskFibras & ~MaskFissura;
MaskMatriz = (zb < 3*LCO) & ~MaskEliminadas & ~MaskFibras & ~MaskFissura;

% Cada barra deve ser classificada em uma e somente uma categoria.
contagemTipos = uint8(MaskEliminadas) + uint8(MaskFibras) + ...
                uint8(MaskFissura) + uint8(MaskAco) + uint8(MaskMatriz);
assert(numel(contagemTipos)==NBT && all(contagemTipos==1), ...
       'Ha barras sem TYPE ou com TYPE duplicado.');

TIPOBARRA = zeros(NBT,1,'uint8');
TIPOBARRA(MaskAco)        = 2;
TIPOBARRA(MaskMatriz)     = 3;
TIPOBARRA(MaskFibras)     = 4;
TIPOBARRA(MaskEliminadas) = 5;
TIPOBARRA(MaskFissura)    = 6;

fprintf('\n===== CONFERENCIA DOS TIPOS =====\n');
for tipo = 1:6
    fprintf('TYPE %d: %d elementos\n', tipo, nnz(TIPOBARRA==tipo));
end
fprintf('Fibras X: %d barras; fibras Y: %d barras\n', ...
        nnz(MaskFibras & MaskFibraX), nnz(MaskFibras & MaskFibraY));
fprintf('Fibras Z=0.5: %d; Z=1.5: %d; Z=2.5: %d\n', ...
        nnz(MaskFibras & MaskPlano05), ...
        nnz(MaskFibras & MaskPlano15), ...
        nnz(MaskFibras & MaskPlano25));
fprintf('Fissura geometrica, antes de prioridades: %d barras\n', ...
        nnz(MaskFissuraGeom));
fprintf('Fissura final (TYPE 6): %d barras\n', nnz(MaskFissura));
fprintf('Fissura sobreposta a regioes enfraquecidas: %d barras\n', ...
        nnz(MaskFissuraGeom & MaskEliminadas));
fprintf('Extremos da fissura (mm): X=[%.3f, %.3f], Y=[%.3f, %.3f]\n', ...
        xmin/LCO, xmax/LCO, ymin/LCO, ymax/LCO);
if ~any(MaskFissura)
    warning('TYPE 6 sem barras. Confira os centros XB, YB, ZB e LCO.');
end

% ========================================================================
% 8. MODIFICACAO DOS TIPOS DOS ELEMENTOS NO ANSYS
% ========================================================================
fprintf(IFILE, '! CLASSIFICACAO POR TIPOS\n');
for I = 1:NBT
    fprintf(IFILE, 'EMODIF,%d,TYPE,%d\n', I, TIPOBARRA(I));
end

% Cores por TYPE: 2 amarelo, 3 verde, 4 vermelho, 5 cinza, 6 azul.
fprintf(IFILE, '/COLOR,NUM,YELL,2\n');
fprintf(IFILE, '/COLOR,NUM,GREE,3\n');
fprintf(IFILE, '/COLOR,NUM,RED,4\n');
fprintf(IFILE, '/COLOR,NUM,LGRA,5\n');
fprintf(IFILE, '/COLOR,NUM,BLUE,6\n');
fprintf(IFILE, '/PNUM,TYPE,1\n');
fprintf(IFILE, '/NUMBER,1\n');
fprintf(IFILE, 'ALLSEL,ALL\n\n');

% ========================================================================
% 9. CONDICOES DE CONTORNO E CARGAS: PRESERVADAS DO ARQUIVO ORIGINAL
% ========================================================================
fprintf(IFILE, '! CONDICOES DE CONTORNO\n');
if TEST == 102
    for I = 1:NNT
        if BDX(I) == 1
            fprintf(IFILE, 'D,%d,UX\n', I);
        end
        if BDY(I) == 1
            fprintf(IFILE, 'D,%d,UY\n', I);
        end
        if BDZ(I) == 1
            fprintf(IFILE, 'D,%d,UZ\n', I);
        end
    end
end

fprintf(IFILE, '\n! CARGAS APLICADAS\n');
if TEST == 102
    fprintf(IFILE, 'F,%d,FY,%d\n', 17664, 100);
end

% ========================================================================
% 10. VISTA E INFORMACOES DE ESTADO
% ========================================================================
fprintf(IFILE, '\n/VIEW,1,0,-1,0\n');
fprintf(IFILE, '/PBC,ALL,1\n');
fprintf(IFILE, 'EPLOT\n\n');

fprintf(IFILE, '! ESTADO DAS BARRAS NO LDEM\n');
fprintf(IFILE, '! BARRAS INTACTAS = %d\n', CONT1+CONT3+CONT6);
fprintf(IFILE, '! BARRAS DANIFICADAS = %d\n', CONT2+CONT4);
fprintf(IFILE, '! BARRAS ROMPIDAS = %d\n\n', CONT5);

% Arquivo-resumo mantido para compatibilidade.
VF = 1.0;
fid30 = fopen('BARR_ROMP.dat','w');
if fid30 < 0
    warning('Nao foi possivel criar BARR_ROMP.dat');
else
    fprintf(fid30, '%12.6E   %12.6E   %7d   %7d   %7d\n', ...
            TIME, VF, CONT1, CONT2, CONT5);
    fclose(fid30);
end

ini = 1;
end