clear all
clc
global INT INB U V W U0 V0 W0 BDX BDY BDZ VLX VLY VLZ 
global ACX ACY ACZ FRX FRY FRZ FX FY FZ MAS MSS CN ENR
global EDG E1 E2 LI DEQB F YNGB GFRB LCRB EPRB EP SS EM
global EN KR KR1 DSS DM FS CAF DT GFR DEQ LCR TEST EST 
global DF DF1 DF2 DMP1 DMP2 TDAF CON M N L M1 N1 L1 MNL
global NNO NNS NNT NBA NBS NBT INCO NCON STEP NSTEP NDAMF
global RUN NRUN LCO TIME TTOT PSN ROH YNG EPR EPMAX EPMIN 
global MEAN  COEV ENEX ENCN ENEL ENGD ENDP ENIN DEN DERIV 
global XB YB ZB UBAR VBAR WBAR REACAOX REACAOY REACAOZ 
global STRAIN STRESS SEED 
global OU2 OU3 OU4 OU5 OU6 OU7 OU8 OU9
global fid_Lei
%
ini=ini_DEM3()  %an initialization of the main variables was done
%
DATNAME=' ';
CONNAME=' ';
ANSNAME=' ';
RESNAME=' ';
DISNAME=' ';
POSNAME=' ';
ENENAME=' ';
%ARMANS(1000)='              ';
CONT1=0;
CONT2=0;
CONT3=0;
CONT4=0;
CONT5=0;
CONT6=0;
%counter of the linked the file writing
PR3=1; %Ansys File
PR67=1; %DIS and POS  files
PR8=1; %ENE file
PR9=1; % SIM01 file inside CALRES
%
      % CHARACTER*14 DATNAME,CONNAME
      % CHARACTER*14 ANSNAME,RESNAME,DISNAME
      % CHARACTER*14 POSNAME,ENENAME
      % CHARACTER*14 ARMANS(10000000)
      % integer cont1,cont2,cont3,cont4,cont5,cont6,B,bb
%
%
% !
%       OPEN (UNIT=10,ACCESS='SEQUENTIAL',FILE='DATA0.DAT',&
%      &     STATUS='OLD')
%   
% !
% !     READAT -----> Lê os dados que estão no arquivo DATA0.DAT
      % CALL READAT0
      ini=read_dat0() % read the input data saved in DATA0.DAT
      'ini=read_dat0()'
      %
      % !ACX
%     Close (10)
% !
% !--------------------------------------------------------------------------
% !  PREPARAÇÃO DAS VARIÁVEIS E CHAMAMENTO DAS SUBROTINAS PRELIMINARES
% !--------------------------------------------------------------------------
% !
% 
% !     PREPVAR ----> setting of initial parameters
% CALL PREPVAR 
 ini=prepvar()  %setting of initial parameters    
'ini=prepvar()'
 
 % !
% !     COORD ------> Gera coordenadas de um prisma
%       CALL COORD
ini=COORD() % Build the nodes coordinate of prism  with MxNxL modules 
'ini=COORD()'

% !
% 
% 
% !     IMPERF -----> Impõe uma aleatoriedade (imperfeição) às coordenadas
%       CALL IMPERF
ini=IMPERF()% !
'ini=IMPERF()'
% !     BOUND ------> Determina as condições de contorno da estrutura
%       CALL BOUND
ini=BOUND()% !
'ini=BOUND()'
% !
% !     NODMASS ----> Atribui um valor de massa a cada nó da estrutura
%       CALL NODMASS
ini=NODMASS()% !
'ini=NODMASS()'
% ========================================================================
% MASSAS NODAIS - ACO, MATRIZ E REGIOES ENFRAQUECIDAS (VAZIOS)
% ========================================================================
% A NODMASS() calculou inicialmente a massa de um prisma todo de aco.
% Aqui, reconstruimos as massas por celulas de dimensao LCO x LCO x LCO.
%
% Hipotese da discretizacao LDEM: cada celula possui 8 nos de vertices
% e 1 no central. Metade da massa fica no centro e metade e repartida
% igualmente entre os 8 vertices (mCel/16 para cada vertice).
%
% IMPORTANTE: a massa equivalente das fibras ainda NAO foi adicionada.
% As barras da fissura tambem nao representam remocao de volume aqui.
% ========================================================================

MassaOriginal = sum(MAS(1:NNT));

rhoAco   = ROH;                    % kg/m3 - aco
ROHm     = 1500;                   % kg/m3 - epoxi real
ROHmr    = (0.4/3)*ROHm;           % kg/m3 - epoxi equivalente
rhoVazio = 1.E-5*rhoAco;          % kg/m3 - massa residual numerica

% Numero de celulas (nao de nos)
nx = M-1;
ny = N-1;
nz = L-1;

assert(nx > 0 && ny > 0 && nz > 0, ...
    'Dimensoes da malha invalidas');
assert(NNT == (nx+1)*(ny+1)*(nz+1) + nx*ny*nz, ...
    'A malha nao corresponde a vertices + centros de celulas');
assert(numel(MAS) >= NNT, 'Vetor MAS menor que NNT');

% Centros das celulas, em unidades de LCO
[Xc,Yc,Zc] = ndgrid(0.5:(nx-0.5), ...
                     0.5:(ny-0.5), ...
                     0.5:(nz-0.5));

% Inicialmente todo o volume e tratado como vazio
rhoCel = rhoVazio*ones(nx,ny,nz);

% Placa continua de aco: Z = 3 a 6 mm
MaskAcoCel = (Zc >= 3);
rhoCel(MaskAcoCel) = rhoAco;

% Camada inferior do reparo: Z = 0 a 1 mm; X = 30..70; Y = 67..109
Regiao1Cel = (Zc < 1) & ...
             (Xc >= 30) & (Xc < 70) & ...
             (Yc >= 67) & (Yc < 109);

% Camada intermediaria: Z = 1 a 2 mm; X = 20..80; Y = 57..119
Regiao2Cel = (Zc >= 1) & (Zc < 2) & ...
             (Xc >= 20) & (Xc < 80) & ...
             (Yc >= 57) & (Yc < 119);

% Camada superior: Z = 2 a 3 mm; X = 10..90; Y = 47..129
Regiao3Cel = (Zc >= 2) & (Zc < 3) & ...
             (Xc >= 10) & (Xc < 90) & ...
             (Yc >= 47) & (Yc < 129);

MaskEpoxiCel = Regiao1Cel | Regiao2Cel | Regiao3Cel;
assert(~any(MaskAcoCel(:) & MaskEpoxiCel(:)), ...
    'Sobreposicao entre celulas de aco e epoxi');
rhoCel(MaskEpoxiCel) = ROHmr;

% Massa total da celula
mCel = rhoCel*LCO^3;

% Coordenadas nominais dos nos, normalizadas por LCO
xnM = U0(1:NNT); xnM = xnM(:)/LCO;
ynM = V0(1:NNT); ynM = ynM(:)/LCO;
znM = W0(1:NNT); znM = znM(:)/LCO;

tolRelMassa = 1e-6;

% Nos de vertices em (i, j, k)
MaskVertice = ...
    (abs(xnM-round(xnM)) < tolRelMassa) & ...
    (abs(ynM-round(ynM)) < tolRelMassa) & ...
    (abs(znM-round(znM)) < tolRelMassa);

% Nos centrais em (i+0.5, j+0.5, k+0.5)
MaskCentro = ...
    (abs(xnM-(floor(xnM)+0.5)) < tolRelMassa) & ...
    (abs(ynM-(floor(ynM)+0.5)) < tolRelMassa) & ...
    (abs(znM-(floor(znM)+0.5)) < tolRelMassa);

assert(all(MaskVertice | MaskCentro) && ...
       ~any(MaskVertice & MaskCentro), ...
    'Existem nos fora dos vertices/centros esperados');
assert(nnz(MaskCentro) == nx*ny*nz, ...
    'Numero de nos centrais incorreto');
assert(nnz(MaskVertice) == (nx+1)*(ny+1)*(nz+1), ...
    'Numero de nos de vertices incorreto');

% Construir mapa de identificadores dos vertices
idV = zeros(nx+1,ny+1,nz+1,'uint32');
nosV = find(MaskVertice);
indV = sub2ind(size(idV), ...
    round(xnM(nosV))+1, ...
    round(ynM(nosV))+1, ...
    round(znM(nosV))+1);
idV(indV) = uint32(nosV);

% Construir mapa de identificadores dos centros
idC = zeros(nx,ny,nz,'uint32');
nosC = find(MaskCentro);
indC = sub2ind(size(idC), ...
    floor(xnM(nosC))+1, ...
    floor(ynM(nosC))+1, ...
    floor(znM(nosC))+1);
idC(indC) = uint32(nosC);

assert(all(idV(:)>0) && all(idC(:)>0), ...
    'Mapa de identificadores dos nos incompleto');

% Montagem das novas massas nodais
MASnova = zeros(NNT,1);

% 50% da massa em cada no central
MASnova(double(idC(:))) = 0.5*mCel(:);

% Outros 50% igualmente divididos entre os oito vertices
for dx = 0:1
    for dy = 0:1
        for dz = 0:1
            ids = idV(1+dx:nx+dx, ...
                      1+dy:ny+dy, ...
                      1+dz:nz+dz);
            ids = double(ids(:));
            MASnova(ids) = MASnova(ids) + mCel(:)/16;
        end
    end
end

% Verificar conservacao da massa
MassaTeorica = sum(mCel(:));
MassaNova = sum(MASnova);

assert(all(isfinite(MASnova)) && all(MASnova > 0), ...
    'Massas nodais nao positivas ou invalidas');
assert(abs(MassaNova-MassaTeorica) < 1e-10, ...
    'Massa nodal e massa volumetrica nao coincidem');

% Aplicar massa corrigida (SEM aplicar a correcao nodal antiga)
MAS(1:NNT) = MASnova;

fprintf('\n===== MASSA CORRIGIDA POR CELULAS =====\n');
fprintf('Nos centrais: %d | Nos de vertices: %d\n', ...
        nnz(MaskCentro),nnz(MaskVertice));
fprintf('Celulas de aco: %d\n', nnz(MaskAcoCel));
fprintf('Celulas de epoxi: %d\n', nnz(MaskEpoxiCel));
fprintf('Celulas vazias: %d\n', ...
        numel(rhoCel)-nnz(MaskAcoCel)-nnz(MaskEpoxiCel));
fprintf('Massa original: %.9f kg\n', MassaOriginal);
fprintf('Massa teorica: %.9f kg\n', MassaTeorica);
fprintf('Massa nodal corrigida: %.9f kg\n', sum(MAS(1:NNT)));
fprintf('Massa retirada: %.9f kg\n', MassaOriginal-MassaNova);

% Evita manter grandes matrizes auxiliares na memoria para a integracao
clear Xc Yc Zc rhoCel mCel idV idC MASnova
clear xnM ynM znM MaskVertice MaskCentro
clear nosV nosC indV indC MaskAcoCel MaskEpoxiCel
clear Regiao1Cel Regiao2Cel Regiao3Cel


% !
% !     CONECTI ----> Gera conectividades das barras
%       CALL CONECTI
ini=CONECTI()% !
'ini=CONECTI()'
% !
% !     CALAUX -----> Cálculos auxiliares e determinação compr. cada barra
%       CALL CALAUX
ini=CALAUX()% !

'ini=CALAUX()'
%
% 
%      OPEN (UNIT=70,ACCESS='SEQUENTIAL',FILE='BARI.txt',&
%     &        STATUS='UNKNOWN')
% 
%     DO I=1, NBT
%       WRITE(70,'(I7,3X,3(E25.12,3X))') I,XB(I),YB(I),ZB(I)
% 
%     ENDDO
% 
%     CLOSE(70)
% Open the file for writing ('w' creates/overwrites, equivalent to STATUS='UNKNOWN')
fid = fopen('BARI.txt', 'w');

% Loop equivalent to DO I=1, NBT
for I = 1:NBT
    % fprintf replaces WRITE with formatted output
    % Format spec: %7d (I7), 3 spaces (3X), 3 floating-point values in scientific notation (E25.12)
    fprintf(fid, '%7d   %25.12E   %25.12E   %25.12E\n', I, XB(I), YB(I), ZB(I));
end


% Close the file
fclose(fid);
'The file BARI.txt was made'



% !   
% !     NAYFES0 ---->BAR STIFNESS DETERMINATION
%       CALL NAYFES0
ini=NAYFES0()
'ini=NAYFES0()'
% !     RMEC ------->BAR STIFNESS DETERMINATION 
%       CALL RMEC
ini=RMEC()
'ini=RMEC()'


% % % % % --- Introduce random field into EP(B) ----------------------------------
% % % % fid = fopen('arquivo_rodrigo', 'r');
% % % % if fid == -1
% % % %     error('Console1: cannot open arquivo_3lc_rnao_cv05');
% % % % %    !arquivo_lc3l_rnao_cv100.txt
% % % % %    !arquivo_3lc_rnao_cv05.txt      
% % % % end
% % % % for I = 1:NBT
% % % %     Phitemp = fscanf(fid, '%f', 1);
% % % %     EP(I) = Phitemp * EP(I);
% % % %     EM(I) = EP(I);            % very important!!!
% % % % end
% % % % fprintf('%g\n', EP(1));
% % % % fclose(fid);

ini=CALCOEV() % computation of Mean value and variation 
% coeficient of EP 
'ini=CALCOEV()' 
%  ini= NEIBOUR()   !WE NEED THIS SBROUTINE TO COMPUTE BARS 
% THAT ARRIVING IN THE NODS i AND J OF THE b BAR
% 'ini=CALCOEV()'  ( We need to implement)
% 
% !--------------------------------------------------------------------------
% !  ABRINDO OS ARQUIVOS DE SAÍDA
% !--------------------------------------------------------------------------
fid_CON = -1; fid_AN1 = -1; fid_RES = -1;
fid_DIS = -1; fid_POS = -1; fid_ENE = -1;
fid_Lei = -1; fid_EnB = -1; fid_RES =-1;

if OU2 > 0, fid_CON = fopen('LET.CON', 'w'); end
if OU3 > 0, fid_AN1 = fopen('LET.AN1', 'r'); end
if OU5 > 0, fid_RES = fopen('LET.RES', 'w'); end
if OU6 > 0, fid_DIS = fopen('LET.DIS', 'w'); end
if OU7 > 0, fid_POS = fopen('LET.POS', 'w'); end
if OU8 > 0, fid_ENE = fopen('LET.ENE', 'w'); end
if OU9 > 0, fid_Lei = fopen('Lei_Barra.TXT',   'w'); end
if OU9 > 0, fid_EnB = fopen('Energia_Barra.TXT','w'); end
if OU9 > 0, fid_RES = fopen('SIM01.DAT','w'); end

% Considerar somente as barras efetivamente geradas
xb = XB(1:NBT);
yb = YB(1:NBT);
zb = ZB(1:NBT);
xb = xb(:);
yb = yb(:);
zb = zb(:);

% Nós inicial e final das barras
NOi = CN(1:2:2*NBT);
NOf = CN(2:2:2*NBT);

NOi = NOi(:);
NOf = NOf(:);

% Coordenadas Z dos nos
zi = W0(NOi);
zf = W0(NOf);

% FORCAR VETORES COLUNA
zi = zi(:);
zf = zf(:);
tol = 1e-9;
cruzaZ = abs(zi - zf) > tol;


%Elasticidade epoxi
EMA=3.5E9;   %aplicação do módulo de elasticidade epoxi
EMAR=(0.4/3)*EMA;

Logic16=(zb<=(2.99*LCO));
Mask13=(Logic16);

E1(Mask13)=(EMAR/YNG)*E1(Mask13);
E2(Mask13)=(EMAR/YNG)*E2(Mask13);

%Reducao da rigidez
% Condições geométricas
Logic1=(zb<=(2.99*LCO)&(zb>=(1.99*LCO)));
Logic2=(yb<(46.99*LCO));
Logic3=(yb>((47.01+82)*LCO));
Logic4=(xb<(10*LCO));
Logic5=(xb>(90*LCO));
Logic17 = (abs(zb - 3.0*LCO) < tol) & cruzaZ;


Mask1=(Logic1 & Logic2);
Mask2=(Logic1 & Logic3);
Mask3=(Logic1 & Logic4);
Mask4=(Logic1 & Logic5);
Mask14=(Logic17 & Logic2);
Mask15=(Logic17 & Logic3);
Mask16=(Logic17 & Logic4);
Mask17=(Logic17 & Logic5);


MaskReducao1 = Mask1 | Mask2 | Mask3 | Mask4 | ...
              Mask14 | Mask15 | Mask16 | Mask17 ;


Logic6=((zb<(1.99*LCO))&(zb>(0.99*LCO)));
Logic7=(yb<(56.99*LCO));
Logic8=(yb>((56.99+62)*LCO));
Logic9=(xb<(20*LCO));
Logic10=(xb>(80*LCO));
Logic18 = (abs(zb - 2.0*LCO) < tol) & cruzaZ;


Mask5=(Logic6 & Logic7);
Mask6=(Logic6 & Logic8);
Mask7=(Logic6 & Logic9);
Mask8=(Logic6 & Logic10);
Mask18=(Logic18 & Logic7);
Mask19=(Logic18 & Logic8);
Mask20=(Logic18 & Logic9);
Mask21=(Logic18 & Logic10);

MaskReducao2 = MaskReducao1 | Mask5 | Mask6 | Mask7 | ...
              Mask8 | Mask18 | Mask19 | Mask20 | ...
              Mask21;

Logic11=(zb<=(0.99*LCO));
Logic12=(yb<(66.99*LCO));
Logic13=(yb>((66.99+42)*LCO));
Logic14=(xb<(30*LCO));
Logic15=(xb>(70*LCO));
Logic19 = (abs(zb - 1.0*LCO) < tol) & cruzaZ;

Mask9=(Logic11 & Logic12);
Mask10=(Logic11 & Logic13);
Mask11=(Logic11 & Logic14);
Mask12=(Logic11 & Logic15);
Mask22=(Logic19 & Logic12);
Mask23=(Logic19 & Logic13);
Mask24=(Logic19 & Logic14);
Mask25=(Logic19 & Logic15);

MaskReducao = MaskReducao2 | ...
            Mask9 | Mask10 | Mask11 | Mask12| ...
            Mask22 | Mask23 | Mask24 | Mask25;


EP(MaskReducao) = 1.E-5 * EP(MaskReducao);

disp(size(zb))
disp(size(zi))
disp(size(zf))
disp(size(Logic17))

% =====================================================
% FIBRAS DE REFORCO - TRES CAMADAS
% =====================================================

Ef = 7.24E10;

tolRel = 1e-8;

x = xb/LCO;
y = yb/LCO;
z = zb/LCO;

% Tres planos das fibras
MaskPlano05 = abs(z - 0.5) < tolRel;
MaskPlano15 = abs(z - 1.5) < tolRel;
MaskPlano25 = abs(z - 2.5) < tolRel;

MaskPlanosFibra = MaskPlano05 | MaskPlano15 | MaskPlano25;

% Coordenadas inteiras
xInteiro = abs(x - round(x)) < tolRel;
yInteiro = abs(y - round(y)) < tolRel;

% Coordenadas semi-inteiras
xMeio = abs(x - (floor(x) + 0.5)) < tolRel;
yMeio = abs(y - (floor(y) + 0.5)) < tolRel;

% Fibras nas direcoes X e Y
MaskFibraX = MaskPlanosFibra & xMeio & yInteiro;
MaskFibraY = MaskPlanosFibra & xInteiro & yMeio;

% Excluir regioes eliminadas
Mask26 = (MaskFibraX | MaskFibraY) & ~MaskReducao;

% Propriedades elasticas das fibras
% O modulo do epoxi ja foi atribuido anteriormente
E1(Mask26) = (Ef/EMAR)*E1(Mask26);
E2(Mask26) = (Ef/EMAR)*E2(Mask26);

% =====================================================
% INTRODUCAO DA FISSURA NO ACO
% =====================================================

% Comprimento da fissura
a = 39*LCO;

% Inclinacao em graus
theta = 45;

% Coordenadas do centro da fissura
xg = 50*LCO;
yg = 88*LCO;

% Extremidades da fissura
xmin = xg - (a/2)*cosd(theta);
xmax = xg + (a/2)*cosd(theta);

ymin = yg - (a/2)*sind(theta);
ymax = yg + (a/2)*sind(theta);

% Somente barras pertencentes ao aco
Logic20 = (zb >= 3*LCO);

% Barras entre xmin e xmax
Logic21 = (xb > xmin) & (xb < xmax);

% Equacao da reta central da fissura
yReta = ymin + (xb - xmin)*tand(theta);

% Faixa de 2*LCO na direcao Y
Logic22 = (yb > yReta - 1*LCO) & ...
          (yb < yReta + 1*LCO);

% Selecao final
MaskFissura = Logic20 & Logic21 & Logic22;

% Reducao da rigidez
EP(MaskFissura) = 1.E-5 * EP(MaskFissura);

% Conferencia
fprintf('Barras da fissura: %d\n', nnz(MaskFissura));

% A massa da matriz, do aco e dos vazios ja foi corrigida acima,
% imediatamente apos NODMASS(). Nao repetir correcao nodal aqui.

% !
%       IF (OU2.GT.0)&
%      &  OPEN (UNIT=12,ACCESS='SEQUENTIAL',FILE='LET.CON',&
%      &        STATUS='UNKNOWN')
% !
%       IF (OU3.GT.0)&
%      &  OPEN (UNIT=13,ACCESS='SEQUENTIAL',FILE='LET.AN1',&
%      &        STATUS='UNKNOWN')
% 
%       IF (OU5.GT.0)&
%      &  OPEN (UNIT=15,ACCESS='SEQUENTIAL',FILE='LET.RES',&
%      &        STATUS='UNKNOWN')
% !
%       IF (OU6.GT.0)&
%      &  OPEN (UNIT=16,ACCESS='SEQUENTIAL',FILE='LET.DIS',&
%      &        STATUS='UNKNOWN')
% !
%       IF (OU7.GT.0)&
%      &  OPEN (UNIT=17,ACCESS='SEQUENTIAL',FILE='LET.POS',&
%      &        STATUS='UNKNOWN')
% !
%       IF (OU8.GT.0)&
%      &  OPEN (UNIT=18,ACCESS='SEQUENTIAL',FILE='LET.ENE',&
%      &        STATUS='UNKNOWN')
% 
%  !                                     
%       IF (OU9.GT.0)&
%      &  OPEN (UNIT=19,ACCESS='SEQUENTIAL',FILE='Lei_Barra.TXT',&
%      &        STATUS='UNKNOWN')
% 
%    !                           
%       IF (OU9.GT.0)&
%      &  OPEN (UNIT=20,ACCESS='SEQUENTIAL',FILE='Energia_Barra.TXT',&
%      &        STATUS='UNKNOWN')
%    !--------------------------------------------------------------------------
%       IF (OU2.GT.0) THEN
% !     WRICON -----> Gravação de informação de entrada no arquivo CONNAME
%       CALL WRICON
%    ! We close(2) in the end of the Main program 
%       END IF
% !--------------------------------------------------------------------------
% !  LEITURA NOMES DE ARQUIVOS DO ANSYS
% !--------------------------------------------------------------------------
% 
%       !
%       IF (OU3.GT.0) THEN
%         DO 30 III=1,(NSTEP/OU3)+1
%           READ (13,*) ARMANS(III)
%    30   CONTINUE
%         CLOSE (13)
% !
% 
% !--------------------------------------------------------------------------
% !  GRAVAÇÃO DO PRIMEIRO ARQUIVO DO ANSYS
% !--------------------------------------------------------------------------
%         I3=1
%         OPEN (UNIT=13,ACCESS='SEQUENTIAL',FILE=ARMANS(I3),&
%      &        STATUS='UNKNOWN')
%         WRITE (13,*) '!ARQUIVO NO TEMPO =',TIME
%         CALL ANSYS (13,CONT1,CONT2,CONT3,CONT4,CONT5,CONT6)
%         CLOSE (13)
%       END IF


% -------------------------------------------------------------------------
if OU2 > 0
    % WRICON -----> writes input info to CONNAME
 ini=WRICON(fid_CON);   % stream passed in (or use global fid)
'ini=WRICON'
end

% -------------------------------------------------------------------------
% READ FILENAMES OF ANSYS FILES
% -------------------------------------------------------------------------
I3 = 1;
if OU3 > 0
ARMANS = cell(floor(NSTEP/OU3) + 1, 1);
for III = 1:floor(NSTEP/OU3) + 1
ARMANS{III} = strtrim(fgetl(fid_AN1));   % lê uma linha e remove espaços
    end
       fclose(fid_AN1);



    % --- Write first ANSYS file -----------------------------------------
    fid_AN1 = fopen(ARMANS{I3}, 'w');
    fprintf(fid_AN1, '!ARQUIVO NO TEMPO = %g\n', TIME);


      ini = ANSYS_Copia(fid_AN1, CONT1, CONT2, CONT3, CONT4, CONT5, CONT6);
     'ini=ANSYS'
   fclose(fid_AN1);
end

% 
% !--------------------------------------------------------------------------
% !  CONCLUIU-SE A MONTAGEM DO MODELO. INICIA-SE A INTEGRAÇÃO PASSO A PASSO
% !--------------------------------------------------------------------------
% !  LOOP PRINCIPAL DE INTEGRAÇÃO
% !--------------------------------------------------------------------------
% 
% 
% 
%       WRITE (*,*)
%       WRITE (*,*) '  *** PERFORMING THE INTEGRATION STEP BY STEP **  '
%       WRITE (*,*)
% !
%       STEPP=1
% !
%       DO 60 STEP = STEPP,NSTEP
% !
% !--------------------------------------------------------------------------
% !  Escreve uma informação de controle na tela em cada passo de integração
% !--------------------------------------------------------------------------
%       WRITE (*,40) STEP
% 40    FORMAT ('   STEP     ======================>  ',I7,' / ',I2)
% 
% !
% !--------------------------------------------------------------------------
% !  Controle de mudança de amortecimento
% !--------------------------------------------------------------------------
%       IF (STEP.EQ.NDAMF) THEN
%         DF1=DF2
%         DMP1=1.0D00-DF2*DT/2.0D00
%         DMP2=1.0D00+DF2*DT/2.0D00
%       END IF
% !
% !--------------------------------------------------------------------------
% !  Incremento do tempo
% !--------------------------------------------------------------------------
%       TIME=TIME+DT
%       WRITE (*,50) TIME
% 50    FORMAT ('   TIME     ======================>  ',F12.8)
%-------------------------------------------------------------------------
% !  Avaliação das forças e determinação da nova configuração
% !--------------------------------------------------------------------------
% !     CALOAD -----> Cálculo das cargas sobre os nós
%       CALL CALOAD
%       write(*,*) 'caload' 
% !
% !  Especifica-se a carga em cada tempo e calcula-se a energia externa
% !  fornecida ao sistema. Em função do exemplo tem-se que modificar esta
% !  subrotina.
% !
% !--------------------------------------------------------------------------
% !     CABAR ------> compute de nodal forces in the structure, in CABAR we have a loop over all the bars and
%       ! we call the constit subroutine to compute the forces in each bars
%       CALL CABAR
%  write(*,*) 'cabar' 
% !
% !--------------------------------------------------------------------------
% !     IGN012 -----> Cálculo das coordenadas no tempo TIME aplicando
% !                   diferenças finitas centrais na equação do movimento
%       CALL IGN012
% !  Nesta subrotina realiza-se a integração via método das diferenças
% !  finitas centrais. Tem-se um loop sobre os nós. Entra-se com as forças
% !  calculadas em CABAR e com a configuração do sistema nos dois tempos
% !  anteriores e determina-se a configuração no tempo novo.
% !

% =====================================================
% MASSA POR PLANO Z - INCLUINDO NOS CENTRAIS
% =====================================================

fprintf('\n===== MASSA POR PLANO Z =====\n');

zn = W0(1:NNT);
zn = zn(:);

mnos = MAS(1:NNT);
mnos = mnos(:);

tolZ = 1e-9;

for k = 0:12

    zPlano = 0.5*k*LCO;

    MaskPlano = abs(zn-zPlano) < tolZ;

    MassaPlano = sum(mnos(MaskPlano));

    fprintf('Z = %.1f mm | Nos = %d | Massa = %.8f kg\n', ...
        0.5*k, nnz(MaskPlano), MassaPlano);

end

fprintf('\nMassa total = %.8f kg\n', sum(mnos));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====================================================
% CONFERENCIA ANTES DA INTEGRACAO
% =====================================================

fprintf('\n===== VERIFICACAO INICIAL =====\n');

fprintf('Nos: %d\n', NNT);
fprintf('Barras: %d\n', NBT);

fprintf('Nos restringidos em X: %d\n', nnz(BDX));
fprintf('Nos restringidos em Y: %d\n', nnz(BDY));
fprintf('Nos restringidos em Z: %d\n', nnz(BDZ));

fprintf('Massa total: %.8e kg\n', sum(MAS(:)));

fprintf('DT: %.8e s\n', DT);
fprintf('TTOT: %.8e s\n', TTOT);
fprintf('NSTEP original: %d\n', NSTEP);

assert(all(isfinite(MAS(:))), ...
    'Existem massas nodais invalidas');

assert(all(isfinite(EP(1:NBT))), ...
    'Existem propriedades EP invalidas');

% Teste inicial com apenas 10 passos
% NSTEP = min(NSTEP,10);

fprintf('NSTEP para teste: %d\n', NSTEP);
return
% -------------------------------------------------------------------------
% MODEL ASSEMBLY COMPLETE. START STEP-BY-STEP INTEGRATION
% -------------------------------------------------------------------------
fprintf('\n  *** PERFORMING THE INTEGRATION STEP BY STEP **  \n\n');

STEPP = 1;

for STEP = STEPP:NSTEP

    fprintf('   STEP     ======================>  %7d / %2d\n', STEP, 0);

    % ---------------------------------------------------------------------
    % Damping change control
    % ---------------------------------------------------------------------
    if STEP == NDAMF
        DF1  = DF2;
        DMP1 = 1.0 - DF2 * DT / 2.0;
        DMP2 = 1.0 + DF2 * DT / 2.0;
    end

    % ---------------------------------------------------------------------
    % Time increment
    % ---------------------------------------------------------------------
    TIME = TIME + DT;
    fprintf('   TIME     ======================>  %12.8f\n', TIME);
    fprintf('%25.12e   %25.12e   %25.12e   %25.12e\n', ...
                TIME, REACAOY, ACZ(NCON(1)), V(2,NCON(2)));
        
    fprintf('\n');
    % ---------------------------------------------------------------------
    % Evaluate forces and determine new configuration
    % ---------------------------------------------------------------------
    ini=CALOAD();
   'ini=CALOAD()'
    
   ini=CABAR(fid_Lei);
   'ini=CABAR()'

   ini=IGN012();     % central finite differences integration
   'ini=IGN012()'
% !--------------------------------------------------------------------------
% !  Geração dos arquivos de saída
% !--------------------------------------------------------------------------
% !  SAÍDA - VISUALIZAÇÃO NO ANSYS
% !--------------------------------------------------------------------------
%      IF (PR3.EQ.OU3) THEN
%         I3=I3+1
%         PR3=1
%         OPEN (UNIT=11,ACCESS='SEQUENTIAL',FILE=ARMANS(I3), &
%       &                STATUS='UNKNOWN')
%         WRITE (11,*) '!ARQUIVO NO TEMPO =',TIME
%         CALL ANSYS (11,CONT1,CONT2,CONT3,CONT4,CONT5,CONT6)
%         CLOSE (11)
%       ELSE
%         PR3=PR3+1
%       END IF
% !
% !--------------------------------------------------------------------------
% !  SAÍDA - REAÇÕES NODAIS
% !--------------------------------------------------------------------------
%       ! write(*,*) 'antes de calres'
% !      CALRES -----> Cálculo de resultados que serão gravados em RESNAME
%        CALL CALRES
%       !  write(*,*) 'calres' 
%       !
% !--------------------------------------------------------------------------
% !  SAÍDA - CONTROLE DE DESLOCAMENTOS, VELOCIDADES E ACELERAÇÕES
% !--------------------------------------------------------------------------
% !
%       IF (PR67.EQ.OU6) THEN
%         CALL WRIDIPO
% !!  Esta subrotina grava os arquivos *.dis e *.pos que tem os valores dos
% !  deslocamentos, velocidades e acelerações nos nós de controle.
% !
%         PR67=1
%       ELSE
%         PR67=PR67+1
%       END IF
% !--------------------------------------------------------------------------
% !  SAÍDA - AVALIAÇÃO DA ENERGIA
% !--------------------------------------------------------------------------
% !
% !     ENERG2 -----> Cálculo das energias que serão gravadas em ENENAME
%       CALL ENERG2
% !  Esta subrotina calcula em cada tempo as energias que estão se
% !  transformando no processo.
% !
%       IF (PR8.EQ.OU8) THEN
%         CALL WRIENE
% !  Esta subrotina grava o arquivo *.ene que tem os valores das energias.
% !
%         PR8=1
%       ELSE
%         PR8=PR8+1
%       END IF
% 
% !
% !--------------------------------------------------------------------------
%    60 CONTINUE
% !--------------------------------------------------------------------------
% !
%       IF (RUN.EQ.NRUN) GO TO 70
% !
%       IF (OU5.GT.0) CLOSE (5)
%       IF (OU6.GT.0) CLOSE (6)
%       IF (OU7.GT.0) CLOSE (7)
%       IF (OU8.GT.0) CLOSE (8)
% !    
% !
%      70 CONTINUE
% !
% 
%       WRITE (*,*)
%       WRITE (*,*) '  ****************** THE END *******************  '
%       WRITE (*,*)
% 
%      close(22)
%      close(2)
%     end program Console1
 % ---------------------------------------------------------------------
    %  Generate output files
    % ---------------------------------------------------------------------

    % Output - ANSYS visualisation
    if PR3 == OU3
        I3 = I3 + 1;
        PR3 = 1;
       fid_AN1 = fopen(ARMANS{I3}, 'w');
               fprintf(fid_AN1, '!ARQUIVO NO TEMPO = %g\n', TIME);
       
         ini=ANSYS(fid_AN1, CONT1, CONT2, CONT3, CONT4, CONT5, CONT6);
        'ini=ANSYS()'
        fclose(fid_AN1);
    else
        PR3 = PR3 + 1;
    end

  
    % ---------------------------------------------------------------------
    %  Output - control of displacements, velocities and accelerations
    % ---------------------------------------------------------------------
    if PR67 == OU6
    % Output - nodal reactions
    % CALRES -----> Compute results to be written to RESNAME
      
        ini=WRIDIPO(fid_DIS, fid_POS);
        'ini=WRIDIPO'
        PR67 = 1;
    else
        PR67 = PR67 + 1;
    end

 if PR9 == OU9
    % Output - nodal reactions
    % CALRES -----> Compute results to be written to Sim01.dat
       ini=CALRES(fid_RES);
       'ini=CALRES()'

        PR9 = 1;
    else
        PR9 = PR9 + 1;
    end



    % ---------------------------------------------------------------------
    %  Output - energy evaluation
    % ---------------------------------------------------------------------
    % ENERG2 -----> Compute energies to be written to ENENAME
    ini=ENERG2(fid_EnB);
     'ini=ENERG2()'
    if PR8 == OU8
        ini=WRIENE(fid_ENE);
        'ini=WRIENE'

        PR8 = 1;
    else
        PR8 = PR8 + 1;
    end

end
% -------------------------------------------------------------------------

if RUN == NRUN
    % GO TO 70
else
    if OU5 > 0 && fid_RES >= 0, fclose(fid_RES); end
    if OU6 > 0 && fid_DIS>= 0, fclose(fid_DIS); end
    if OU7 > 0 && fid_POS >= 0, fclose(fid_POS); end
    if OU8 > 0 && fid_ENE >= 0, fclose(fid_ENE); end
if OU9 > 0, fid_Lei = fopen('Lei_Barra.TXT',   'w'); end
if OU9 > 0, fid_EnB = fopen('Energia_Barra.TXT','w'); end
if OU9 > 0, fid_RES = fopen('SIM01.DAT','w'); end
end

% 70 CONTINUE

fprintf('\n');
fprintf('  ****************** THE END *******************  \n');
fprintf('\n');

% Close remaining open files
%if fid19 >= 0, fclose(fid19); end
%if fid20 >= 0, fclose(fid20); end
%if fid12 >= 0, fclose(fid12); end

% End of program Console1




