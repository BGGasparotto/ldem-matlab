function ini = ANSYS_Copia(IFILE,CONT1,CONT2,CONT3,CONT4,CONT5,CONT6)
global  INT INB U V W U0 V0 W0 BDX BDY BDZ VLX VLY VLZ 
global ACX ACY ACZ FRX FRY FRZ FX FY FZ MAS MSS CN ENR
global EDG E1 E2 LI DEQB F YNGB GFRB LCRB EPRB EP SS EM
global EN KR KR1 DSS DM FS CAF DT GFR DEQ LCR TEST EST 
global DF DF1 DF2 DMP1 DMP2 TDAF CON M N L M1 N1 L1 MNL
global NNO NNS NNT NBA NBS NBT INCO NCON STEP NSTEP NDAMF
global RUN NRUN LCO TIME TTOT PSN ROH YNG EPR EPMAX EPMIN 
global MEAN  COEV ENEX ENCN ENEL ENGD ENDP ENIN DEN DERIV 
global XB YB ZB UBAR VBAR WBAR REACAOX REACAOY REACAOZ 
global STRAIN STRESS  DM DS

%==========================================================================
% HEADER / SETUP
%==========================================================================
fprintf(IFILE, '/PREP7\n');
fprintf(IFILE, '/TITLE, TEMPO =%g\n', TIME);
% LINK4: TYPE=1 temporario; 2=aco; 3=epoxi; 4=fibras; 5=regioes enfraquecidas.
for tipo = 1:5
    fprintf(IFILE, 'ET,%d,4\n', tipo);
end
fprintf(IFILE, 'MP,EX,1,2E10\n');
fprintf(IFILE, 'MP,PRXY,1,0.25\n');
fprintf(IFILE, 'MP,DENS,1,2500\n');
fprintf(IFILE, 'R,1,0.13,0.0014,0.36\n');
fprintf(IFILE, '\n');

%--------------------------------------------------------------------------
% NODAL COORDINATES
%--------------------------------------------------------------------------
fprintf(IFILE, '!COORDENADAS NODAIS\n');

for I = 1:NNT
    fprintf(IFILE, 'N,%d,%g,%g,%g\n', ...
            I, single(U(3,I)), single(V(3,I)), single(W(3,I)));
end
fprintf(IFILE, '\n');

%--------------------------------------------------------------------------
% ELEMENT CONNECTIVITIES
%--------------------------------------------------------------------------
fprintf(IFILE, '!CONECTIVIDADES DOS ELEMENTOS\n');

CONT1 = 0;
CONT2 = 0;
CONT3 = 0;
CONT4 = 0;
CONT5 = 0;
CONT6 = 0;
CONT7 = 0;

% Na criacao todos os elementos sao temporariamente TYPE=1.
% O TYPE definitivo sera atribuido pelas mascaras geometricas abaixo.
% As contagens DM/FS originais continuam independentes dos tipos geometricos.
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

%==========================================================================
% REGIOES GEOMETRICAS - MESMAS MASCARAS DO MAIN1.M
%==========================================================================
% Trabalhar somente com as NBT barras existentes e vetores coluna.
xb = XB(1:NBT); xb = xb(:);
yb = YB(1:NBT); yb = yb(:);
zb = ZB(1:NBT); zb = zb(:);

% Nos inicial e final de cada barra, na ordem do vetor CN.
NOi = CN(1:2:2*NBT); NOi = NOi(:);
NOf = CN(2:2:2*NBT); NOf = NOf(:);

NOi = NOi(:);
NOf = NOf(:);

% Coordenadas Z dos nos
zi = W0(NOi);
zf = W0(NOf);

zi = zi(:);
zf = zf(:);

% Uma barra que cruza um plano Z possui nos com cotas diferentes.

tol = 1e-9;  % tolerancia em metros para comparar centros das barras
cruzaZ = abs(zi - zf) > tol;

% Degrau superior do reparo: Z entre 2 e 3 mm
Logic1  = (zb <= 2.99*LCO) & (zb >= 1.99*LCO);
Logic2  = (yb < 46.99*LCO);
Logic3  = (yb > (47.01+82)*LCO);
Logic4  = (xb < 10*LCO);
Logic5  = (xb > 90*LCO);
Logic17 = (abs(zb - 3.0*LCO) < tol) & cruzaZ;


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

% Degrau intermediario: Z entre 1 e 2 mm
Logic6  = (zb < 1.99*LCO) & (zb > 0.99*LCO);
Logic7  = (yb < 56.99*LCO);
Logic8  = (yb > (56.99+62)*LCO);
Logic9  = (xb < 20*LCO);
Logic10 = (xb > 80*LCO);
Logic18 = (abs(zb - 2.0*LCO) < tol) & cruzaZ;


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

% Degrau inferior: Z entre 0 e 1 mm
Logic11 = (zb <= 0.99*LCO);
Logic12 = (yb < 66.99*LCO);
Logic13 = (yb > (66.99+42)*LCO);
Logic14 = (xb < 30*LCO);
Logic15 = (xb > 70*LCO);
Logic19 = (abs(zb - 1.0*LCO) < tol) & cruzaZ;

Mask9  = Logic11 & Logic12;
Mask10 = Logic11 & Logic13;
Mask11 = Logic11 & Logic14;
Mask12 = Logic11 & Logic15;
Mask22 = Logic19 & Logic12;
Mask23 = Logic19 & Logic13;
Mask24 = Logic19 & Logic14;
Mask25 = Logic19 & Logic15;

% Exatamente a mesma uniao de barras enfraquecidas do main1.m.
MaskReducao = MaskReducao2 | Mask9 | Mask10 | Mask11 | Mask12|...
                Mask22 | Mask23 | Mask24 | Mask25;

% ========================================================================
% FIBRAS NAS TRES CAMADAS DO REPARO
% Centros das fibras: Z/LCO = 0.5, 1.5 e 2.5
% Direcao X: centro x semi-inteiro e y inteiro
% Direcao Y: centro x inteiro e y semi-inteiro
% Os centros (x,y) inteiros em ambas as coordenadas nao sao fibras.
% ========================================================================
tolRel = 1e-8;
x = xb/LCO;
y = yb/LCO;
z = zb/LCO;

% Selecionar os centros nas tres cotas das fibras.
MaskPlano05 = abs(z - 0.5) < tolRel;
MaskPlano15 = abs(z - 1.5) < tolRel;
MaskPlano25 = abs(z - 2.5) < tolRel;
MaskPlanosFibra = MaskPlano05 | MaskPlano15 | MaskPlano25;

% Geometria dos centros nas direcoes X e Y.
xInteiro = abs(x - round(x)) < tolRel;
yInteiro = abs(y - round(y)) < tolRel;
xMeio = abs(x - (floor(x) + 0.5)) < tolRel;
yMeio = abs(y - (floor(y) + 0.5)) < tolRel;

% Mesma selecao de fibra para cada plano.
MaskFibraX = MaskPlanosFibra & xMeio & yInteiro;
MaskFibraY = MaskPlanosFibra & xInteiro & yMeio;

% Regioes enfraquecidas prevalecem sobre as fibras.
MaskEliminadas = MaskReducao;
MaskFibras = (MaskFibraX | MaskFibraY) & ~MaskEliminadas;
MaskAco = (zb >= 3*LCO) & ~MaskEliminadas & ~MaskFibras;
MaskMatriz = (zb < 3*LCO) & ~MaskEliminadas & ~MaskFibras;

% Verificar que cada barra esta classificada em exatamente um grupo.
assert(numel(MaskEliminadas)==NBT && numel(MaskFibras)==NBT && ...
       numel(MaskAco)==NBT && numel(MaskMatriz)==NBT, ...
       'Dimensoes das mascaras incompatíveis com NBT');
assert(all((double(MaskEliminadas) + double(MaskFibras) + ...
            double(MaskAco) + double(MaskMatriz)) == 1), ...
       'Barras sem classificacao ou classificadas em mais de uma regiao');

fprintf('Aco (TYPE=2): %d barras\n', nnz(MaskAco));
fprintf('Matriz (TYPE=3): %d barras\n', nnz(MaskMatriz));
fprintf('Fibras X (TYPE=4): %d barras\n', nnz(MaskFibraX & MaskFibras));
fprintf('Fibras Y (TYPE=4): %d barras\n', nnz(MaskFibraY & MaskFibras));
fprintf('Fibras em Z=0.5*LCO: %d barras\n', nnz(MaskPlano05 & MaskFibras));
fprintf('Fibras em Z=1.5*LCO: %d barras\n', nnz(MaskPlano15 & MaskFibras));
fprintf('Fibras em Z=2.5*LCO: %d barras\n', nnz(MaskPlano25 & MaskFibras));
fprintf('Enfraquecidas (TYPE=5): %d barras\n', nnz(MaskEliminadas));

% ==========================================================
% TIPOS DOS ELEMENTOS NO ANSYS
% ==========================================================
TIPOBARRA = zeros(NBT,1);
TIPOBARRA(MaskAco)        = 2;  % Amarelo
TIPOBARRA(MaskMatriz)     = 3;  % Verde
TIPOBARRA(MaskFibras)     = 4;  % Vermelho
TIPOBARRA(MaskEliminadas) = 5;  % Cinza
assert(all(TIPOBARRA >= 2 & TIPOBARRA <= 5), 'TIPOBARRA contem valores invalidos');

% Elementos sao numerados na ordem de criacao (I = 1:NBT).
for I = 1:NBT
    fprintf(IFILE, 'EMODIF,%d,TYPE,%d\n', I, TIPOBARRA(I));
end
fprintf(IFILE, 'ALLSEL,ALL\n');

% ==========================================================
% COLORIR POR TIPO (nao por selecao momentanea de elementos)
% /PNUM,TYPE,1 ativa cores por TYPE; /NUMBER,1 oculta numeros.
% /COLOR,NUM define a cor do numero de tipo em cada mapa discreto.
% ==========================================================
fprintf(IFILE, '/COLOR,NUM,YELL,2\n'); % Aco
fprintf(IFILE, '/COLOR,NUM,GREE,3\n'); % Matriz
fprintf(IFILE, '/COLOR,NUM,RED,4\n');  % Fibras
fprintf(IFILE, '/COLOR,NUM,LGRA,5\n'); % Enfraquecidas
fprintf(IFILE, '/PNUM,TYPE,1\n');
fprintf(IFILE, '/NUMBER,1\n');
fprintf(IFILE, 'ALLSEL,ALL\n');

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%----------------------------------------------------------
%--------------------------------------------------------------------------
% BOUNDARY CONDITIONS
%--------------------------------------------------------------------------
fprintf(IFILE, '!CONDIÇÕES DE CONTORNO\n');

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

fprintf(IFILE, '\n');
% --- Commented-out D-constraints from original code kept here for reference
% fprintf(IFILE, 'D,%d,UY\n', ncon(1));
% ... (etc.)

%--------------------------------------------------------------------------
% APPLIED LOADS
%--------------------------------------------------------------------------
fprintf(IFILE, '!CARGAS APLICADAS\n');
 if TEST == 102
      fprintf(IFILE, 'F,%d,FY,%d\n', 17664, 100);
           end
% %     %VF2 = 400.0e-5;   % m/s
%       VF2 = 400.0e-4 ;
%      for I = 1:NNT
%          LOGIC1 = (W0(I) < (0.51) * LCO) && ...
%                  (W0(I) > (0.49) * LCO);
%          LOGIC2 = (W0(I) < (L1 - 0.49)  * LCO) && ...
%                   (W0(I) > (L1 - 0.51)  * LCO);
%          LOGIC3 = (V0(I) < (N1-0.49) * LCO) && ...   
%                   (V0(I) > (N1-0.51) * LCO);
%           LOGIC4 = (V0(I) < ((N1-1)-0.49) * LCO) && ...   
%                   (V0(I) > ((N1-1)-0.51) * LCO);
%          LOGIC5 = (V0(I) < ((N1-2)-0.49) * LCO) && ...   
%                   (V0(I) > ((N1-2)-0.51) * LCO);
%           LOGIC6 = (V0(I) < ((N1-3)-0.49) * LCO) && ...   
%                   (V0(I) > ((N1-3)-0.51) * LCO);
%  
%  if LOGIC1 && LOGIC3
%               fprintf(IFILE, 'F,%d,FZ,%d\n', I, -100);
%           end
%     
%  if LOGIC2 && LOGIC3
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, 100);
%           end
%     
%  if LOGIC1 && LOGIC4
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, -50);
%           end
%     
%  if LOGIC2 && LOGIC4
%              fprintf(IFILE, 'F,%d,FZ,%d\n', I, 50);
%  end
% if LOGIC1 && LOGIC5
%              fprintf(IFILE, 'F,%d,FZ,%d\n', I, 50);
%           end
%     
%  if LOGIC2 && LOGIC5
%              fprintf(IFILE, 'F,%d,FZ,%d\n', I,-50);
%  end
% if LOGIC1 && LOGIC6
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, 100);
%           end
%     
%  if LOGIC2 && LOGIC6
%              fprintf(IFILE, 'F,%d,FZ,%d\n', I, -100);
% 
%  end 
%  end 
% end       
% if TEST == 102
%     for I = 1:NNT
%         % Filtra as extremidades (Esquerda e Direita)
%         LOGIC1 = (W0(I) < 0.51 * LCO) && (W0(I) > 0.49 * LCO);
%         LOGIC2 = (W0(I) < (L1 - 0.49) * LCO) && (W0(I) > (L1 - 0.51) * LCO);
%         
%         % Filtra as fibras superior e inferior do conjunto tracionado
%         LOGIC3 = (V0(I) < (N1 - 0.49) * LCO) && (V0(I) > (N1 - 0.51) * LCO);
%         LOGIC6 = (V0(I) < ((N1 - 3) - 0.49) * LCO) && (V0(I) > ((N1 - 3) - 0.51) * LCO);
% 
%         % Imprime as setas para visualização do Momento Fletor nas pontas
%         if LOGIC1 && LOGIC3
%             % Ponta esquerda, fibra superior (puxando para baixo)
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, -100); 
%         elseif LOGIC2 && LOGIC3
%             % Ponta direita, fibra superior (puxando para cima)
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, 100);  
%         elseif LOGIC1 && LOGIC6
%             % Ponta esquerda, fibra interna (empurrando para cima)
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, 100);  
%         elseif LOGIC2 && LOGIC6
%             % Ponta direita, fibra interna (empurrando para baixo)
%             fprintf(IFILE, 'F,%d,FZ,%d\n', I, -100); 
%         end
%     end
% end

fprintf(IFILE, '\n');

%--------------------------------------------------------------------------
% VIEW / PLOT
%--------------------------------------------------------------------------
fprintf(IFILE, '/VIEW,1,0,-1,0\n');
fprintf(IFILE, '/PBC,ALL,1\n');
fprintf(IFILE, 'EPLOT\n');
fprintf(IFILE, '\n');

%--------------------------------------------------------------------------
% BAR STATE SUMMARY
%--------------------------------------------------------------------------
fprintf(IFILE, '!ESTADO DAS BARRAS\n');
fprintf(IFILE, '!BARRAS INTACTAS =    %d\n',    CONT1 + CONT3 + CONT6);
fprintf(IFILE, '!BARRAS DANIFICADAS = %d\n',   CONT2 + CONT4);
fprintf(IFILE, '!BARRAS ROMPIDAS =    %d\n',    CONT5);
fprintf(IFILE, '\n');

%--------------------------------------------------------------------------
% Separate bar-state summary file
% Fortran: OPEN (UNIT=30,FILE='BARR_ROMP.dat',STATUS='UNKNOWN')
% FORMAT (E12.6,3X,E12.6,3X,I7,3X,I7,3X,I7)
%--------------------------------------------------------------------------
VF=1.0
fid30 = fopen('BARR_ROMP.dat', 'w');
fprintf(fid30, '%12.6E   %12.6E   %7d   %7d   %7d\n', ...
        TIME, VF, CONT1, CONT2, CONT5);
fclose(fid30);
 ini=1;
 end