%% homogenization.m
% Compute homogenized compliance tensor of lattice of Bernoulli's beam 
% see https://doi.org/10.1016/j.commatsci.2011.07.014
% for detailed explanations
% Dos Reis F.
% 23.01.2021
function [ME1,MS4]=homogenization(Tb)
    %%initiation of variables
    global L1 L2 Y1 Y2 nbeams nnodes nodes Ob Eb delta1 delta2 Elast
    ME1=zeros(1,10);
    Lb=zeros(nbeams,1);   % length of beams
    ka=zeros(nbeams,1);   % traction modulus
    kf=zeros(nbeams,1);   % bending modulus
    e=zeros(nbeams,2);   % matrix nbeams x 2 of director's beam
    eT=zeros(nbeams,2);   % matrix nbeams x 2 of perpendicular beam's vector
    % order components in depli : (dUx/dx dUx/dy dUy/dx dUy/dy)
    deplN=zeros(nbeams,4); % displacements imposed on top and right of the cell 
    deplT=zeros(nbeams,4);
    deplMo=zeros(nbeams,4);
    deplMe=zeros(nbeams,4);
    % forces
    Nvs=zeros(nbeams,nnodes*3); % axial force
    Tvs=zeros(nbeams,nnodes*3); % shearing force
    N=zeros(nbeams,4);
    T=zeros(nbeams,4);
    Mvos=zeros(nbeams,nnodes*3); % moment force origin node
    Mves=zeros(nbeams,nnodes*3); % moment force end node
    K=zeros(nnodes*3,nnodes*3);
    B=zeros(nnodes*3,4);
    X=zeros(nnodes*3,4);
    omega_rigid=zeros(nnodes*3,4);
    S1=zeros(2,4);
    S2=zeros(2,4);
    Mat=zeros(4,4);
    MS=zeros(4,4);
    MatVo=zeros(3,3);
    MS2=zeros(3,3);
    MS3=zeros(3,3);
    MS4=zeros(1,6);
    %% DU Generation
    P1=[L1*Y1(1) L1*Y1(2)]; 
    P2=[L2*Y2(1) L2*Y2(2)];
    dU1=[P1(1) P1(2) 0 0 ; 0 0 P1(1) P1(2)];
    dU2=[P2(1) P2(2) 0 0; 0 0 P2(1) P2(2)];
    dU=[dU1;dU2]; %réorganiser cela sous la forme d'une seule matrice ?

    %% Building forces for each beam
    for i=1:nbeams 
        % initial parameters
        indexn1=Ob(i);  
        indexn2=Eb(i);
        Posn1=nodes(indexn1,1:2);
        Posn2=nodes(indexn2,1:2)+delta1(i)*Y1*L1+delta2(i)*Y2*L2;
        ObEb=Posn2-Posn1; % ObEb vector
        Lb(i)=sqrt(ObEb(1,1)^2+ObEb(1,2)^2); % length of the beam
        e(i,1:2)=ObEb/Lb(i);
        ka(i)=Elast(i)*Tb(i)/Lb(i);
        kf(i)=Elast(i)*(Tb(i)/Lb(i))^3;
        % N(beam)
        kae=ka(i)*e(i,1:2);
        deplN(i,1:4)=delta1(i)*kae*dU1+delta2(i)*kae*dU2; 
        indexn1K=3*indexn1-2;    % position index node 1 in global stiffness matrix K
        indexn2K=3*indexn2-2;
        % (il faut que j'utilise un appel vectorisé des matrices)
        Nvs(i,indexn2K:indexn2K+1)=kae;
        Nvs(i,indexn1K:indexn1K+1)=Nvs(i,indexn1K:indexn1K+1)-kae;
        % T(beam)
        eT(i,1)=-e(i,2);    % calculation of perpendicular beam's vector
        eT(i,2)=e(i,1);
        kfe=kf(i)*eT(i,1:2);    % modulus of displacements in shear forces
        deplT(i,1:4)=delta1(i)*kfe*dU1+delta2(i)*kfe*dU2; 
        % (il faut que j'utilise un appel vectorisé des matrices)
        Tvs(i,indexn2K:indexn2K+1)=kfe;
        Tvs(i,indexn1K:indexn1K+1)=Tvs(i,indexn1K:indexn1K+1)-kfe;
        kbL=-kf(i)*Lb(i)/2;   % modulus of rotations in shear force (mettre l'equation)
        Tvs(i,indexn2K+2)=kbL;
        Tvs(i,indexn1K+2)=Tvs(i,indexn1K+2)+kbL;
        %Mo & Me
        kfe2=-kf(i)*eT(i,1:2)*Lb(i)/2;
        kbL2=kf(i)*Lb(i)^2/6;
        kbL3=-Lb(i)/2;
        Mves(i,indexn2K:indexn2K+1)=kfe2;
        Mves(i,indexn1K:indexn1K+1)=Mves(i,indexn1K:indexn1K+1)-kfe2;
        Mves(i,indexn1K+2)=kbL2;
        Mves(i,indexn2K+2)=Mves(i,indexn2K+2)+2*kbL2;
        deplMe(i,1:4)=kbL3*deplT(i,:);

        Mvos(i,indexn2K:indexn2K+1)=kfe2;
        Mvos(i,indexn1K:indexn1K+1)=Mvos(i,indexn1K:indexn1K+1)-kfe2;
        Mvos(i,indexn1K+2)=2*kbL2;
        Mvos(i,indexn2K+2)=Mvos(i,indexn2K+2)+kbL2;
        deplMo(i,1:4)=kbL3*deplT(i,:);
        % Build of stiffness matrix [K] & displacement [B]
        equ1=reshape(e(i,:),[2,1])*Nvs(i,:)+reshape(eT(i,:),[2,1])*Tvs(i,:);
        K(indexn1K:indexn1K+1,:)=K(indexn1K:indexn1K+1,:)-equ1;
        K(indexn2K:indexn2K+1,:)=K(indexn2K:indexn2K+1,:)+equ1;
        K(indexn1K+2,:)=K(indexn1K+2,:)+Mvos(i,:);
        K(indexn2K+2,:)=K(indexn2K+2,:)+Mves(i,:);
        equ2=reshape(e(i,:),[2,1])*deplN(i,:)+reshape(eT(i,:),[2,1])*deplT(i,:);
        B(indexn1K:indexn1K+1,:)=B(indexn1K:indexn1K+1,:)+equ2;
        B(indexn2K:indexn2K+1,:)=B(indexn2K:indexn2K+1,:)-equ2;    
        B(indexn1K+2,:)=B(indexn1K+2,:)-deplMo(i,:);
        B(indexn2K+2,:)=B(indexn2K+2,:)-deplMe(i,:);
    end;
    %% Solve [X] as [K][X]=[B]
    % [K] is singular
    % we assume dx1=0 dy1=0
    X(3:3*nnodes,:)=K(3:3*nnodes,3:3*nnodes)\B(3:3*nnodes,:);
    for i=3:3:3*nnodes    % il faudrait le faire qu'une seule fois
        omega_rigid(i,2:3)=[1/2 -1/2];
    end;
    X=X+omega_rigid;
    X(1,:)=[0 0 0 0];
    X(2,:)=[0 0 0 0];
    %% Evaluation of forces
    N=Nvs*X+deplN;
    T=Tvs*X+deplT;
    %% Evaluation of stiffness vectors sum Si
    for i=1:nbeams % il faut vectoriser cette partie 
        S1=S1+delta1(i)*(reshape(e(i,:),[2,1])*N(i,:)+reshape(eT(i,:),[2,1])*T(i,:));
        S2=S2+delta2(i)*(reshape(e(i,:),[2,1])*N(i,:)+reshape(eT(i,:),[2,1])*T(i,:));
    end;
    %% Evaluation of stiffness tensor
    g=[P1;P2];
    detg=det(g);
    Mat=reshape(1/detg*(reshape(P1,[2,1])*reshape(S1,[1,8])+...
        reshape(P2,[2,1])*reshape(S2,[1,8])),[4,4]);
    % rearrange (dUx/dx dUx/dy dUy/dx dUy/dy) -> (dUx/dx, dUy/dy, dUx/dy, dUy/dx)
    % to obtain Matrix in natural basis notation
    L=Mat(4,:);
    Mat(4,:)=Mat(3,:);
    Mat(3,:)=Mat(2,:);
    Mat(2,:)=L;
    C=Mat(:,4);
    Mat(:,4)=Mat(:,3);
    Mat(:,3)=Mat(:,2);
    Mat(:,2)=C;
    % Matrix in Voigt notation
    MatVo(1:3,1:2)=Mat(1:3,1:2);
    MatVo(1,3)=(Mat(1,3)+Mat(1,4))/2;
    MatVo(2,3)=(Mat(2,3)+Mat(2,4))/2;
    MatVo(3,3)=(Mat(3,3)+Mat(3,4))/2;
    %% evaluation of compliance matrix
    % MS in natural form
    MS=inv(Mat);
    % in Voigt form
    MS2(1:2,1:2)=MS(1:2,1:2);
    MS2(1,3)=MS(1,3)+MS(1,4);
    MS2(2,3)=MS(2,3)+MS(2,4);
    MS2(3,1)=MS(3,1)+MS(4,1);
    MS2(3,2)=MS(3,2)+MS(4,2);
    MS2(3,3)=MS(3,3)+MS(3,4)+MS(4,3)+MS(4,4);
    % in basis tensor notation
    MS3=MS2;
    MS3(1,3)=MS3(1,3)/sqrt(2);
    MS3(2,3)=MS3(2,3)/sqrt(2);
    MS3(3,1)=MS3(3,1)/sqrt(2);
    MS3(3,2)=MS3(3,2)/sqrt(2);
    MS3(3,3)=MS3(3,3)/2;
    %% evaluation of mechanical moduli
    % mise sous forme d'un vecteur linéaire (1,6) plutot qu'un tenseur pour la
    % valeur de retour : [S11, S22, S66, S26, S16, S12]
    MS4(1,1)=MS3(1,1);MS4(1,2)=MS3(2,2);MS4(1,3)=MS3(3,3);
    MS4(1,4)=MS3(2,3);MS4(1,5)=MS3(1,3);MS4(1,6)=MS3(1,2);
    
    %K: = simplify(1 / (MS3[1, 1] + MS3[2, 2] + MS3[1, 2] + MS3[2, 1]));
    ME1(1,1)=1/(MS3(1,1)+MS3(2,2)+MS3(1,2)+MS3(2,1));

	%Ex: = simplify(1 / MS3[1, 1]);
    ME1(1,2)=1/MS3(1,1);

	%Ey: = simplify(1 / MS3[2, 2]);
    ME1(1,3)=1/MS3(2,2);

	%nuyx: = simplify(-MS3[1, 2] * Ey); 
    ME1(1,4)=-MS3(1,2)*ME1(1,3);

	%nuxy: = simplify(-MS3[2, 1] * Ex);
    ME1(1,5)=-MS3(2,1)*ME1(1,2);

	%muxy: = simplify(1 / (2 * MS3[3, 3]));
    ME1(1,6)=1/(2*MS3(3,3));

	%etaxxy: = simplify(2 * MS3[1, 3] * muxy); 
    ME1(1,7)=2*MS3(1,3)*ME1(1,6);

	%etayxy: = simplify(2 * MS3[2, 3] * muxy); 
    ME1(1,8)=2*MS3(2,3)*ME1(1,6);

	%etaxyx: = simplify(MS3[3, 1] * Ex); 
    ME1(1,9)=MS3(3,1)*ME1(1,2);

	%etaxyy: = simplify(MS3[3, 2] * Ey);
    ME1(1,10)=MS3(3,2)*ME1(1,3);
end

