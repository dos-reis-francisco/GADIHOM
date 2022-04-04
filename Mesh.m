%% Mesh
% create the mesh of the lattice 
% sub module for inverse homogenization code
% Dos Reis F.
% 03.02.2021

function Mesh(seed)
    global L1 L2 nbeams nnodes nodes Ob Eb delta1 delta2 Lb;
    NL=seed; NC=seed;
    L1s=L1/seed;L2s=L2/seed;
    nbeams=6*NL*NC;
    nnodes=NL*NC*2;
    nodes=zeros(nnodes,2);
    Ob=zeros(nbeams,1);
    Eb=zeros(nbeams,1);
    delta1=zeros(nbeams,1);
    delta2=zeros(nbeams,1);
    Lb=zeros(1,nbeams);
    Ldiag=1/2*sqrt(L1s^2+L2s^2);
    Lhor=L1s; Lvert=L2s;
    beam=1;
    for L=1:NL
        for C=1:NC
            type=1;
            [numbernode,d1,d2,node]=NumNode(type,NC,NL,C,L,L1s,L2s);
            n1=numbernode;nodes(numbernode,1:2)=node;
            type=2;
            [numbernode,d1,d2,node]=NumNode(type,NC,NL,C,L,L1s,L2s);
            n2=numbernode;nodes(numbernode,1:2)=node;
			% first beam
            Ob(beam)=n1;Eb(beam)=n2;
            Lb(1,beam)=Ldiag;
            beam=beam+1;
            type=3;
            [numbernode,d1,d2,node]=NumNode(type,NC,NL,C,L,L1s,L2s);
            n3=numbernode;
            delta1(beam)=d1;delta2(beam)=d2;
			% second beam
            Ob(beam)=n2;Eb(beam)=n3;
            Lb(1,beam)=Ldiag;
            beam=beam+1;
            delta1(beam)=d1;delta2(beam)=d2;
			% third beam
            Ob(beam)=n1;Eb(beam)=n3;
            Lb(1,beam)=Lhor;
            beam=beam+1;
            type=4;
            [numbernode,d1,d2,node]=NumNode(type,NC,NL,C,L,L1s,L2s);
            n4=numbernode;
            delta1(beam)=d1;delta2(beam)=d2;
			% thourth beam
            Ob(beam)=n2;Eb(beam)=n4;
            Lb(1,beam)=Ldiag;
            beam=beam+1;
            type=5;
            [numbernode,d1,d2,node]=NumNode(type,NC,NL,C,L,L1s,L2s);
            n5=numbernode;
            delta1(beam)=d1;delta2(beam)=d2;
			% fith beam
            Ob(beam)=n2;Eb(beam)=n5;
            Lb(1,beam)=Ldiag;
            beam=beam+1;
            delta1(beam)=d1;delta2(beam)=d2;
			% sixth beam
            Ob(beam)=n1;Eb(beam)=n5;
            Lb(1,beam)=Lvert;
            beam=beam+1;
        end
    end
end

%% function NumNode
% return the number node and deltai values 
% type : 1 ou 2 suivant le numéro du noeuds dans la cellule
% NC : number of cell in columns
% type is a code : 5    4
%                     2
%                  1    3
% L1s, L2s : increment x, y for one cell
function [numbernode,delta1,delta2,node]=NumNode(type,NC,NL,C,L,L1s,L2s)
    node=zeros(2,1);
    if type==1 
        numbernode=(L-1)*NC*2+C*2-1;
        delta1=0;
        delta2=0;
        node(1,1)=(C-1)*L1s;node(2,1)=(L-1)*L2s;
    elseif type==2
        numbernode=(L-1)*NC*2+C*2;
        delta1=0;
        delta2=0;       
        node(1,1)=(C-0.5)*L1s;node(2,1)=(L-0.5)*L2s;
    elseif type==3
        if (C+1)<=NC
            numbernode=(L-1)*NC*2+(C+1)*2-1;
            delta1=0;
            delta2=0; 
            node(1,1)=C*L1s;node(2,1)=(L-1)*L2s;
        else
            numbernode=(L-1)*NC*2+1;
            delta1=1;
            delta2=0; 
            node(1,1)=0;node(2,1)=(L-1)*L2s;
        end
    elseif type==4
        if (C+1)<=NC
            C=C+1;
            delta1=0;
            node(1,1)=C*L1s;
        else
            C=1;
            delta1=1;
            node(1,1)=0;
        end
        if (L+1)<=NL
            L=L+1;
            delta2=0;
            node(2,1)=L*L2s;
        else
            L=1;
            delta2=1;
            node(2,1)=0;
        end
        numbernode=(L-1)*NC*2+C*2-1;
    elseif type==5
        delta1=0;
        if (L+1)<=NL
            L=L+1;
            delta2=0;
            node(2,1)=L*L2s;
        else
            L=1;
            delta2=1;
            node(2,1)=0;
        end
        node(1,1)=(C-1)*L1s;
        numbernode=(L-1)*NC*2+C*2-1;        
    end
end