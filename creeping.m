%Cavalleri Andrea 
%1853313
function Tesina_Creeping
%Modello di traffico eterogeneo con due classi di veicoli che permette
%sorpassi e creeping in flussi di traffico altamente eterogenei.
clc
close all
clear
%%
%Impostazione problema
modello=input('\nScegli il modello da utilizzare:\n1. Modello con creeping\n2. Modello senza creeping\n');
while modello~=1 && modello~=2
            fprintf('\nInput non valido! \nDigitare il numero a sinistra del modello scelto\n')
            modello=input('\nScegli il modello da utilizzare:\n1. Modello con creeping\n2. Modello senza creeping\n');
end
r1max=1.8; %densità massima della prima classe di veicoli (veicoli piccoli)
if modello==1
    r2max=1.0; %densità massima della seconda classe di veicoli (veicoli lunghi) nel modello con creeping
elseif modello==2
    r2max=1.8; %densità massima della seconda classe di veicoli (veicoli lunghi) nel modello senza creeping
end
dato=input('\nScegli uno tra i tre esempi:\n1. Overtaking\n2. Creeping\n3. Overtaking+Creeping\n');
while dato~=1 && dato~=2 && dato~=3
            fprintf('\nInput non valido! \nDigitare il numero a sinistra dell esempio scelto\n')
            dato=input('\nScegli uno tra i tre esempi:\n1. Overtaking\n2. Creeping\n3. Overtaking+Creeping\n');
end
if modello==1
    v1max=1.8;
    v2max=1.8;
elseif modello==2
    if dato==1
        v1max=1.8;
        v2max=1.0;
    elseif dato==2
        v1max=1.8;
        v2max=1.0;
    elseif dato==3
        v1max=1.8;
        v2max=1.0;
    end
end
vmax=max(v1max,v2max); %velocità massima
Q1= @(r1,r2) r1.*v1max.*(1-(r1+r2)./r1max); %Flusso della prima classe di veicoli
if modello==1
    Q2= @(r1,r2) max(r2.*v2max.*(1-(r1+r2)./r2max),0); %Flusso della seconda classe di veicoli
elseif modello==2
    Q2= @(r1,r2) r2.*v2max.*(1-(r1+r2)./r2max);
end

%dQ1= @(r1,r2) (vmax/r1max).*(r1max-2.*(r1+r2)); %Derivata (rispetto ad r1) di Q1
%dQ2= @(r1,r2) max((vmax/r2max).*(r2max-2.*(r1+r2)),0); %Derivata (rispetto a r2) di Q2
switch dato
    case 1 %overtaking
        r10= @(x) 0.9*(x<=10);
        r20= @(x) 0.9*((x<=20) & (x>=11));
    case 2 %Creeping
        r10= @(x) 0.7*(x<=19);
        r20= @(x) 0.7*(x>=20);
    case 3 
        r10= @(x) 0.9*(x<=10);
        r20= @(x) 0.9*((11<=x)&(x<=25));
end
x0=0; xf=50;
y0=0; yf=2;
%%
%Discretizzazione in spazio
Dx=0.05; %passo spaziale
M=abs(xf-x0)/Dx; %Numero di celle in spazio
nodi=x0+(0:M)'.*Dx;
C=(nodi(1:end-1)+nodi(2:end))/2; %M centri-cella
%%
%Verifica CFL
lambda=1/vmax; 
Dt=0.7*lambda*Dx; %passo temporale
T=65;
%%
for i=1:T
%Scelta e calcolo soluzione numerica
[U1,U2]=Godunov2classi(r10(C),r20(C),Q1,Q2,r1max,r2max,Dx,Dt,i,0,dato);
[~,NT1]=size(U1); NT1=NT1-1; %nodi in tempo (aggiorno NT perché Dt potrebbe essere stato modificato)
[~,NT2]=size(U2); NT2=NT2-1; 
%%
%Discretizzazione stato iniziale (per il grafico)
U10=[r10(C),zeros(M,NT1)]; %inizializzazione matrice soluzione esatta
U20=[r20(C),zeros(M,NT2)];
for n=1:NT1
    U10(:,n+1)=r10(centri,sum(U(end,1:n+1))); %l'ultima riga di U contiene i Dt
end
for n=1:NT2
    U20(:,n+1)=r20(centri,sum(U(end,1:n+1)));
end
%%
%Rappresentazione grafico
grafico(U1,U2,C,dato,modello,U10,U20,x0,xf,y0,yf,NT2);
end
end
function [U1,U2]=Godunov2classi(U1old,U2old,Q1,Q2,r1max,r2max,Dx,Dt,T,t,dato)
%Approssima ricorsivamente la soluzione della legge di conservazione.
%associata a Q1 e Q2 tramite il metodo di Godunov con passi DX eDT
%Alla prima iterazione t=0 e U1old e U2old sono due vettori colonna contenenti 
%le approssimazioni dei dati iniziali.
%t è il tempo al quale si approssima la soluzione ad ogni iterazione.
%Restituisce due matrici U1 e U2 in cui ogni colonna corrisponde alla soluzione
%approssimata a un certo un passo temporale.
%Le condizioni al bordo vengono scelte a seconda di quale esempio è preso
%in considerazione indicato dal valore di dato.
%U1 è la soluzione relativa alla prima classe di veicoli e U2 quella
%relativa alla seconda classe di veicoli.
    if t>T
        U1=U1old;
        U2=U2old;
        return;
    end
    U1new=zeros(size(U1old)); 
    U2new=zeros(size(U2old));
    
    U1new(2:end-1)=U1old(2:end-1)-Dt/Dx*(Fg(U1old(2:end-1),U1old(3:end),U2old(2:end-1),U2old(3:end),Q1,Q2,r1max,r2max,1)-Fg(U1old(1:end-2),U1old(2:end-1),U2old(1:end-2),U2old(2:end-1),Q1,Q2,r1max,r2max,1));
    U2new(2:end-1)=U2old(2:end-1)-Dt/Dx*(Fg(U1old(2:end-1),U1old(3:end),U2old(2:end-1),U2old(3:end),Q1,Q2,r1max,r2max,2)-Fg(U1old(1:end-2),U1old(2:end-1),U2old(1:end-2),U2old(2:end-1),Q1,Q2,r1max,r2max,2));
    
    if dato==1
    %Condizioni al bordo per l'overtaking
        U1new(1) = U1old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,1)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,1)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U1new(end) = U1old(end)-Dt/Dx*(Fg(U1old(end-1),0,U2old(end-1),0,Q1,Q2,r1max,r2max,1)-Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,1)); % Downstream della regione di studio vuoto (veicoli liberi di uscire dal dominio)
        U2new(1) = U2old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,2)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,2)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U2new(end) = U2old(end)-Dt/Dx*(Fg(U1old(end-1),0,U2old(end-1),0,Q1,Q2,r1max,r2max,2)-Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,2)); % Downstream della regione di studio vuoto (veicoli liberi di uscire dal dominio)
    elseif dato==2
    %condizioni al bordo per il creeping
        U1new(1) = U1old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,1)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,1)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U1new(end) = U1old(end)+Dt/Dx*(Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,1)); % Downstream outflow uguale a zero (veicoli fermati dal semaforo rosso)
        U2new(1) = U2old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,2)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,2)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U2new(end) = U2old(end)+Dt/Dx*(Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,2)); % Downstream outflow uguale a zero (veicoli fermati dal semaforo rosso)
    elseif dato==3
        U1new(1) = U1old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,1)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,1)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U1new(end) = U1old(end)+Dt/Dx*(Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,1)); % Downstream outflow uguale a zero (veicoli fermati dal semaforo rosso)
        U2new(1) = U2old(1)-Dt/Dx*(Fg(U1old(2),U1old(3),U2old(2),U2old(3),Q1,Q2,r1max,r2max,2)-Fg(0,U1old(2),0,U2old(2),Q1,Q2,r1max,r2max,2)); % Upstream inflow uguale a zero (nessun nuovo veicolo entra)
        U2new(end) = U2old(end)+Dt/Dx*(Fg(U1old(end-2),U1old(end-1),U2old(end-2),U2old(end-1),Q1,Q2,r1max,r2max,2)); % Downstream outflow uguale a zero (veicoli fermati dal semaforo rosso)
    end
    %Chiamata ricorsiva
    [U1, U2] = Godunov2classi(U1new, U2new, Q1, Q2, r1max, r2max, Dx, Dt, T, t + Dt,dato);
end
function F=Fg(Ul,Ur,Ul2,Ur2,Q1,Q2,r1max,r2max,class)
    if class==1 
        S=Q1(Ul,Ul2).*(Ul<=(r1max-Ul2)/2)+Q1((r1max-Ul2)/2,Ul2).*(Ul>(r1max-Ul2)/2);
        R=Q1(Ur,Ur2).*(Ur>(r1max-Ur2)/2)+Q1((r1max-Ur2)/2,Ur2).*(Ur<=(r1max-Ur2)/2);
        F=min(S,R);
    elseif class==2 
        S=Q2(Ul,Ul2).*(Ul2<=(r2max-Ul)/2)+Q2(Ul,(r2max-Ul)/2).*(Ul2>(r2max-Ul)/2);
        R=Q2(Ur,Ur2).*(Ur2>(r2max-Ur)/2)+Q2(Ur,(r2max-Ur)/2).*(Ur2<=(r2max-Ur)/2);
        F=min(S,R);
    end
end
%%
%Function 3 (plot)
function grafico(U1,U2,centri,dato,modello,U1in,U2in,x0,xf,y0,yf,NT)
%Plotta la soluzione numerica della classe di veicoli class di dato memorizzata nella matrice U.
%centri è un vettore con i centri delle celle della griglia spaziale.
%NT è uguale al numero di nodi temporali in cui la soluzione è approssimata.
%il dato iniziale è memorizzato nella matrice Uin (altrimenti Uex=[]).
%x0,xf,y0,yf sono i valori per gli estremi del grafico.
for n=1:NT+1
    subplot(2,1,1);
    %Grafico veicoli leggeri
    plot(centri,U1(:,n),'b-*',centri,U1in(:,n)); %grafico soluzione approssimata (U1) e dati iniziali (U1in)
    %Titolo grafico (a seconda dell'esempio preso in considerazione
    if modello==1
        if dato==1
            title('Creeping Model for Heterogeneous Traffic Flow: Overtaking')
        elseif dato==2
            title('Creeping Model for Heterogeneous Traffic Flow: Creeping')
        elseif dato==3
            title('Creeping Model for Heterogeneous Traffic Flow: Overtaking+Creeping')
        end
    elseif modello==2
            if dato==1
                title('n-Population Model for Heterogeneous Traffic Flow: Overtaking')
            elseif dato==2
                title('n-Population Model for Heterogeneous Traffic Flow: Creeping')
            elseif dato==3
                title('n-Population Model for Heterogeneous Traffic Flow: Overtaking+Creeping')
            end
    end
    
    
    %Grafico veicoli leggeri (titolo, legenda e assi)
    subtitle('Veicoli Piccoli')
    legend('Densità dei veicoli piccoli','Stato iniziale','Location','northwest')
    xlabel('Spazio x');
    xlim([x0 xf]) %estremi delle x
    ylabel('\rho_1')
    ylim([y0 yf])%estremi delle y
    drawnow
    %Grafico veicoli pesanti
    subplot(2,1,2);
    plot(centri,U2(:,n),'b-*',centri,U2in(:,n)); %grafico soluzione approssimata (U2) e dati iniziali (U2in)
    %Grafico veicoli pesanti (titolo, legenda e assi)
    subtitle('Veicoli Grandi')
    legend('Densità veicoli grandi','Stato iniziale','Location','northwest')
    xlabel('Spazio x');
    xlim([x0 xf]) %estremi delle x
    ylabel('\rho_2')
    ylim([y0 yf])%estremi delle y
    drawnow
end
end
