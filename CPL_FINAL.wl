(* ::Package:: *)

SetDirectory["C:\Users\Diego\Desktop\TFG\Ajustes"];
(*//////////////////////////////////////////////*)
dataSNIa=Import["uniondos.dat"];          (* Archivo de datos de supernovas*)
ndatsnia=557;                               (* N\[UAcute]mero de datos *)
Do[zsnia[i]=dataSNIa[[i,2]],{i,1,ndatsnia}];     (* Lee la columna del redshift*)
Do[mo[i]=dataSNIa[[i,3]],{i,1,ndatsnia}];    (* Lee la columna del m\[OAcute]dulo de la distancia *)
Do[errmo[i]=dataSNIa[[i,4]],{i,1,ndatsnia}]; (* Lee el error de la distancia *)
CC=Sum[1/errmo[i]^2,{i,1,ndatsnia}];
(*//////////////////////////////////////////////*)
dataBAO=Import["BAO.dat"];          (* Archivo de datos*)
ndatbao=17;                               (* N\[UAcute]mero de datos *)
Do[zbao[i]=dataBAO[[i,2]],{i,1,ndatbao}];     (* Lee la columna del redshift*)
Do[dz[i]=dataBAO[[i,3]],{i,1,ndatbao}];    (* Lee la columna del par\[AAcute]metro dz(z) *)
Do[errdz[i]=dataBAO[[i,4]],{i,1,ndatbao}]; (* Lee el error de dz *)
Do[Az[i]=dataBAO[[i,5]],{i,1,ndatbao}];    (* Lee la columna del par\[AAcute]metro A(z) *)
Do[errAz[i]=dataBAO[[i,6]],{i,1,ndatbao}]; (* Lee el error de A *)
(* Matrices de covariancia *)
covarinvdz=Import["covariancedz.txt","Table"];
covarinvaz=Import["covarianceaz.txt","Table"];
covarinv=Import["covarianceCMB.txt","Table"];
(*//////////////////////////////////////////////*)
dataHZ=Import["Hdatos.dat"];          (* Archivo de datos*)
ndathz=31;                               (* N\[UAcute]mero de datos *)
Do[zhz[i]=dataHZ[[i,1]],{i,1,ndathz}];   (* Lee la columna del redshift*)
Do[H[i]=dataHZ[[i,2]],{i,1,ndathz}];    (* Lee la columna del m\[OAcute]dulo de la distancia *)
Do[errH[i]=dataHZ[[i,3]],{i,1,ndathz}]; (* Lee el error de la distancia *)
(*//////////////////////////////////////////////*)
rpl=1.7448;
lapl=301.46;
\[Omega]bpl=0.0224;
c:=300000; (*Velocidad de la luz*)
cs:=173205 (* Velocidad del sonido en la recombinaci\[OAcute]n*);
acmb:=1/(1089.90+1);


Clear[jfinal,ff1,chisqsnU,n,b,d,a,w,H0,chisqsnUbis,\[Delta]1bis,\[Delta]2bis,\[Delta]3bis,\[Delta]5bis,\[CapitalOmega]m,c,s,FF1,\[Delta]4bis,R,nc,ncbis,integralw,dzteor]
wx:=1 (* N\[UAcute]mero de cadenas de Markov *)
ModuloMCMC2[jfinalp_,np_]:=Module[{jfinal=jfinalp,n=np}, (* jfinal numero de iteraciones, np numero de cadenas de Markov *)
Clear[ff1];
For[j=1,j<jfinal+1,
Clear[Aleatorio,hh,hh1,hh2,Aleatorio1,Aleatorio2,Aleatorio3,sol,hubble,h,solution,sol1,Eq,chi,int,integralw];
Aleatorio=RandomReal[{0,1}];
chisqsnUbis[1]=1;
(* Valores iniciales*)
Initial[1]:=0.30;
Initial[2]:=-0.9;
Initial[3]:= -0.1;
 
 (* Distribuciones iniciales *)

l[j-1]=If[j>1,RandomVariate[NormalDistribution[\[Delta]3bis[j-1],0.02]],Initial[1]]; (* Saltos aleatorios en cada iteraci\[OAcute]n, valor dentro de una gaussiana con mu y sigma*)
ll[j-1]=If[j>1,RandomVariate[NormalDistribution[\[Delta]4bis[j-1],0.02]],Initial[2]]; (* Saltos aleatorios en cada iteraci\[OAcute]n, valor dentro de una gaussiana con mu y sigma*)
lll[j-1]=If[j>1,RandomVariate[NormalDistribution[\[Delta]5bis[j-1],0.02]],Initial[3]]; (* Saltos aleatorios en cada iteraci\[OAcute]n, valor dentro de una gaussiana con mu y sigma*)
\[CapitalOmega]m[j]=l[j-1];
w[j]=ll[j-1];
wprima[j]=lll[j-1];

c:=300000; (*Velocidad de la luz*)
Clear[sol,hubble,h,solution,solutionhz,sol1,Eq,chi,chiBAO,chiCMB,chihz,chisnia,int,intcmb,integralw,compchi,jtotal,dzteor,azteor,\[Omega]bteor,rteor,lateor,soundhorizon,solutiona,solutionz,solutionacmb,solutionzcmb];

\[CapitalOmega]r=9.23*10^(-5);
aeq[j]:=\[CapitalOmega]r/\[CapitalOmega]m[j];

solutionz[z_]:=(\[CapitalOmega]r (1+z)^4 +  \[CapitalOmega]m[j] (1+z)^3+(1-\[CapitalOmega]m[j]-\[CapitalOmega]r) (1+z)^(3(1+w[j]+wprima[j])Exp[-3wprima[j]z/(1+z)]))^(1/2);   (* Ecuaciones de Friedmann, calcula el par\[AAcute]metro de Hubble para cada valor de \[CapitalOmega]m *)
solutiona[a_]:= (\[CapitalOmega]m[j](a+aeq[j]) a^(-4)+(1-\[CapitalOmega]m[j]-\[CapitalOmega]r) a^(-3(1+w[j])))^(1/2);

(*///////////////////// SNe Ia /////////////////////////*)

Do[intsnia[i]=NIntegrate[1/solutionz[x],{x,0,zsnia[i]},Method->{Automatic, "SymbolicProcessing" -> False}],{i,1,ndatsnia}];
chisnia=Sum[(mo[i]-(5*Log[10,(1+zsnia[i])*intsnia[i]]))^2/errmo[i]^2,{i,1,ndatsnia}]-(Sum[(mo[i]-(5*Log[10,(1+zsnia[i])*intsnia[i]]))/errmo[i]^2,{i,1,ndatsnia}])^2/CC;

(*///////////////////// H(z) /////////////////////////*)

solutionhz[z_]:=67.3*(\[CapitalOmega]m[j] (1+z)^3+(1-\[CapitalOmega]m[j]) (1+z)^(3(1+w[j]+wprima[j])Exp[-3wprima[j]z/(1+z)] ))^(1/2);  (* Ecuaciones de Friedmann, calcula el par\[AAcute]metro de Hubble para cada valor de \[CapitalOmega]m *)
   
chihz= Sum[(solutionhz[zhz[i]]-H[i])^2/(errH[i])^2,{i,1,ndathz}]; (* F\[OAcute]rmula del chi 2*)


(*///////////////////// BAO /////////////////////////*)

Do[int[i]=NIntegrate[1/solutionz[x],{x,0,zbao[i]},Method->{Automatic, "SymbolicProcessing" -> False}],{i,1,ndatbao}]; (* Integra 1/H entre 0 y cada uno de los redshifts*)

Do[If[dz[i]==0, dzteor[i]=0, dzteor[i]= 10457/c*(zbao[i]*(int[i])^2/(solutionz[zbao[i]]))^(-1/3),{i,1,ndatbao}],{i,1,ndatbao}];
Do[If[Az[i]==0, azteor[i]=0 , azteor[i]= Sqrt[\[CapitalOmega]m[j]]/(c*zbao[i]) * (c*zbao[i]*(c*(1+zbao[i])*int[i])^2/((1+zbao[i])^2*solutionz[zbao[i]]))^(1/3)],{i,1,ndatbao}];

chiBAO=(Array[dz,ndatbao]-Array[dzteor,ndatbao]).covarinvdz.(Array[dz,ndatbao]-Array[dzteor,ndatbao])+(Array[Az,ndatbao]-Array[azteor,ndatbao]).covarinvaz.(Array[Az,ndatbao]-Array[azteor,ndatbao]);


(*///////////////////// CMB /////////////////////////*)
solutionzcmb[z_]:=67.3*(\[CapitalOmega]r (1+z)^4 +  \[CapitalOmega]m[j] (1+z)^3+(1-\[CapitalOmega]m[j]-\[CapitalOmega]r)(1+z)^(3(1+w[j]+wprima[j])Exp[-3wprima[j]z/(1+z)]))^(1/2);  (* Ecuaciones de Friedmann, calcula el par\[AAcute]metro de Hubble para cada valor de \[CapitalOmega]m *)
solutionacmb[a_]:= 67.3*(\[CapitalOmega]r a^(-4)  + \[CapitalOmega]m[j] a^(-3)+(1-\[CapitalOmega]m[j]-\[CapitalOmega]r) a^(-3(1+w[j])))^(1/2);

intcmb=NIntegrate[1/solutionzcmb[x],{x,0,1089.9},Method->{Automatic, "SymbolicProcessing" -> False}];
soundhorizon=1/Sqrt[3]* NIntegrate[1/(x^2 solutionacmb[x]*Sqrt[1+x*5.687*10^(-10)]),{x,0,acmb},Method->{Automatic, "SymbolicProcessing" -> False}];

rteor=67.3*Sqrt[\[CapitalOmega]m[j]]*intcmb;
lateor=\[Pi]*intcmb/soundhorizon;
\[Omega]bteor=\[Omega]bpl-0.002109*(rteor-rpl)-0.0004341*(lateor-lapl);


chiCMB=({rteor,lateor,\[Omega]bteor}-{rpl,lapl,\[Omega]bpl}).covarinv.({rteor,lateor,\[Omega]bteor}-{rpl,lapl,\[Omega]bpl});
 
(*////////////////////////////////////////////////////////////////////////////////////////////////////*)
chi= chihz+chisnia+chiBAO;

chisqsnU[j]=If[MatchQ[chi,_Real],chi,10^10]; (* Descarta valores no reales y les da un chi gigante*)

compchi=Exp[-chisqsnU[j]/2]/(Exp[-chisqsnUbis[j]/2]); (* Para comparar con el likelyhood *)

chisqsnUbis[j+1]=If[j>1,If[compchi>Aleatorio,chisqsnU[j],chisqsnUbis[j]],chisqsnU[j]]; 

\[Delta]3bis[j]=If[j>1,If[compchi>Aleatorio,\[CapitalOmega]m[j],\[Delta]3bis[j-1]],\[CapitalOmega]m[j]]; (* Para solo coger los mas peque\[NTilde]os significativamente *)
\[Delta]4bis[j]=If[j>1,If[compchi>Aleatorio,w[j],\[Delta]4bis[j-1]],w[j]]; 
\[Delta]5bis[j]=If[j>1,If[compchi>Aleatorio,wprima[j],\[Delta]5bis[j-1]],wprima[j]]; 

ncbis[j]=If[jfinal>j>1,If[compchi>Aleatorio,1,0],0];
ncbis[jfinal]:=1;

nc[j]=If[j>1,If[compchi>Aleatorio||j==jfinal,Sum[ncbis[l],{l,1,j}],0],0];
jj[nc[j]]=If[N[nc[j]]>1,j-1-Sum[jj[l],{l,1,N[nc[j]]-1}],j-1];

Print[{j,\[Delta]3bis[j],\[Delta]4bis[j],\[Delta]5bis[j],chisqsnUbis[j+1],nc[j],jj[nc[j]]}];

If[ncbis[j]>0||j==jfinal,ff1[nc[j],1]=N[jj[nc[j]]];ff1[nc[j],2]=N[chisqsnUbis[j]];ff1[nc[j],3]=N[\[Delta]3bis[j-1]];ff1[nc[j],4]=N[\[Delta]4bis[j-1]];ff1[nc[j],5]=N[\[Delta]5bis[j-1]]];


Clear[Aleatorio,hh,hh1,hh2,Aleatorio1,Aleatorio2,Aleatorio3,sol,hubble,h,solution,sol1,Eq,chi,int];
j++];
jtotal=Count[Table[ncbis[ww],{ww,1,jfinal}],1];
FF1=Array[ff1,{jtotal,5}];
If[n==1,Export["AjusteCPL.dat",FF1],If[n==2,Export["Pruebaw001E1Bbis.dat",FF1],Export["Pruebaw001E1Cbis.dat",FF1]]];
Clear[FF1,ff1,q,l,m,b,\[CapitalOmega]m,s,\[Alpha],\[Delta]2bis,\[Delta]3bis,\[Delta]4bis,\[Delta]5bis,chisqsnU,chisqsnUbis,ncbis,nc,jtotal]
];


Do[ModuloMCMC2[30000,wx],{wx,1,1,1}](* Se suelen correr mas de 10000*)





