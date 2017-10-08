#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

//Domi Box
struct Box{
	float xmin,ymin,zmin;  //ta x,y,z tou koutiou (min)
	int specid,genid,neighbors,*neighborID,qpoints,cpoints,*qpointsID,*cpointsID;  //EidikoID, genikoID, arithmos geitonwn, ID geitonwn, arithmos shmeiwn Q kai C kai ids twn Q kai C
};
struct Box *Boxes,*neighBoxes;                       //Ta Boxes kai ta geitonika Boxes gia kathe process

//Domi Point
struct Point{
	float x,y,z;                      //Oi syntetagmenes tou shmeiou x,y,z
        float minDist;                    //H apostasi apo to kontinotero C (gia ta Q)
	int id,genid,boxid,nearestid;     //eidikoID, genikoID, genikoID koutiou sto opoio anhkei, genikoID kontinotero shmeiou C (gia ta Q)
};
struct Point *Q,*C,*neighC;                          //Ta shmeia Q kai C gia kathe process

double **Qin,**Cin,**recQ,**recC,**sendC,**recC2;    //Oi buffers twn Q kai C pou tha xrhsimopoihthoun gia tis epikoinwnies metaksy process
int neighbors,*neighborID,neighFlag[26],*neighFlagID; //Plithos Geitonwn process, ID geitonwn, Flags geitonwn(apo poies kateuthynseis exw geitona), FlagID geitonwn
int n,m,k,pn,pm,pk,n_2,m_2,k_2;       //Diastaseis Grid Boxes, Grid process, Boxes/process
int NoP,BoxNum;                       //Arithmos processes, Boxes
int qCnt=0,cCnt=0,cCnt2,cSendCnt[26],cRecCnt[26],neighBoxesCnt=0;  //Plithos shmeiwn Q,C,geitonikwn C, megethi buffers, arithmos geitonikwn Boxes
int processX,processY,processZ;       //Oi syntetagmenes tou process sto Grid twn process
double pstartX,pstartY,pstartZ;       //To xmin,ymin,zmin tou kathe process
int selfPID;     //To PID mou
double t1,t2;    //Metavlites gia ti metrisi tou xronou ekteleshs


void pGridDimensions(int num);           //Vriskei tis diastaseis tou Grid twn processes
void BoxGridDimensions(int num);         //Vriskei tis diastaseis tou Grid twn koutiwn
int BoxIdentify(struct Point p);         //Vriskei se poio kouti anhkei to kathe shmeio
void BoxCreate();                        //Dhmiourgei ta koutia
void BoxAssign();                        //Kanei assign ta shmeia mou sta antistoixa Boxes
double distance(struct Point p1,struct Point p2); //Ypologizei tin apostasi metaksy 2 shmeiwn ston 3d xwro
int ProcIdentify(double p[3]);           //Vriskei se poio process anhkei to shmeio
void findNearest(int boxid);             //Synarthsh ypologismou tou kontinoterou shmeiou C gia kathe shmeio Q tou koutiou boxid 
void findPNeighbors(int x,int y,int z);  //Ypologizei tis geitonikes process gia ti dothousa process
int spec2gen(int spec);                  //Metatropi tou eidikou ID tou box enos process se geniko ID oloklirou tou Grid
int gen2spec(int gen);                   //Metatropi tou genikou ID tou box tou Grid se eidiko ID gia to sygkekrimeno process
void findSendPoints();                   //Vriskei ta shmeia C pou tha steilei stis geitonikes process
void BoxNeighCreate(int cpoints);        //Dimiourgei ta geitonika oriaka koutia pou eginan receive
int ifMine(int boxid);                   //Elegxei an to boxid einai box tis process
int findIndexNeigh(int box);             //Vriskei to index ston pinaka neighBoxes gia to sygkekrimeno box
void checkAndCorrect();                  //Elegxei kai diorthonei alloiwseis kata tin lipsi twn prwtwn minimatwn
void checkAndCorrect2();                 //Elegxei kai diorthonei alloiwseis kata tin lipsi twn deuterwn minimatwn

int main(int argc,char *argv[]){
        int i,j,o;                       //Genikoi counters gia ta for loops
	int Nq,Nc;                       //O arithmos ton shmeiwn mou
        
        Nq=atoi(argv[1]);                //Pernao ta 3 arguments mou, Arithmos shmeion Nq, arithmos koution BoxNum, arithmos diergasiwn NoP
        BoxNum=atoi(argv[2]);
        NoP=atoi(argv[3]);

        MPI_Request mpireq[2*(NoP-1)];  //Ftiaxno ton Communicator mou, ksekinao ta processes, kai dimiourgw voithitika MPIrequests kai MPIstats
        int reqcnt=0;                   //Counter twn MPIrequest mou
        MPI_Status mpistat[2*(NoP-1)];
        MPI_Init( &argc, &argv);
        MPI_Comm_size( MPI_COMM_WORLD, &NoP );
        MPI_Comm_rank( MPI_COMM_WORLD, &selfPID );
       
        
   
        pGridDimensions(NoP);           //Dimiourgw to Grid twn processes (pn,pm,pk)
        processX=(selfPID%pn);          //Ypologizw tis syntetagmenes tou process sto Grid twn processes (processX,processY,processZ),   
        pstartX=(double) processX/pn;   //kathos kai ta actual x, y, z tou process mou (pstartX,pstartY,pstartZ)
        processY=((selfPID/pn)%pm);
        pstartY=(double) processY/pm;
        processZ=(selfPID/pn)/pm;    
        pstartZ=(double) processZ/pk;  
	
        Nq=pow(2,Nq);                   //Ypswnw to 2 eis tin Nq gia na parw gia ton synoliko arithmo twn shmeiwn Q mou
        Nc=Nq;                          //Idios arithmos gia to synolo C
    
	Qin=(double **)malloc(NoP*sizeof(double *)); //Kanw allocate mnhmhs gia tous arxikous pinakes twn Q kai C pou dimiourgw
        Cin=(double **)malloc(NoP*sizeof(double *)); //Kathe seira antistoixei sto antistoixo process
        for (i=0;i<NoP;i++){
             Qin[i]=(double *)malloc(((Nq/NoP)/NoP*4)*sizeof(double)); 
             Cin[i]=(double *)malloc(((Nc/NoP)/NoP*4)*sizeof(double));
        }
        
        srand(time(NULL)/(selfPID+1)); //ftiaxno to seed mou
        int procBelong;                //prosorini metavliti pou tha anaferetai sthn process pou anhkei to kathe shmeio
        int NumPointsQ[NoP],NumPointsC[NoP],recNumQ[NoP],recNumC[NoP]; //Arithmos shmeiwn pou stelnw/lamvanw gia kathe diergasia
        for (i=0;i<NoP;i++){           //Ta mhdenizw
              NumPointsQ[i]=0;
              NumPointsC[i]=0;
              recNumQ[i]=0;
              recNumC[i]=0;
        }
        double t[3];                   //temp double gia tis syntetagmenes twn shmeiwn
        for (i=0;i<(Nq/NoP);i++){      //Kathe process dhmiourgei Nq/NoP shmeia C kai Q
                for (j=0;j<3;j++){     //Meso tis rand()/RAND_MAX pairnw syntetagmenh sto [0,1)
	              t[j]=((double) rand() / (RAND_MAX));
                      if (t[j]==1)  t[j]=((double) rand() / (RAND_MAX));
                }
                procBelong=ProcIdentify(t);  //Vrisko se poio process anhkei
                if (NumPointsQ[procBelong]>=((Nq/NoP)/NoP*4)){  //Elegxw an xreiazomai Reallocate mnhmhs
                     Qin[procBelong]=(double *)realloc(Qin[procBelong],(NumPointsQ[procBelong]+4)*sizeof(double)); 
                }
                for (j=0;j<3;j++){     //Pernaw tis 3 syntetagmenes x,y,z ws synexomenes sti thesi Qin[procbelong][j]
                     Qin[procBelong][NumPointsQ[procBelong]++]=t[j];
                }
                Qin[procBelong][NumPointsQ[procBelong]++]=i+selfPID*(Nq/NoP); //Sth synexeia pernaw kai to genikoID tou shmeiou (monadiko)
           	for (j=0;j<3;j++){     //Kanw akrivws to idio gia ton pinaka Cin pou anaferetai sto synolo C 
	              t[j]=((double) rand() / (RAND_MAX));
                      if (t[j]==1)  t[j]=((double) rand() / (RAND_MAX));
                }
                procBelong=ProcIdentify(t);
                if (NumPointsC[procBelong]>=((Nc/NoP)/NoP)*4){
                     Cin[procBelong]=(double *)realloc(Cin[procBelong],(NumPointsC[procBelong]+4)*sizeof(double)); 
                }
                for (j=0;j<3;j++){
                     Cin[procBelong][NumPointsC[procBelong]++]=t[j];
                }
                Cin[procBelong][NumPointsC[procBelong]++]=i+selfPID*(Nc/NoP);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);  //Kanw MPI_Barrier gia na metrisw mesw ths MPI_Wtime() ton xrono ekkinishs tou algorithmou
        t1=MPI_Wtime();
        
        recQ=(double **)malloc(NoP*sizeof(double *)); //Kanw allocate theseis gia tous buffers tou receive apo kathe process
        recC=(double **)malloc(NoP*sizeof(double *));
        for (i=0;i<NoP;i++){           //Stelnw se kathe process prwta to megethos tou buffer kai sth synexeia ton buffer me ta shmeia Q kai C pou ths anhkoun
              if (i!=selfPID){
                  MPI_Isend(&NumPointsQ[i],1,MPI_INT,i,50,MPI_COMM_WORLD, &mpireq[i]);  //Diaforetiko tag se kathe mynhma (edw den einai aparaithto)
                  MPI_Isend(&NumPointsC[i],1,MPI_INT,i,150,MPI_COMM_WORLD,&mpireq[i]);
                  MPI_Isend(&(Qin[i][0]),NumPointsQ[i],MPI_DOUBLE,i,100,MPI_COMM_WORLD,&mpireq[i]);
                  MPI_Isend(&(Cin[i][0]),NumPointsC[i],MPI_DOUBLE,i,200,MPI_COMM_WORLD,&mpireq[i]);
              }
         }
         for (i=0;i<NoP;i++){            //Prwta kanw receive to megethos tou buffer
                  if (i!=selfPID){
                     MPI_Irecv(&recNumQ[i],1,MPI_INT,i,50,MPI_COMM_WORLD,&mpireq[reqcnt++]);
                     MPI_Irecv(&recNumC[i],1,MPI_INT,i,150,MPI_COMM_WORLD,&mpireq[reqcnt++]);
                  }
         }
         
         MPI_Waitall(2*(NoP-1),mpireq,mpistat); //Perimenw na mou erthoun ola ta mhkh twn buffers
         
         int psum[2]={0,0};              //Krataw to synoliko mhkos olwn twn buffer gia ta Q kai C
         for (i=0;i<NoP;i++){
             if (i!=selfPID){
                 psum[0]=psum[0]+recNumQ[i];
                 psum[1]=psum[1]+recNumC[i];
                
             }
         }
         reqcnt=0;
         for (i=0;i<NoP;i++){
                  if (i!=selfPID){        
                     recQ[i]=(double *)malloc(recNumQ[i]*sizeof(double));    //Kanw prwta allocate gnorizontas to megethos tou buffer
                     recC[i]=(double *)malloc(recNumC[i]*sizeof(double));
                     MPI_Irecv(&(recQ[i][0]),recNumQ[i]*4,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&mpireq[reqcnt++]);  //Kanw receive ton buffer me ta shmeia apo kathe process
                     MPI_Irecv(&(recC[i][0]),recNumC[i]*4,MPI_DOUBLE,i,200,MPI_COMM_WORLD,&mpireq[reqcnt++]);
                 }
        }   
        
        MPI_Waitall(2*(NoP-1),mpireq,mpistat);    //Perimenw na erthoun ola ta shmeia
        
        Q=(struct Point *)malloc(((NumPointsQ[selfPID]+psum[0])/4)*sizeof(struct Point));   //Ftiaxnw ton pinaka me ta shmeia mou Q kai C
        C=(struct Point *)malloc(((NumPointsC[selfPID]+psum[1])/4)*sizeof(struct Point));   //Synoliko megethos = Dika mou shmeia (NumPoints[selfPID]/4) + shmeia pou eginan receive (psum/4)
        for (qCnt=0;qCnt<NumPointsQ[selfPID]/4;qCnt++){  //Prwta pernaw ta dika mou shmeia
             Q[qCnt].x=Qin[selfPID][4*qCnt];             //Gnwrizw pws sth 1h thesi yparxei to x, meta to y, z, id
             Q[qCnt].y=Qin[selfPID][4*qCnt+1];           //O pinakas mou einai apo struct Point pleon     
             Q[qCnt].z=Qin[selfPID][4*qCnt+2];
             Q[qCnt].genid=(int) Qin[selfPID][4*qCnt+3];
             Q[qCnt].id=qCnt;                            //San eidiko id vazw exw to index tou shmeiou ston pinaka
        }
        
        for (cCnt=0;cCnt<NumPointsC[selfPID]/4;cCnt++){  //Idio gia ton pinaka C
             C[cCnt].x=Cin[selfPID][4*cCnt];
             C[cCnt].y=Cin[selfPID][4*cCnt+1];
             C[cCnt].z=Cin[selfPID][4*cCnt+2];
             C[cCnt].genid=(int) Cin[selfPID][4*cCnt+3];
             C[cCnt].id=cCnt;
        }
        
        for (i=0;i<NoP;i++){                        //Twra pernaw ta shmeia pou elava apo kathe diergasia
             if (i!=selfPID){                       
                 for (j=0;j<recNumQ[i]/4;j++){
                    Q[qCnt].x=recQ[i][4*j];
                    Q[qCnt].y=recQ[i][4*j+1];
                    Q[qCnt].z=recQ[i][4*j+2];
                    Q[qCnt].genid=(int) recQ[i][4*j+3];
                    Q[qCnt].id=qCnt++;
                 }
                 for (j=0;j<recNumC[i]/4;j++){
                    C[cCnt].x=recC[i][4*j];
                    C[cCnt].y=recC[i][4*j+1];
                    C[cCnt].z=recC[i][4*j+2];
                    C[cCnt].genid=(int) recC[i][4*j+3];
                    C[cCnt].id=cCnt++;
                 }
                 
              }
        }
        checkAndCorrect();                      //Elegxw an eixa alloiwsh kata tin lipsi

       	BoxNum=pow(2,BoxNum);                   //Synolikos arithmos koutiwn = 2 eis tin BoxNum         
       	BoxGridDimensions(BoxNum);              //Ypologizw tis diastaseis tou Grid twn boxes
	Boxes=(struct Box *)malloc((BoxNum/NoP) * sizeof(struct Box));  //Kanw allocate BoxNum/NoP boxes (gia kathe process)
	BoxCreate();                            //Dimiourgw ta Boxes mou
        BoxAssign();                            //Anathetw ta shmeia mou sta Boxes
        findPNeighbors(processX,processY,processZ);        //Vriskw ta geitonika process
        for (i=0;i<26;i++){              //Midenizw tous counters me ta shmeia C pou tha steilw kai tha lavw
            cSendCnt[i]=0;
            cRecCnt[i]=0;
        }
        findSendPoints();                                  //Vriskw ta shmeia C pou tha steilw sta geitonika process kai ta topothetw ston buffer sendC
        int sendTo;       //Prosorini metavliti pou deixnei ton paralipti
        for (i=0;i<neighbors;i++){                  //Stelnw shmeia C se olous tous geitones mou
             sendTo=neighborID[i];       //to PID tou paralipti
             MPI_Isend(&(cSendCnt[neighFlagID[i]]),1,MPI_INT,sendTo,50,MPI_COMM_WORLD, &mpireq[i]);    //Stelnw to megethos tou buffer kai ton buffer
             MPI_Isend(&(sendC[neighFlagID[i]][0]),cSendCnt[neighFlagID[i]],MPI_DOUBLE,sendTo,100,MPI_COMM_WORLD, &mpireq[i]);
        }
        reqcnt=0;   
        int receiveFrom;  //Prosorini metavliti pou deixnei ton apostolea
        for (i=0;i<neighbors;i++){
              receiveFrom=neighborID[i]; //PID apostolea
              MPI_Irecv(&(cRecCnt[neighFlagID[i]]),1,MPI_INT,receiveFrom,50,MPI_COMM_WORLD, &mpireq[reqcnt++]);   //Apo kathe geitona pairnw to megethos tou buffer
         }
            
        MPI_Waitall(neighbors,mpireq,mpistat);   //Perimenw ola ta megethi twn buffer pou tha lavw
            
        int sum2=0;  //Edw tha kratisw to synoliko megethos twn buffer pou tha lavw
        recC2=(double **)malloc(26*sizeof(double *)); //Kanw allocate ton buffer lipsis gia kathe geitona
        reqcnt=0;
        for (i=0;i<neighbors;i++){
               receiveFrom=neighborID[i];  //PID apostolea
               sum2=sum2+cRecCnt[neighFlagID[i]];     
               recC2[neighFlagID[i]]=(double *)malloc(cRecCnt[neighFlagID[i]]*sizeof(double));   //Kanw allocate gnwrizotas to megethos tou buffer
               MPI_Irecv(&(recC2[neighFlagID[i]][0]),cRecCnt[neighFlagID[i]],MPI_DOUBLE,receiveFrom,100,MPI_COMM_WORLD, &mpireq[reqcnt++]);    //Kanw receive ton buffer
                
        }
        MPI_Waitall(neighbors,mpireq,mpistat);    //Perimenw ola ta geitonika shmeia C

        neighC=(struct Point *)malloc(sum2/4*sizeof(struct Point));  //Kanw allocate ton pinaka neighC pou periexei ta geitonika shmeia C pou tha xreiastw
        for (i=0;i<neighbors;i++){
                for (j=0;j<cRecCnt[neighFlagID[i]]/4;j++){
                    neighC[cCnt2].x=recC2[neighFlagID[i]][4*j];          //Pernaw se struct Point ola ta shmeia pou elava (x,y,z,id)
                    neighC[cCnt2].y=recC2[neighFlagID[i]][4*j+1];
                    neighC[cCnt2].z=recC2[neighFlagID[i]][4*j+2];
                    neighC[cCnt2].genid=(int) recC2[neighFlagID[i]][4*j+3];
                    neighC[cCnt2].id=cCnt2;
                    neighC[cCnt2].boxid=BoxIdentify(neighC[cCnt2]);  //To general id tou box tou geitonikou process pou anhkei to shmeio
                    cCnt2++;
    
                 }
        }
        checkAndCorrect2();             //Elegxw tyxousa alloiwsh twn shmeiwn
        BoxNeighCreate(cCnt2);          //Dhmiourgw ta geitonika Boxes mou
        for (i=0;i<BoxNum/NoP;i++){     //Gia kathe Box, trexw thn synarthsh eureshw tou kontinoterou C gia kathe Q tou Box 
              findNearest(i);
        }
        MPI_Barrier(MPI_COMM_WORLD);    //Mesw tou Barrier metrw swsta ton xrono termatismou tou algorithmou
        t2=MPI_Wtime();
        if (selfPID==0) printf(" Elapsed  time = %f sec\n",t2-t1);   //Ektypwnw th diarkeia
        MPI_Finalize();                 //Termatismos twn processes
	return 0;
}

void checkAndCorrect(){
         int i;
         double stepx=(double) 1/pn;  //Apostasi metaksy geitonikwn processes
         double stepy=(double) 1/pm;
         double stepz=(double) 1/pk;
         for (i=0;i<qCnt;i++){
              if (Q[i].x==pstartX+stepx) Q[i].x=Q[i].x-0.000001;      //An logw alloiosis sti lipsi to shmeio anhkei (oriaka) se allo process afairw 0.000001 kai apokathistw th zhmia 
              if (Q[i].y==pstartY+stepy) Q[i].y=Q[i].y-0.000001;
              if (Q[i].z==pstartZ+stepz) Q[i].z=Q[i].z-0.000001;
         }
         for (i=0;i<cCnt;i++){
              if (C[i].x==pstartX+stepx) C[i].x=C[i].x-0.000001;
              if (C[i].y==pstartY+stepy) C[i].y=C[i].y-0.000001;
              if (C[i].z==pstartZ+stepz) C[i].z=C[i].z-0.000001;
         }
         
}

void checkAndCorrect2(){      //Idia ylopoihsh opws sthn prwti lhpsh aplws elegxw ta 2 akra kai analogws prosthafairw
         int i;
         double stepx=(double) 1/pn;
         double stepy=(double) 1/pm;
         double stepz=(double) 1/pk;
         for (i=0;i<cCnt2;i++){
              if (neighC[i].x==pstartX) neighC[i].x=neighC[i].x-0.000001;
              if (neighC[i].x==pstartX+stepx) neighC[i].x=neighC[i].x+0.000001;
              if (neighC[i].y==pstartY) neighC[i].y=neighC[i].y-0.000001;
              if (neighC[i].y==pstartY+stepy) neighC[i].y=neighC[i].y+0.000001;
              if (neighC[i].z==pstartZ) neighC[i].z=neighC[i].z-0.000001;
              if (neighC[i].z==pstartZ+stepz) neighC[i].z=neighC[i].z+0.000001;
         }
}

void findPNeighbors(int x,int y,int z){    //Ypologizw tous geitones gia to process me syntetagmenes (x,y,z) sto Grid twn process
        int i,j,o;
        int cnt=0;
        for (i=0;i<26;i++){
             neighFlag[i]=0;      //An Flag=1 exw geitona se auth thn kateuthinsi
        }
        neighbors=0;
        neighborID=(int *)malloc(7*sizeof(int));  //Allocate gia 7 geitones 
        neighFlagID=(int *)malloc(7*sizeof(int));
        int boundaries[6]={x!=0,x!=pn-1,y!=0,y!=pm-1,z!=0,z!=pk-1};  //Elegxos gia ta boundaries ths process kata +/-x,+/-y,+/-z
        int idDiff[6]={-1,+1,-pn,+pn,-(pn*pm),+(pn*pm)};             //Diafora id gia kata periptwsi pou den exw boundary kata +/-x,+/-y,+/-z
        for (i=0;i<6;i++){                                           //Elegxw prwta tous varikous geitones stis epifaneies ths process
            if (boundaries[i]){
                neighborID[neighbors]=selfPID+idDiff[i];         //An exw geitona, prosthetw to id tou, to flag tou, kai auksanw tous geitones
                neighFlagID[neighbors]=cnt;
                neighbors++;
                neighFlag[cnt]=1;
            }
            cnt++;
        }
        for (i=0;i<2;i++){                                           //Kanw to idio gia tous geitones sto idio epipedo z, me tous opoious exw koines akmes (kata x-y)
           for (j=2;j<4;j++){
               if (boundaries[i] && boundaries[j]){
                  if (neighbors>=7){
                      neighborID=(int *)realloc(neighborID,(neighbors+1)*sizeof(int));    //An kseperasw to arxiko allocate, kanw reallocate
                      neighFlagID=(int *)realloc(neighFlagID,(neighbors+1)*sizeof(int));
                  }
                  neighborID[neighbors]=selfPID+idDiff[i]+idDiff[j];
                  neighFlagID[neighbors]=cnt;
                  neighbors++; 
                  neighFlag[cnt]=1;
               }
               cnt++;
           }
        }
        for (i=0;i<2;i++){                                           //Kanw to idio gia tous geitones me tous opoious exw koines akmes (kata x-z)
           for (j=4;j<6;j++){
               if (boundaries[i] && boundaries[j]){
                  if (neighbors>=7){
                      neighborID=(int *)realloc(neighborID,(neighbors+1)*sizeof(int));    
                      neighFlagID=(int *)realloc(neighFlagID,(neighbors+1)*sizeof(int));
                  }
                  neighborID[neighbors]=selfPID+idDiff[i]+idDiff[j];
                  neighFlagID[neighbors]=cnt;
                  neighbors++; 
                  neighFlag[cnt]=1;
               }
               cnt++;
           }
        }
        for (i=2;i<4;i++){                                           //Kanw to idio gia tous geitones me tous opoious exw koines akmes (kata y-z)
           for (j=4;j<6;j++){
               if (boundaries[i] && boundaries[j]){
                  if (neighbors>=7){
                      neighborID=(int *)realloc(neighborID,(neighbors+1)*sizeof(int));    
                      neighFlagID=(int *)realloc(neighFlagID,(neighbors+1)*sizeof(int));
                  }
                  neighborID[neighbors]=selfPID+idDiff[i]+idDiff[j];
                  neighFlagID[neighbors]=cnt;
                  neighbors++; 
                  neighFlag[cnt]=1;
               }
               cnt++;
           }
        }
        for (i=0;i<2;i++){                                            //Telos, elegxw tous geitones sta epipeda z+1,z-1, me tous opoious exw koines gwnies
            for (j=2;j<4;j++){
                for (o=4;o<6;o++){
                   if (boundaries[i] && boundaries[j] && boundaries[o]){
                      if (neighbors>=7){
                          neighborID=(int *)realloc(neighborID,(neighbors+1)*sizeof(int));    
                          neighFlagID=(int *)realloc(neighFlagID,(neighbors+1)*sizeof(int));
                      }
                      neighborID[neighbors]=selfPID+idDiff[i]+idDiff[j]+idDiff[o];
                      neighFlagID[neighbors]=cnt;
                      neighbors++;
                      neighFlag[cnt]=1;
                   }
                   cnt++;
                }
            }
        }
}

void findSendPoints(){                        //Etoimazei ton buffer sendC gia ta shmeia C pou tha steilei to kathe process tous geitones tou 
     int i,i2,j2,o2,cnt;                                   
     sendC=(double **)malloc(26*sizeof(double *));         //Kanw allocate ton pinaka
     for (i=0;i<neighbors;i++){
            sendC[neighFlagID[i]]=(double *)malloc(4*cCnt*sizeof(double));
     }
     double stepx=(double) 1/n,stepy=(double) 1/m,stepz=(double) 1/k;              //apostaseis metaksy boxes kata aksones
     double stepPX=(double) 1/pn,stepPY=(double) 1/pm,stepPZ=(double) 1/pk;        //apostaseis metaksy processes kata aksones
     int limitPosition[6];
     for (i=0;i<cCnt;i++){
       limitPosition[0]=C[i].x>=pstartX && C[i].x<pstartX+stepx;
       limitPosition[1]=C[i].x>=pstartX+stepPX-stepx && C[i].x<pstartX+stepPX;
       limitPosition[2]=C[i].y>=pstartY && C[i].y<pstartY+stepy;
       limitPosition[3]=C[i].y>=pstartY+stepPY-stepy && C[i].y<pstartY+stepPY;
       limitPosition[4]=C[i].z>=pstartZ && C[i].z<pstartZ+stepz;
       limitPosition[5]=C[i].z>=pstartZ+stepPZ-stepz && C[i].z<pstartZ+stepPZ;
       cnt=0;
       for (i2=0;i2<6;i2++){                                               //Gia ola ta C shmeia elegxw an einai "oriaka", prepei dhladh na ta steilw se geitoniko process
            if (limitPosition[i2]){                                        //Gia na to kanw auto, elegxw tis syntetagmenes tou, kai an vrisketai se "oriako" box to stelnw 
               if (neighFlag[cnt]==1){                                     //sto antistoixo geitoniko process
                  sendC[cnt][cSendCnt[cnt]++]=C[i].x;                      //Stelnw tis syntetagmenes x,y,z kai to genid tou shmeiou
                  sendC[cnt][cSendCnt[cnt]++]=C[i].y;                      //Elegxw prwta ta shmeia pou einai se oriakh epifaneia tou box mou
                  sendC[cnt][cSendCnt[cnt]++]=C[i].z;
                  sendC[cnt][cSendCnt[cnt]++]=(double) C[i].genid;
               }
            }
            cnt++;
       }
       for (i2=0;i2<2;i2++){
            for (j2=2;j2<4;j2++){                                            //Sth synexeia anazhtw ta shmeia pou vriskontai sthw akmes tou box mou (x & y)
               if (limitPosition[i2] && limitPosition[j2]){
                  if (neighFlag[cnt]==1){                                     
                     sendC[cnt][cSendCnt[cnt]++]=C[i].x;                      
                     sendC[cnt][cSendCnt[cnt]++]=C[i].y;
                     sendC[cnt][cSendCnt[cnt]++]=C[i].z;
                     sendC[cnt][cSendCnt[cnt]++]=(double) C[i].genid;
                   }
                }
                cnt++;
             }
       }
        for (i2=0;i2<2;i2++){
            for (j2=4;j2<6;j2++){                                          //Sth synexeia anazhtw ta shmeia pou vriskontai sthw akmes tou box mou (x & z)
               if (limitPosition[i2] && limitPosition[j2]){
                  if (neighFlag[cnt]==1){                                     
                     sendC[cnt][cSendCnt[cnt]++]=C[i].x;                      
                     sendC[cnt][cSendCnt[cnt]++]=C[i].y;
                     sendC[cnt][cSendCnt[cnt]++]=C[i].z;
                     sendC[cnt][cSendCnt[cnt]++]=(double) C[i].genid;
                   }
                }
                cnt++;
             }
        }
        for (i2=2;i2<4;i2++){                                             //Sth synexeia anazhtw ta shmeia pou vriskontai sthw akmes tou box mou (y & z)
            for (j2=4;j2<6;j2++){
               if (limitPosition[i2] && limitPosition[j2]){
                  if (neighFlag[cnt]==1){                                     
                     sendC[cnt][cSendCnt[cnt]++]=C[i].x;                      
                     sendC[cnt][cSendCnt[cnt]++]=C[i].y;
                     sendC[cnt][cSendCnt[cnt]++]=C[i].z;
                     sendC[cnt][cSendCnt[cnt]++]=(double) C[i].genid;
                   }
                }
                cnt++;
             }
        }
        for (i2=0;i2<2;i2++){                                              //Telos anazhtw ta shmeia pou vriskontai stis 8 gonies tou box mou (x & y & z)
            for (j2=2;j2<4;j2++){
               for (o2=4;o2<6;o2++){
                  if (limitPosition[i2] && limitPosition[j2] && limitPosition[o2]){
                     if (neighFlag[cnt]==1){                              
                         sendC[cnt][cSendCnt[cnt]++]=C[i].x;                      
                         sendC[cnt][cSendCnt[cnt]++]=C[i].y;
                         sendC[cnt][cSendCnt[cnt]++]=C[i].z;
                         sendC[cnt][cSendCnt[cnt]++]=(double) C[i].genid;                           
                      }
                   }
                   cnt++;
                }
           }
        }                   
      }
}                     
             
void pGridDimensions(int num) {       
        pn=1;                       //Ypogizw tis diastaseis tou Grid twn process
        pm=1;                       //To pn antistoixei sta process pou exw kata x, to pm kata y kai to pk kata z
        pk=1;
        int i=0;
        while (pn*pm*pk!=num){
              switch (i%3){
                 case 0:
                   pn=2*pn;
                   break;
                 case 1:
                   pm=2*pm;
                   break;
                 case 2:
                   pk=2*pk;
                   break;
              }
        i++;
        }
}
       
void BoxGridDimensions(int num){   //Ypologizw tis diastaseis tou Grid twn Boxes
	n=1;                       // n = Boxes kata x, m = Boxes kata y, k = Boxes kata z
        m=1;
        k=1;
        int i=0;
        while (n*m*k!=num){
          switch (i%3){
                case 0:
                   k=k*2;
                   break;
                case 1:
                   m=m*2;
                   break;
                case 2:
                   n=n*2;
                   break;
              }
         i++;
         }
n_2=n/pn;                        //n_2 = Boxes / process kata x, m_2 = Boxes / process kata y, k_2 = Boxes / process kata z
m_2=m/pm;
k_2=k/pk;
}

void BoxAssign(){                //Anathetw ola ta shmeia Q kai C sto Box pou anhkoun
        int i=0;
        for (i=0;i<qCnt;i++){          
            Q[i].boxid=gen2spec(BoxIdentify(Q[i]));         //Mesw ths BoxIdentify vriskw to genid tou Box. To metatrepw se special gia na exw amesh prosvash ston pinaka Boxes
            if (Boxes[Q[i].boxid].qpoints>=qCnt/(BoxNum/NoP)){   //An kseperasa to arxiko allocate prosthetw mia thesh Box
                 Boxes[Q[i].boxid].qpointsID=(int *)realloc(Boxes[Q[i].boxid].qpointsID,(Boxes[Q[i].boxid].qpoints+1)*sizeof(int));
            }
            Boxes[Q[i].boxid].qpointsID[Boxes[Q[i].boxid].qpoints]=Q[i].id;  //prosthetw to id tou shmeiou sta qID tou Box
            Boxes[Q[i].boxid].qpoints++;                                     //Auksanw ta qPoints tou Box
        }
        for (i=0;i<cCnt;i++){       //Kanw antistoixh douleia gia ta shmeia C
            C[i].boxid=gen2spec(BoxIdentify(C[i])); 
            if (Boxes[C[i].boxid].cpoints>=cCnt){
                 Boxes[C[i].boxid].cpointsID=(int *)realloc(Boxes[C[i].boxid].cpointsID,(Boxes[C[i].boxid].cpoints+1)*sizeof(int));
            }
            Boxes[C[i].boxid].cpointsID[Boxes[C[i].boxid].cpoints]=C[i].id;
            Boxes[C[i].boxid].cpoints++;
        }
}

int ProcIdentify(double p[3]){     //Vriskw to process sto opoio anhkei to shmeio p[3]
        double stepx=(double) 1/pn;  
        double stepy=(double) 1/pm;
        double stepz=(double) 1/pk;
        int x=p[0]/stepx;  //Ypologizw ta x,y,z tou process
        int y=p[1]/stepy;
        int z=p[2]/stepz;
        int process=x+y*pn+z*(pn*pm);  //Afou ypologisw ta x,y,z tou process to id tou einai auto
        return process;
}

int BoxIdentify(struct Point p){     //Ypologizw to generalID tou Box sto opoio anhkei to shmeio p
        float stepx=(float) 1/n;
        float stepy=(float) 1/m;
        float stepz=(float) 1/k;
        int x,y,z;
        x=p.x/stepx;
        y=p.y/stepy;
        z=p.z/stepz;
        int id=x+y*n+z*(n*m);  //Opos stin periptwsi tou process, me ton idio tropo prokyptei to genID tou Box
        return id;
}
int gen2spec(int gen){   //Metatropi general ID se special gia Boxes
        int idStartX=processX*n_2;  //ta idStartX,Y,Z einai san threshold gia ton ypologismo tou ID
        int idStartY=processY*m_2*n;
        int idStartZ=processZ*k_2*n*m;
        int spec=gen-idStartX-idStartY-idStartZ;
        int x=spec%n;  
        int y=(spec/n)%m;
        int z=spec/(n*m);
        spec=x+y*n_2+z*(n_2*m_2);   //Meta apo tous ypologismous prokyptei to special ID
        return spec;
        
}

int spec2gen(int spec){   //Antistrofi diadikasia me prin
        int idStartX=processX*n_2;
        int idStartY=processY*m_2*n;
        int idStartZ=processZ*k_2*n*m;
        int x=spec%n_2;
        int y=(spec/n_2)%m_2;
        int z=spec/(n_2*m_2);
        int gen=x+y*n+z*(n*m)+idStartX+idStartY+idStartZ;  //Ypologizw to general ID
        return gen;
}

int ifMine(int boxid){                                //Anagnorizo an einai diko mou ena kouti
         int i,flag=0;
         int spec=gen2spec(boxid);                    //To metatrepo se special id vlepw an einai ta oria twn special id kathe process
         if (spec>=0 && spec<BoxNum/NoP) flag=1;      //To ksanakanw gen kai vlepw to flag. Flag=1 exei mono h process sthn opoia anhkei
         if (spec2gen(spec)!=boxid) flag=0;
         return flag;
}
     
void BoxCreate(){                                       //Synartisi dimiourgias twn koutiwn
        double stepx=(double) 1/n;                      //Vimata (apostaseis) kata x, y, z twn koutiwn
	double stepy=(double) 1/m;
	double stepz=(double) 1/k;
	int i,j,o,i2,j2,o2;                             //General Counters
	int ids=0;                                      //counter twn koutiwn
        int neigh;                                      //arithmos geitonwn kathe koutiou
        int boundaries[6];                              //Pinakas bool me ta boundaries. Elegxei an yparxoun boundaries kata +/-x,+/-y,+/-z. An yparxei boundary = 0
        int idDiff[6]={-1,+1,-n,+n,-(n*m),+(n*m)};      //Pinakas me ta id Differencies gia kathe periptwsi boundary
        for (o=0;o<k_2;o++){
	     for (j=0;j<m_2;j++){
                  for (i=0;i<n_2;i++){
			Boxes[ids].xmin=pstartX+i*stepx;  //Anathetw se kathe kouti to xmin, ymin kai zmin, kathos kai to eidiko kai geniko id
                        Boxes[ids].ymin=pstartY+j*stepy;
			Boxes[ids].zmin=pstartZ+o*stepz;
			Boxes[ids].specid=ids;
                        Boxes[ids].genid=spec2gen(ids);
                        
                        Boxes[ids].qpoints=0;                  //Arxikopoiw me 0 ta Q kai C points tou box, kai kanw allocate ton pinaka me ta ID twn Q kai C
                        Boxes[ids].cpoints=0;
                        Boxes[ids].qpointsID=(int *)malloc(qCnt*sizeof(int));
                        Boxes[ids].cpointsID=(int *)malloc(cCnt*sizeof(int));
                        Boxes[ids].neighbors=0;                //Arxikopoiw me 0 tous geitones tou box kai kanw allocate gia tous geitones tou
                        Boxes[ids].neighborID=(int *)malloc(7*sizeof(int));
                        
                        neigh=0;
                        boundaries[0]=Boxes[ids].xmin!=0;                     //Elegxw ta boundaries tou box, an diladi einai oriako sto Grid twn Boxes
                        boundaries[1]=Boxes[ids].xmin!=1-stepx;
                        boundaries[2]=Boxes[ids].ymin!=0;
                        boundaries[3]=Boxes[ids].ymin!=1-stepy;
                        boundaries[4]=Boxes[ids].zmin!=0;
                        boundaries[5]=Boxes[ids].zmin!=1-stepz;
                        for (i2=0;i2<6;i2++){                                  //Anagnorizw prwta tous vasikous geitones, dhladh autous pou vriskontai stis epifaneies tou
                            if (boundaries[i2]){                             //Elegxw ta boundaries (kata x,y,z)
                                Boxes[ids].neighbors++;
                                Boxes[ids].neighborID[neigh]=spec2gen(ids)+idDiff[i2];
                                neigh++;
                            }
                        }
                        for (i2=0;i2<2;i2++){                                //Sth sunexeia anagnorizw tous geitones pou vriskontai stis akmes tou sto idio epipedo (x,y)
                            for (j2=2;j2<4;j2++){
                                if (boundaries[i2] && boundaries[j2]){       //Elegxw ta boundaries (x & y)
                                    if (neigh>=7){
                                        Boxes[ids].neighborID=(int*)realloc(Boxes[ids].neighborID,(neigh+1)*sizeof(int));   //An kseperasw to arxiko allocate, kano reallocate +1 thesi
                                    }
                                    Boxes[ids].neighbors++;
                                    Boxes[ids].neighborID[neigh]=spec2gen(ids)+idDiff[i2]+idDiff[j2];
                                    neigh++;
                                }
                            }
                        }
                        for (i2=0;i2<2;i2++){                                //Sth sunexeia anagnorizw tous geitones pou vriskontai stis akmes tou sto idio epipedo (x,z)
                            for (j2=4;j2<6;j2++){
                                if (boundaries[i2] && boundaries[j2]){       //Elegxw ta boundaries (x & z)
                                    if (neigh>=7){
                                        Boxes[ids].neighborID=(int*)realloc(Boxes[ids].neighborID,(neigh+1)*sizeof(int));   
                                    }
                                    Boxes[ids].neighbors++;
                                    Boxes[ids].neighborID[neigh]=spec2gen(ids)+idDiff[i2]+idDiff[j2];
                                    neigh++;
                                }
                            }
                        } 
                        for (i2=2;i2<4;i2++){                                //Sth sunexeia anagnorizw tous geitones pou vriskontai stis akmes tou sto idio epipedo (x,y)
                            for (j2=4;j2<6;j2++){
                                if (boundaries[i2] && boundaries[j2]){       //Elegxw ta boundaries (y & z)
                                    if (neigh>=7){
                                        Boxes[ids].neighborID=(int*)realloc(Boxes[ids].neighborID,(neigh+1)*sizeof(int));  
                                    }                                   
                                    Boxes[ids].neighbors++;
                                    Boxes[ids].neighborID[neigh]=spec2gen(ids)+idDiff[i2]+idDiff[j2];
                                    neigh++;
                                }
                            }
                        }   
                        for (i2=0;i2<2;i2++){                                //Telos elegxw tous geitones pou vriskontai sto pio panw kai pio katw epipedo z, akmes kai gonies
                            for (j2=2;j2<4;j2++){
                                for (o2=4;o2<6;o2++){
                                    if (boundaries[i2] && boundaries[j2] && boundaries[o2]){
                                       if (neigh>=7){                           //Elegxw ta boundaries (x & y & z)
                                          Boxes[ids].neighborID=(int*)realloc(Boxes[ids].neighborID,(neigh+1)*sizeof(int));   
                                       }                                     
                                       Boxes[ids].neighbors++;
                                       Boxes[ids].neighborID[neigh]=spec2gen(ids)+idDiff[i2]+idDiff[j2]+idDiff[o2];
                                       neigh++;
                                    }
                                }
                            }
                       }
                       ids++;       //Gia kathe box auksano kata 1 ton counter twn boxes
                  }
              }
	  }
}

void BoxNeighCreate(int cpoints){                            //Dhmiourgei ton pinaka twn geitonikwn boxes, kai anathetei ta geitonika C pou elava stis theseis tou pinaka
           int i,j,box,size=2*(n_2*m_2+n_2*n_2+k_2*m_2+k_2)+4*(n_2+m_2+k_2)+8;     //size = arithmos geitonikwn koutiwn (einai o megistos arithmos geitonikwn boxes pou mporei na exei kathe box)
           int boxesId[size],cpointsCnt[size],flag;
           neighBoxes=(struct Box *)malloc(size*sizeof(struct Box));      //Kanw allocate ton pinaka
           for (i=0;i<size;i++){                                          //Arxikopoiw me 0 tous counters kai me -1 ton pinaka me ta id
                boxesId[i]=-1;
                cpointsCnt[i]=0;
                neighBoxes[i].cpoints=0;
           }
           for (i=0;i<cpoints;i++){              
                box=BoxIdentify(neighC[i]);                               //Gia kathe C pou elava, vriskw se poio box anhkei
                flag=0;
                j=0;
                while (flag==0 && j<neighBoxesCnt) {                      //Psaxnw se poia thesh tou pinaka me ta id yparxei to sygkekrimeno box, An yparxei
                     if (boxesId[j]==box){
                         flag=1;
                     }
                     j++;
                }
                if (flag==0){                                            //An den yparxei, dhmiourgw to Box, kanw allocate gia ta shmeia C tou, pernaw to id tou kai auksanw kata 1 ta Cpoints tou
                     boxesId[neighBoxesCnt]=box;
                     neighBoxes[neighBoxesCnt].genid=box;
                     neighBoxes[neighBoxesCnt].cpointsID=(int *)malloc(cpoints*sizeof(int));
                     neighBoxes[neighBoxesCnt].cpointsID[cpointsCnt[neighBoxesCnt]++]=neighC[i].id;
                     neighBoxes[neighBoxesCnt++].cpoints++;
                }
                else if (flag==1){
                     neighBoxes[j-1].cpointsID[cpointsCnt[j-1]++]=neighC[i].id;   //An yparxei, aplws pigainw sth thesi auth kai prosthetw to shmeio C
                     neighBoxes[j-1].cpoints++;
                }  
           }
}

int findIndexNeigh(int box){                          //Vriskei se poia thesi tou pinaka neighBoxes vrisketai to box
       int i=0,flag=0;                                
       while (flag==0 && i<neighBoxesCnt){
           if (neighBoxes[i].genid==box) {
              flag=1;
           }
       i++;
       }
       return i-1;
}

double distance(struct Point p1,struct Point p2){                    //Vriskei thn apostasi 2 Points ston 3-D xwro
           double dist=sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
           return dist;
}

void findNearest(int boxid){                               //Gia kathe Q tou boxid Box vriskoume ton kontinotero geitona tou synolou C
       double minDist,dist;                                //O kontinoteros geitonas mporei na vrisketai sto idio kouti, h sta geitonika
       int i,j,o,nearestID=-1,belong;
       int index;
       for (i=0;i<Boxes[boxid].qpoints;i++){              
           minDist=sqrt(3);                                
           for (j=0;j<Boxes[boxid].cpoints;j++){
               dist=distance(Q[Boxes[boxid].qpointsID[i]],C[Boxes[boxid].cpointsID[j]]); //Prwta kanw search sta C tou idiou koutiou, vriskw thn eukleidia apostasi kai th sygkrinw me thn minimum Distance
               if (dist<minDist){                                                        //Ean einai mikroteri, exw nea MinDist
                    minDist=dist;
                    nearestID=Boxes[boxid].cpointsID[j];
               }
           }
           Q[Boxes[boxid].qpointsID[i]].minDist=minDist;
           Q[Boxes[boxid].qpointsID[i]].nearestid=nearestID;
        }
        for (o=0;o<Boxes[boxid].neighbors;o++){                     //Sth synexeia kanw to idio gia ta geitonika boxes
              belong=ifMine(Boxes[boxid].neighborID[o]);            //Vriskw prwta an to geitoniko Box einai diko mou, an dhladh prepei na paw ston pinaka Boxes h neighBoxes
              if (belong==1){                                       //An belong = 1 anhkei se mena to box
                 for (i=0;i<Boxes[boxid].qpoints;i++){
                   minDist=Q[Boxes[boxid].qpointsID[i]].minDist;      //Anaktw thn timh minDist kai nearestID 
                   nearestID=Q[Boxes[boxid].qpointsID[i]].nearestid;
                   for (j=0;j<Boxes[gen2spec(Boxes[boxid].neighborID[o])].cpoints;j++){       //Psaxnw gia ola ta shmeia C tou geitonikou box
                       dist=distance(Q[Boxes[boxid].qpointsID[i]],C[Boxes[gen2spec(Boxes[boxid].neighborID[o])].cpointsID[j]]);
                       if (dist<minDist){
                           minDist=dist;
                           nearestID=Boxes[gen2spec(Boxes[boxid].neighborID[o])].cpointsID[j];
                       }
                   }
                   Q[Boxes[boxid].qpointsID[i]].minDist=minDist;
                   Q[Boxes[boxid].qpointsID[i]].nearestid=nearestID;
                 }
              }
              if (belong==0){                                         //belong = 0 den anhkei se mena
                   index=findIndexNeigh(Boxes[boxid].neighborID[o]);  //Vriskw se poia thesh tou neighBoxes vrisketai to box
                   for (i=0;i<Boxes[boxid].qpoints;i++){
                      minDist=Q[Boxes[boxid].qpointsID[i]].minDist;
                      nearestID=Q[Boxes[boxid].qpointsID[i]].nearestid;
                      for (j=0;j<neighBoxes[index].cpoints;j++){      //Psaxnw ksana ola ta shmeia C tou box kai vriskw to kontinotero
                           dist=distance(Q[Boxes[boxid].qpointsID[i]],neighC[neighBoxes[index].cpointsID[j]]);
                           if (dist<minDist){
                                minDist=dist;
                                nearestID=neighBoxes[index].cpointsID[j];
                            }
                      }
                      Q[Boxes[boxid].qpointsID[i]].minDist=minDist;            //Pernaw thn teliki timi tou minDist kai tou nearestID    
                      Q[Boxes[boxid].qpointsID[i]].nearestid=nearestID;
                   }
              }
        }
}



