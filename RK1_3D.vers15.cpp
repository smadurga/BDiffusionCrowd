/// Reaction Kinetics in a macromolecular
// overcrowded environment

// g++
// 3D version
// PBC, off-latice
// Generation of new positions as random
// Quadratic potential
// 
// RMSD calculation and D linear fit implemented.
// Microscopic rate constants calculation implemented!
// Implementation in progress of thermodynamic rate constants

// Canvi formula de calcul de les microconstants


 #include <stdlib.h>  // Per compatibilitat 
 #include <time.h>    // Per fer idum depenent del temps
 #include <iostream>
 using namespace std;    // Per compatibilitat

 #define _USE_MATH_DEFINES // Pel nombre PI (M_PI)
 #include <math.h>
 #include <fstream>
 #include "mpi.h"
 #include <string.h>
 #include <stdio.h>
 #define MAXPART 100000
 #define MAXTEMPS 1000000
 #define MAXDESP 0.2

////////////////
// SUBRUTINES //
////////////////


 int main(int narg, char **arg);
 void   lectura_input(void);
 int    control_parametres(void);
 void   inicialitzacio_sistema(void);
 void   interaccions_sistema(void);
 void   reaccions (void); 
 double noves_coordenades(void); 
 void   trajectoria(void);
 void   vmd_script(void);
 void	r2_acumuladors(void);
 void	r2_averages(void);
 void	Dif_coef(void);
 void	constants_calc (void);
 void	constants_each50(void);
 void   output_variables(void);
 void   create_output(void);
 void   create_output_r2(void);
 void	kmichaelis_calc(void);
 float  ran2(long *idum); 
 void coordpbc3D(double xnopbc, double ynopbc, double znopbc, double *xsipbc, double *ysipbc, double *zsipbc, int celln);
 void difpbc3D(double xnopbc, double ynopbc, double znopbc, double *xsipbc, double *ysipbc, double *zsipbc);
 void  intercanvi_infor(void); 

//////////////
// FUNCIONS //
//////////////
 double DstokesEinstein(double radiSE);
// double PotencialQuadratic (double diffpbcx, double diffpbcy, double diffpbcz, double R_n, double R_n1, double *dVx,  double *dVy, double *dVz);
 double r2_sample(double x, double y, double z, double px, double py, double pz, int xparet, int yparet, int zparet);
// double PotencialAntic ();
 double D_RegLineal (int i_max, double ts, double r2_average []);
 double microconstantrate (double concC, double concA, double concB, double ts);
double PotencialCSW (double r_x, double r_y, double r_z, double S_h, double S_c, double *dVx,  double *dVy, double *dVz);


//////////////////////
// GLOBAL VARIABLES //
//////////////////////

 double px [MAXPART+1];
 double py [MAXPART+1];
 double pz [MAXPART+1];
 double x_ini [MAXPART+1];
 double y_ini [MAXPART+1];
 double z_ini [MAXPART+1];
 int typ [MAXPART+1];
 int parella [MAXPART+1];
 double grad_x [MAXPART+1];
 double grad_y [MAXPART+1];
 double grad_z [MAXPART+1];
 double radi_h [MAXPART+1];
 double radi_c [MAXPART+1];
 double Dif [MAXPART+1];
 int xparet[MAXPART+1];
 int yparet[MAXPART+1];
 int zparet[MAXPART+1];


 double cE [MAXTEMPS+1];
 double cS [MAXTEMPS+1];
 double cC [MAXTEMPS+1];
 double cP [MAXTEMPS+1];
 double cES [MAXTEMPS+1];
 double Cd1 [MAXTEMPS+1];
 double Cd2 [MAXTEMPS+1];
 double acum_r2_E [MAXTEMPS+1];
 double acum_r2_S [MAXTEMPS+1];
 double acum_r2_C [MAXTEMPS+1];
 double acum_r2_P [MAXTEMPS+1];
 double acum_r2_O1 [MAXTEMPS+1];
 double acum_r2_O2 [MAXTEMPS+1];
 double acum_i_E [MAXTEMPS+1];
 double acum_i_S [MAXTEMPS+1];
 double acum_i_C [MAXTEMPS+1];
 double acum_i_P [MAXTEMPS+1];
 double acum_i_O1 [MAXTEMPS+1];
 double acum_i_O2 [MAXTEMPS+1];
 double ac [MAXTEMPS+1];
 double r2_average_E [MAXTEMPS+1];
 double r2_average_S [MAXTEMPS+1];
 double r2_average_C [MAXTEMPS+1];
 double r2_average_P [MAXTEMPS+1];
 double r2_average_O1 [MAXTEMPS+1];
 double r2_average_O2 [MAXTEMPS+1];
 double k1_i [MAXTEMPS+1];
 double k_1_i [MAXTEMPS+1];
 double k2_i [MAXTEMPS+1];
 double kmichaelis_concs_i [MAXTEMPS+1];
 double kmichaelis_kts_i [MAXTEMPS+1];
 double kmeach50 [MAXTEMPS+1];
 double ctsacum [MAXTEMPS+1][50][3];
 double aveach50 [MAXTEMPS+1][3];


 const double R = 8.3144621; // J mol-¹K-¹
 const double kb = 1.3806488e-23; // m2 kg s-1 K-1
 const double visc = 0.891e-3; // kg s-1 m-1	Dynamic viscosity of water at 25 C


 int i,j,k,m,celln,celln1,c_op,freqtr,i_eq,mmax,ttt;
 double sum_x,sum_y,sum_z,tp,t,ts,a,b,c,r,r1,r2,r3,g,grandom,h,hrandom;
 double D_E,D_S,D_P,D_O1,D_O2,Dp,Beta,alpha,fr,fr_o,vx,vy,vz;
 double x1,x2,wy,d,e,p1,p2,p3,Vo,kpair,dVx,dVy,dVz;
 double rad_h_E,rad_h_S,rad_h_P,rad_h_O1,rad_h_O2;
 double rad_c_E,rad_c_S,rad_c_P,rad_c_O1,rad_c_O2;
 double T;
 double r2_i;
 double D_E_RMSD, D_S_RMSD, D_C_RMSD, D_P_RMSD, D_O1_RMSD, D_O2_RMSD;
 int sample,sample_max,i_max;
 int cell_max,count_E,count_S,count_C,count_P,count_ES; 
 int cell_new;
 int count_Cd1, count_Cd2;
 int Num_Eini,Num_Sini,Num_Oini1,Num_Oini2;
 double costat_box_ES,costat_box_Obstacles;
 double k1_av, k1, k2_av, k2, k_1_av, k_1;
 double kmichaelis_concs_av, kmichaelis_kts_av;
 double Ur, S_h, S_c;


 long idum = -1234;

 // variables paralel

 int totalcpu,nrank;
 int MASTER=0;
 MPI_Request request;
 int ndiv,res,icpu,Nini,Nfin;
 double acum_aux [MAXTEMPS+1];


 
 /////////////////////////
 //	Outputs		//
 /////////////////////////

 ofstream posOut("rk_output.txt");	// Output with number of particles vs time
 ofstream poskts("microk.txt"); 	// Output with macrokts vs time, and mean kt's
 ofstream poskts50("microk_each50.txt"); // Output with macrokts vs time averaged each 50 its.
 ofstream posRAVE("r2_output.txt");	// Output with the average RMSD of all the particles of the system
 ofstream posKM("michaelisct.txt");	// Output with K_eq vs time and mean K_eq's
 ofstream posKM50("michaelisct_each50.txt"); 
 ofstream posDiff("Diffusion_Coeff.txt");
 ofstream posDisp("Displacement.txt");
 // ofstream pot("dades_pot.txt");
					// The output file (cout) contains the input parameters and
					// the linear regression of RMSD to obtain an approximate
					// initial guess of D.



 ///////////////////////////////////
 // Variables d'entrada i unitats //
 ///////////////////////////////////
 
 // rad_i [nm]		// i = E,S,O1,O2,P
 // ts [ns]		// time step
 // V0 [nm2/ns]		// Parametre d'interaccio pel potencial no quadratic
 // D_i [nm2/ns] 	// i = E,S,O1,O2,P
 // p1			// Probabilitat de reaccio E + S -> C
 // p2			// Probabilitat de reaccio inversa C - > E + S
 // p3			// Probabilitat de reaccio C -> E + P
 // Num_Eini		// Partícules d'E inicial
 // Num_Sini		// Partícules de S inicial
 // Num_Oini1		// Partícules de O inicial tipus 1
 // Num_Oini2		// Partícules de O inicialstipus 2
 // idum		// Generació d'un nombre aleatori
 // i_max		// Nombre maxim d'iteracions 
 // freqtr [ns]		// Cada quan s'escriu la trajectoria
 // kpair [J/mol nm2]	// Constant del potencial quadratic
 

/*********************************************************************/


 int main (int narg, char **arg) 
{

 int k;


 MPI_Init(&narg,&arg);  // Inicializa MPI




 MPI_Comm_size(MPI_COMM_WORLD,&totalcpu);
 MPI_Comm_rank(MPI_COMM_WORLD,&nrank);
 if(nrank==MASTER) {
    printf("Total cpus = %d, proc ID = %d\n",totalcpu,nrank);
    }



 lectura_input();
 k=control_parametres();
 if(k==-1) {return(-2);}



 // Nomes printeja el MASTER

 if(nrank==MASTER) {
 output_variables();
	 cout << "Calculation performed at T= " << T << endl;
	 posOut<<"#it temps[ns]  Enzyme  Substrate  Complex Product  rES  Cd1  Cd2"<<endl;
	 poskts<<"#temps[ns]	k1[nm3part-1ns-1]	   k-1[ns-1]       k2[ns-1]" <<endl;
	 poskts50<<"#temps[ns]	k1[nm3part-1ns-1]	   k-1[ns-1]       k2[ns-1]" <<endl;
	 posRAVE<<"#temps[ns]    RMSD_E    RMSD_S    RMSD_C    RMSD_P    RMSD_O1    RMSD_O2" << endl;
	 posKM<<"#temps	KM(concs) [p nm-3]	KM(kts) [p nm-3]"<<endl;
         posKM50<<"#temps[ns]	KM(kts, cada50) [p nm-3]"<<endl;
         posDisp<<"#iteracio    despla�ament_r(nm)" << endl;
}

 div_t divresult_a;



// divisio del treball entre els procesadors. Si la divisio de treball no es exacte, l'ultim procesador s'encarrega del treball extra

ndiv=sample_max/totalcpu;
res=sample_max%totalcpu;

if (nrank != totalcpu-1){
Nini=nrank*ndiv+1;
Nfin=(nrank+1)*ndiv;


}
else
{
Nini=nrank*ndiv+1;
Nfin=sample_max;


}




// Inicialitzador dels acumuladors
 for (i=1; i<=i_max; i++) {acum_r2_E [i] = 0.0; acum_r2_S [i] = 0.0; acum_r2_P [i] = 0.0; acum_r2_C [i] = 0.0; acum_r2_O1 [i] = 0.0; acum_r2_O2 [i] = 0.0; k1_i[i] = 0.0;}

/*********************/
/* Inici sample loop */
/*********************/
 for (sample=Nini; sample<=Nfin; sample++) {
	for (celln=1; celln<=cell_max; celln++) {xparet[celln] = 0; yparet[celln] = 0; zparet[celln] = 0;}
 	inicialitzacio_sistema();
        for (celln=1; celln<=cell_max; celln++){

        }

/*******************/
/* Inici time loop */
/*******************/

	 for (i=1; i<=i_max; i++) {

 		interaccions_sistema();
		if ((i>i_eq)&&(p1 != 0.)) {reaccions();}	// Nomes poden reaccionar despres de la fase d'equilibratge
 		noves_coordenades();
	
		// Initial positions (for the calculation of the RMSD) //
		if(i<=(i_eq)) {
			for (celln=1; celln<=cell_max; celln++) {
			x_ini [celln] = px [celln];
			y_ini [celln] = py [celln];
			z_ini [celln] = pz [celln];
			xparet[celln] = 0;
			yparet[celln] = 0;
			zparet[celln] = 0;
			if ((D_O1==0.0)&&(typ[celln]==3)) {Dif[celln]=D_O1;}
			if ((D_O2==0.0)&&(typ[celln]==6)) {Dif[celln]=D_O2;}
		}
		}

		// Imprimeix la trajectoria pel primer sample cada freqtr.
		if((sample==1)&&(i>=(i_eq+1))) {if (i%freqtr==0) { trajectoria(); } } 
		// Calcula el MSD despres de fer el procediment d'equilibracio del sistema.
		if(i>=(i_eq+1)) {r2_acumuladors();}
 

	 } // end of time loop //


 } // end of sample loop //
 

// Els procesadors envien al MASTER la informacio sobre els promitjos de r2 que han realitzat per cada particula


intercanvi_infor();



 // el procesado de datos y el output lo realiza el MASTER










 if (nrank==MASTER) {
 r2_averages();
 Dif_coef();
 create_output_r2();
 create_output();
 vmd_script();
 if (Num_Eini!=0. && Num_Sini!=0. && p1!=0.) {constants_calc(); constants_each50(); kmichaelis_calc();}
}


 MPI_Finalize();


 return 0;
 } // end MAIN();







//////////////////////////////////////////////////////////////////////////////////
//		///// SUBROUTINES AND FUNCITONS /////				//
/////////////////////////////////////////////////////////////////////////////////


/******************************************/
void constants_calc(void)
{

 for (i=(i_eq+1); i<=i_max; i++) {
 if ((cE[i]==0)||(cS[i]==0)) {k1_i[i]=0.0;} else {k1_i[i] = microconstantrate (cES[i]/sample_max, cE[i]/sample_max, cS[i]/sample_max, ts);}
 if (cC[i]==0) {k_1_i[i]=0.0;} else {k_1_i[i] = microconstantrate (Cd1[i]/sample_max, cC[i]/sample_max, 1.0, ts);}
 if (cC[i]==0) {k2_i[i]=0.0;}  else {k2_i[i]  = microconstantrate (Cd2[i]/sample_max, cC[i]/sample_max, 1.0, ts);}
 poskts << (i-i_eq-1)*ts << "  " << k1_i[i]*pow(costat_box_ES,3) << "  " << k_1_i[i] << "  " << k2_i[i] <<  endl;

 }

 k1_av = 0.0;
 k2_av = 0.0;
 k_1_av = 0.0;

 for (i=(i_eq+1); i<=i_max; i++) {
 k1_av += k1_i[i];
 k2_av += k2_i[i];
 k_1_av += k_1_i[i];
 }

 k1 = k1_av*pow(costat_box_ES,3)/(i_max-i_eq);
 k2 = k2_av/(i_max-i_eq);
 k_1 = k_1_av/(i_max-i_eq);


	poskts << " " << endl;
	poskts << " " << endl;
	poskts <<"#+++++++++++++++++++++++++++++++++++++" << endl;
	poskts <<"#+Microscopic reaction rate constants+" << endl;
	poskts <<"#+++++++++++++++++++++++++++++++++++++" << endl;
 	poskts << "#k1 = " << k1 << " nm3 particles-1 ns-1;	"<< k1*1E-24*R*1E9/kb <<" L mol-1 s-1"  << endl;
 	poskts << "#k-1 = " << k_1 << " ns-1;	" << k_1*1E9 <<" s-1" << endl;
 	poskts << "#k2 = " << k2 << " ns-1;	" << k2*1E9 << " s-1" << endl;

}

/******************************************/
void constants_each50(void)
{
 mmax=i_max/50;
 ttt=0;
 for (k=1; k<=mmax; k++) {
	for (m=1; m<=50; m++) {
		ctsacum[k][m][1] = k1_i[m+(50*ttt)]*pow(costat_box_ES,3);
		ctsacum[k][m][2] = k_1_i[m+(50*ttt)];
		ctsacum[k][m][3] = k2_i[m+(50*ttt)];
	}
	ttt += 1;
 }


 for (k=1; k<=mmax; k++) {
	ac[1] = 0.0;
	ac[2]=0.0;
	ac[3]=0.0;
	for (m=1; m<=50; m++) {
		for (j=1; j<=3; j++) {
			ac[j] += ctsacum[k][m][j];
		}
	}
	for (j=1; j<=3; j++) {
		aveach50[k][j] = ac[j]/50.0;
	}
 }
 for (k=1; k<=mmax; k++) {
 poskts50 << ((k-1)*5) << "  " << aveach50[k][1]  << "  " << aveach50[k][2] << "  " << aveach50[k][3] <<  endl;
 if (aveach50[k][1] == 0.0) {kmeach50[k] = 0.0;} else {kmeach50[k] = (aveach50[k][2]+aveach50[k][3])/aveach50[k][1];}
 posKM50 << ((k-1)*5) << "	" <<kmeach50[k]  << endl;	//Michaelis constant calculated with k's averaged each 50 iters
 }

}

/******************************************/
 void kmichaelis_calc(void)
 {

 for (i=(i_eq+1); i<=i_max; i++) {
 if (cC[i]==0) {kmichaelis_concs_i[i]=0;} else {kmichaelis_concs_i[i]= ( (cE[i]/sample_max) * (cS[i]/sample_max) ) / ((cC[i]/sample_max)*pow(costat_box_ES,3) );}
 if (k1_i[i]!=0) {kmichaelis_kts_i[i]=(k_1_i[i]+k2_i[i])/k1_i[i];} else {kmichaelis_kts_i[i]=0.0;}
 posKM << (i-i_eq-1)*ts << "  " << kmichaelis_concs_i[i] << "	" << kmichaelis_kts_i[i]<<endl;
}

 kmichaelis_concs_av =0.0;
 kmichaelis_kts_av = 0.0;
 for (i=(i_eq+1); i<=i_max; i++) {
 	kmichaelis_concs_av += kmichaelis_concs_i[i];
	kmichaelis_kts_av += kmichaelis_kts_i[i];
 }

 kmichaelis_concs_av = kmichaelis_concs_av/(i_max-i_eq);
 kmichaelis_kts_av = kmichaelis_kts_av/(i_max-i_eq);


 posKM << endl;
 posKM << endl;
 posKM <<"#Average KM constant (as a dissociation constant, by means of [E][S]/[C])"<<endl;
 posKM <<"#KM(concs) = " << kmichaelis_concs_av << endl;
 posKM << endl;
 posKM <<"#Average KM constant (as KM definition, by means of (k-1+k2)/k1)"<<endl;
 posKM <<"#KM(kts) = " << kmichaelis_kts_av << endl;
 posKM << endl;

}

/******************************************/
 void coordpbc3D(double xnopbc,double ynopbc, double znopbc, double *xsipbc,double *ysipbc, double *zsipbc, int celln)
 {
	
  *xsipbc=xnopbc;
  *ysipbc=ynopbc;
  *zsipbc=znopbc;


     if(xnopbc<0) { *xsipbc+=costat_box_ES; xparet[celln]=xparet[celln]-1;}
     if(xnopbc>costat_box_ES) { *xsipbc-=costat_box_ES; xparet[celln]=xparet[celln]+1;}   

     if(ynopbc<0) { *ysipbc+=costat_box_ES; yparet[celln]=yparet[celln]-1;}
     if(ynopbc>costat_box_ES) { *ysipbc-=costat_box_ES; yparet[celln]=yparet[celln]+1;}   

     if(znopbc<0) { *zsipbc+=costat_box_ES; zparet[celln]=zparet[celln]-1;}
     if(znopbc>costat_box_ES) { *zsipbc-=costat_box_ES; zparet[celln]=zparet[celln]+1;}   

}

/***********************************************************/
void difpbc3D(double xnopbc,double ynopbc,double znopbc,double *xsipbc,double *ysipbc, double *zsipbc)
{
   *xsipbc=xnopbc-costat_box_ES*rint(xnopbc/costat_box_ES); 
   *ysipbc=ynopbc-costat_box_ES*rint(ynopbc/costat_box_ES); 
   *zsipbc=znopbc-costat_box_ES*rint(znopbc/costat_box_ES); 

}



/******************************************/
void trajectoria(void)
{
int linies,lin;
float funits;
ofstream trajOut("trajectoria.xyz", std::ios_base::app);  // Fitxer sortida de trajectoria


// UNITATS EN AGSTROMS !!
funits=10;

linies=Num_Eini*2+Num_Sini*2+Num_Oini1+Num_Oini2;

trajOut << linies << endl;
trajOut << "trajectoria RK, iter:" << i << " temps: " << i*ts << " ns" << endl;

// Coordenades de E
for(lin=1;lin<=Num_Eini;lin++)
 { if(typ[lin]==1)
   { trajOut << "E " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "E  -100. -100. -100. " << endl; } 
 }


// Coordenades de S
for(lin=Num_Eini+1;lin<=Num_Eini+Num_Sini;lin++)
 { if(typ[lin]==2)
   { trajOut << "S " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "S  -100. -100. -100. " << endl; } 
 }


// Coordenades de O1
for(lin=Num_Eini+Num_Sini+1;lin<=Num_Eini+Num_Sini+Num_Oini1;lin++)
 { if(typ[lin]==3)
   { trajOut << "O1 " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "O1  -100. -100. -100. " << endl; } 
 }


// Coordenades de O2
for(lin=Num_Eini+Num_Sini+Num_Oini1+1;lin<=Num_Eini+Num_Sini+Num_Oini1+Num_Oini2;lin++)
 { if(typ[lin]==6)
   { trajOut << "O2 " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "O2  -100. -100. -100. " << endl; } 
 }


// Coordenades de C
for(lin=1;lin<=Num_Eini;lin++)
 { if(typ[lin]==4)
   { trajOut << "C " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "C  -100. -100. -100. " << endl; } 
 }


// Coordenades de P
for(lin=Num_Eini+1;lin<=Num_Eini+Num_Sini;lin++)
 { if(typ[lin]==5)
   { trajOut << "P " << px[lin]*funits << " " << py[lin]*funits << " " << pz[lin]*funits << " 0" << endl; } 
   else
   { trajOut << "P  -100. -100. -100. " << endl; } 
 }

}
/******************************************/
void vmd_script(void)
{
ofstream vmdOut("radis3D.tcl");  // Fitxer parametres per VMD

// Carrega el fitxer d'output de la trajectoria
// vmd -e radis3D.tcl per tal de visualitzar-lo

vmdOut << " " << endl;
vmdOut << "mol delete top" << endl;
vmdOut << "mol load xyz trajectoria.xyz" << endl;
vmdOut << "mol delrep 0 top" << endl;
vmdOut << "display resetview" << endl; 
vmdOut << " " << endl;

// Canvia les propietats de cada tipus d'atoms adaptant-les al nostres sistema

vmdOut << "puts \"Radis E S O1 O2 C P\" " << endl;
vmdOut << "puts \"Carrega amb: vmd -e radis3D.tcl\" " << endl;
vmdOut << " " << endl;

vmdOut << "set selE [atomselect top \"name E\"]" << endl;
vmdOut << "$selE set radius " << rad_h_E*10 << endl;

vmdOut << "set selC [atomselect top \"name C\"]" << endl;
vmdOut << "$selC set radius " << rad_h_E*10 << endl;

vmdOut << "set selS [atomselect top \"name S\"]" << endl;
vmdOut << "$selS set radius " << rad_h_S*10 << endl;

vmdOut << "set selP [atomselect top \"name P\"]" << endl;
vmdOut << "$selP set radius " << rad_h_P*10 << endl;

vmdOut << "set selO1 [atomselect top \"name O1\"]" << endl;
vmdOut << "$selO1 set radius " << rad_h_O1*10 << endl;

vmdOut << "set selO2 [atomselect top \"name O2\"]" << endl;
vmdOut << "$selO2 set radius " << rad_h_O2*10 << endl;
vmdOut << " " << endl;

vmdOut << "# Substrat (Color Blau)" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name S" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 0" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;

vmdOut << "# Enzim (Color Vermell)" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name E" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 1" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;

vmdOut << "# Complex" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name C" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 2" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;

vmdOut << "# Producte" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name P" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 3" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;

vmdOut << "# Obstacle 1 (Color Groc)" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name O1" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 4" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;

vmdOut << "# Obstacle 2 (Color Gris)" << endl;
vmdOut << "mol representation VDW 1.000000 16.000000" << endl;
vmdOut << "mol selection name O2" << endl;
vmdOut << "mol material Opaque" << endl;
vmdOut << "mol color ColorID 5" << endl;
vmdOut << "mol addrep top" << endl;
vmdOut << " " << endl;



//LIMITS CAPSA:

vmdOut << "graphics top cylinder {0 0 0} {0 " << costat_box_ES*10 << " 0} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 0 0} {" << costat_box_ES*10 << " 0 0} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 0 0} {0 0 " << costat_box_ES*10 <<"} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 " << costat_box_ES*10 << " 0} {"<< costat_box_ES*10 << " " << costat_box_ES*10 << " 0} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 0 " << costat_box_ES*10 << "} {0 "<< costat_box_ES*10 << " " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {" << costat_box_ES*10 << " 0 0} {" << costat_box_ES*10 << " " << costat_box_ES*10 << " 0} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {" << costat_box_ES*10 << " " << costat_box_ES*10 << " 0} {" << costat_box_ES*10 << " " << costat_box_ES*10 << " " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 " << costat_box_ES*10 <<" " << costat_box_ES*10 << "} {" << costat_box_ES*10 << " " << costat_box_ES*10 << " " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {" << costat_box_ES*10 << " 0 " << costat_box_ES*10 << "} {" << costat_box_ES*10 << " " << costat_box_ES*10 << " " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 " << costat_box_ES*10 << " 0} {0 " << costat_box_ES*10 << " " << costat_box_ES*10 << "}  radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {0 0 " << costat_box_ES*10 << "} {" << costat_box_ES*10 << " 0 " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << "graphics top cylinder {" << costat_box_ES*10 << " 0 0} {" << costat_box_ES*10 << " 0 " << costat_box_ES*10 << "} radius 5 resolution 60 filled yes" << endl;
vmdOut << " " << endl;

// Ves al primer step de la trajectoria
vmdOut << " " << endl;
vmdOut << "animate goto 0" << endl;
vmdOut << " " << endl;

// Canvi de color del fons a blanc

vmdOut << " " << endl;
vmdOut << "color Display Background white" << endl;
}
/******************************************/
void lectura_input(void)
{
char cadena[500];
int num;
time_t rawtime;
struct tm * timeinfo;
long idum_aux;

ifstream file_input ("RK-INPUT.TXT");

if (file_input)
   { }  // FILE OK
   else
   { cout << "ERROR!!! Check: RK-INPUT.TXT"; exit(1);  } // FILE MAL


/* costat_box[nm]= : costat_box_ES */
file_input >> cadena;
file_input >>  costat_box_ES;
file_input.getline(cadena,500);

/* costat_box_Obstacles */

costat_box_Obstacles=costat_box_ES;

/* Num_Eini */
file_input >> cadena;
file_input >>  Num_Eini;
file_input.getline(cadena,500);

/* Num_Sini */
file_input >> cadena;
file_input >>  Num_Sini;
file_input.getline(cadena,500);

/* Num_Oini1 */
file_input >> cadena;
file_input >>  Num_Oini1;
file_input.getline(cadena,500);

/* Num_Oini2 */
file_input >> cadena;
file_input >>  Num_Oini2;
file_input.getline(cadena,500);

/* radi_hidrodinamic E[nm] */
file_input >> cadena;
file_input >>  rad_h_E;
file_input.getline(cadena,500);

/* radi_compacte E[nm] */
file_input >> cadena;
file_input >>  rad_c_E;
file_input.getline(cadena,500);


/* radi_hidrodinamic S[nm] */
file_input >> cadena;
file_input >>  rad_h_S;
file_input.getline(cadena,500);

/* radi_compacte S[nm] */
file_input >> cadena;
file_input >>  rad_c_S;
file_input.getline(cadena,500);


/* radi_hidrodinamic O1[nm] */
file_input >> cadena;
file_input >>  rad_h_O1;
file_input.getline(cadena,500);

/* radi_compacte O1[nm] */
file_input >> cadena;
file_input >>  rad_c_O1;
file_input.getline(cadena,500);

/* radi_hidrodinamic O2[nm] */
file_input >> cadena;
file_input >>  rad_h_O2;
file_input.getline(cadena,500);

/* radi_compacte O2[nm] */
file_input >> cadena;
file_input >>  rad_c_O2;
file_input.getline(cadena,500);

/* radi_hidrodinamic P[nm] */
file_input >> cadena;
file_input >>  rad_h_P;
file_input.getline(cadena,500);

/* radi_compacte P[nm] */
file_input >> cadena;
file_input >>  rad_c_P;
file_input.getline(cadena,500);


/* iteracions: i_max */  
file_input >> cadena;
file_input >>  i_max;
file_input.getline(cadena,500);


/* time step[ns]  */
file_input >> cadena;
file_input >>  ts;
file_input.getline(cadena,500);

/* temps d'equilibraci� del sistema[ns]  */
file_input >> cadena;
file_input >>  i_eq;
file_input.getline(cadena,500);

/* D_E[nm2/ns]   */ 
file_input >> cadena;
file_input >>  D_E;
file_input.getline(cadena,500);

/* D_S[nm2/ns]   */ 
file_input >> cadena;
file_input >>  D_S;
file_input.getline(cadena,500);

/* D_O1[nm2/ns]   */ 
file_input >> cadena;
file_input >>  D_O1;
file_input.getline(cadena,500);

/* D_O2[nm2/ns]   */ 
file_input >> cadena;
file_input >>  D_O2;
file_input.getline(cadena,500);

/* D_P[nm2/ns]   */ 
file_input >> cadena;
file_input >>  D_P;
file_input.getline(cadena,500);

/* Vo[nm2/ns]  */  
file_input >> cadena;
file_input >>  Vo;
file_input.getline(cadena,500);

/* p1 */  // probability of E and S reacting to give C
file_input >> cadena;
file_input >>  p1;
file_input.getline(cadena,500);

/* p2 */ // probability of C decaying into S and E
file_input >> cadena;
file_input >>  p2;
file_input.getline(cadena,500);

/* p3 */ // probability of C decaying into E and P
file_input >> cadena;
file_input >>  p3;
file_input.getline(cadena,500);

/* Repetions: sample_max  */  
file_input >> cadena;
file_input >>  sample_max;
file_input.getline(cadena,500);

/* IDUM */ // Generació de nombres aleatoris
file_input >> cadena;
file_input >> idum;
file_input.getline(cadena,500);

/* freqtr */ // Cada quan imprimir la trajectòria
file_input >> cadena;
file_input >> freqtr;
file_input.getline(cadena,500);

/* Temperatura */ //
file_input >> cadena;
file_input >> T;
file_input.getline(cadena,500);


/* Ur */ // Parametre d'interaccio entre dues particules. Al�ada del shoulder del potencial de interaccio.
file_input >> cadena;
file_input >> Ur;
file_input.getline(cadena,500);


/* Si l'usuari introdueix un 0, l'idum es genera aleatòriament (i la llavor també) */

if(idum==0)
  {

  if (nrank == MASTER)
  {
  time ( &rawtime );
  idum=-abs(rawtime);
  for (icpu=1;icpu<totalcpu;icpu++) {
  idum_aux=idum+icpu;
  MPI_Isend(&idum_aux,sizeof(idum_aux),MPI_CHAR,icpu,1,MPI_COMM_WORLD,&request);  
   }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (nrank != MASTER)
  {
  MPI_Recv(&idum,sizeof(idum),MPI_CHAR,MASTER,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }

  cout << "cpu= " << nrank << " idum= " << idum << endl;
  }
//FINAL DE LA LECTURA//

file_input.close(); 


}


/******************************************/
void output_variables(void)
{


/* costat_box[nm]= : costat_box_ES */
posOut << "#costat_box[nm]= ";
posOut <<  costat_box_ES << endl;

/* Num_Eini */
posOut << "#Enzims=          ";
posOut <<  Num_Eini << endl;

/* Num_Sini */
posOut << "#Substracts=     ";
posOut <<  Num_Sini << endl;

/* Num_Oini1 */
posOut << "#Obstacles1=      ";
posOut <<  Num_Oini1 << endl;

/* Num_Oini2 */
posOut << "#Obstacles2=      ";
posOut <<  Num_Oini2 << endl;

/* radi_E[nm] */
posOut << "#radi_E[nm]=      "; 
posOut <<  rad_h_E << endl;

/* radi_S[nm] */
posOut << "#radi_S[nm]=      "; 
posOut <<  rad_h_S << endl;

/* radi_O1[nm] */
posOut << "#radi_O1[nm]=      "; 
posOut <<  rad_h_O1 << endl;

/* radi_O2[nm] */
posOut << "#radi_O2[nm]=      "; 
posOut <<  rad_h_O2 << endl;

/* radi_P[nm] */
posOut << "#radi_P[nm]=      "; 
posOut <<  rad_h_P << endl;

/* iteracions: i_max */  
posOut << "#iteracions=     ";
posOut <<  i_max << endl;

/* ts[ns]  */
posOut << "#ts[ns]=        ";
posOut <<  ts << endl;

/* temps d'equilibracio del sistema */
posOut << "#i_eq [ns]=        ";
posOut <<  i_eq << endl;

/* D_E[nm2/ns]   */ 
posOut << "#D_E[nm2/ns]=     ";
posOut <<  D_E << endl;

/* D_S[nm2/ns]   */ 
posOut << "#D_S[nm2/ns]=     ";
posOut <<  D_S << endl;

/* D_O1[nm2/ns]   */ 
posOut << "#D_O1[nm2/ns]=     ";
posOut <<  D_O1 << endl;

/* D_O2[nm2/ns]   */ 
posOut << "#D_O2[nm2/ns]=     ";
posOut <<  D_O2 << endl;

/* D_P[nm2/ns]   */ 
posOut << "#D_P[nm2/ns]=     ";
posOut <<  D_P << endl;

/* Vo[nm2/ns]  */  
posOut << "#Vo[nm2/ns]=    ";
posOut <<  Vo << endl;

/* p1 */  // probability of E and S reacting to give C
posOut << "#p1(E+S->C)=    ";
posOut <<  p1 << endl;

/* p2 */ // probability of C decaying into S and E
posOut << "#p2(C->E+S)=    ";
posOut <<  p2 << endl;

/* p3 */ // probability of C decaying into E and P
posOut << "#p3(C->E+P)=    ";
posOut <<  p3 << endl;

/* Repections: sample_max  */  
posOut << "#repeticions=    ";
posOut <<  sample_max << endl;

/* IDUM */
posOut << "#idum=    ";
posOut <<  idum << endl;

/*Frequencia d'escriptura de la trajectoria*/
posOut << "#freqtr=    ";
posOut <<  freqtr << endl;

/* Temperatura */
posOut << "#temperatura=    ";
posOut << T << endl;

/* Constant del potencial quadratic */
posOut << "#kpair=    ";
posOut << kpair << endl;

/* costat_box_Obstacles */

posOut << endl << "#costat_box_Obstacles: " << costat_box_Obstacles << endl;
posOut <<  "#costat_box_ES: " << costat_box_ES << endl << endl;

}
/*******************************************************/
int control_parametres()

{

	if(costat_box_Obstacles<costat_box_ES) 
	{ printf("ERROR: costat_box_Obstacles<costat_box_ES\n"); return(-1); }

	/* control de MAXTEMPS */

	if(i_max>MAXTEMPS)
	{ printf("ERROR: el temps de simulació establert excedeix el límit\n"); return(-1); }

	/* control de MAXPART */

	if((Num_Eini*2+Num_Sini*2+Num_Oini1+Num_Oini2)>MAXPART)
	{ printf("ERROR: el nombre de partícules establert excedeix el límit\n"); return(-1); }

	/* Si D<0 (per tots els casos) D per Stokes-Einstein */

	if(D_E<0)
	{
	D_E=DstokesEinstein(rad_h_E);
	cout<< "The D_E coefficient calculated by Stokes-Einstein is : "<<D_E<< " [nm2/ns]"<<endl;
  	cout << "at T= " << T << "[K] and dynamic viscosity= " << visc << "[kg/s*m]" << endl;
	cout << " " << endl;
	cout << "The average displacement of the Enzyme particles should be " << sqrt(D_E*ts*2)<< " nm"  << endl;
	cout << " " << endl;
    }
 
if(D_S<0)
    {
	D_S=DstokesEinstein(rad_h_S);
	cout<< "The D_S coefficient calculated by Stokes-Einstein is : "<<D_S<< " [nm2/ns]"<<endl;
  	cout << "at T= " << T << "[K] and dynamic viscosity= " << visc << "[kg/s*m]" << endl;
	cout << " " << endl;
        cout << "The average displacement of the Substrate particles should be " << sqrt(D_S*ts*2)<< " nm"  << endl;
	cout << " " << endl;
    }

if(D_P<0)
    {
	D_P=DstokesEinstein(rad_h_P);
	cout<< "The D_P coefficient calculated by Stokes-Einstein is : "<<D_P<< " [nm2/ns]"<<endl;
  	cout << "at T= " << T << "[K] and dynamic viscosity= " << visc << "[kg/s*m]" << endl;
	cout << " " << endl;
        cout << "The average displacement of the Product particles should be " << sqrt(D_P*ts*2) << " nm" << endl;
	cout << " " << endl;
    }
  
if(D_O1<0)
    {
	D_O1=DstokesEinstein(rad_h_O1);
	cout<< "The D_O1 coefficient calculated by Stokes-Einstein is : "<<D_O1<< " [nm2/ns]"<<endl;
  	cout << "at T= " << T << "[K] and dynamic viscosity= " << visc << "[kg/s*m]." << endl;
	cout << " " << endl;
        cout << "The average displacement of the Obstacle 1  particles should be " << sqrt(D_O1*ts*2) << " nm" << endl;
	cout << " " << endl;
    }

if(D_O2<0)
    {
	D_O2=DstokesEinstein(rad_h_O2);
	cout<< "The D_O2 coefficient calculated by Stokes-Einstein is : "<<D_O2<< " [nm2/ns]"<<endl;
  	cout << "at T= " << T << "[K] and dynamic viscosity= " << visc << "[kg/s*m]." << endl;
	cout << " " << endl;
        cout << "The average displacement of the Obstacle 2 particles should be " << sqrt(D_O2*ts*2) << " nm" << endl;
	cout << " " << endl;
    }

return(0);
}
/****************************************************/

 double DstokesEinstein(double radiSE)
{
 double DSE;
 


        radiSE = radiSE*1e-9;
        DSE=kb*T/(6*M_PI*visc*radiSE);
        DSE = DSE*1e9;  /* D en [nm2/ns] */


 return(DSE);
}

/*********************************************/

// double PotencialQuadratic (double diffpbcx, double diffpbcy, double diffpbcz,  double R_n, double R_n1, double *dVx,  double *dVy, double *dVz)
//{

// Calcula el gradient del Potencial Quadratic de cadascuna de les coordenades

//	 *dVx=((kpair * diffpbcx* (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)) -(R_n+R_n1))) / (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2) )));

//         *dVy=((kpair * diffpbcy* (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)) -(R_n+R_n1))) / (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2) )));

//         *dVz=((kpair * diffpbcz* (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)) -(R_n+R_n1))) / (sqrt( pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2) )));

//}


/*********************************************/


/* Potencial CSW antic


double PotencialCSW(double diffpbcx, double diffpbcy, double diffpbcz,  double R_h, double R_c, double *dVx,  double *dVy, double *dVz)
{

// Calcula el gradient del potencial CSW modificat per nom�s ser repulsiu

*dVx= -delta*Ur*diffpbcx*exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))/(Rc*sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))*pow((1+exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))),2))-24*pow(Rc,24)*diffpbcx/pow((pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)),13);

*dVy= -delta*Ur*diffpbcy*exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))/(Rc*sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))*pow((1+exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))),2))-24*pow(Rc,24)*diffpbcy/pow((pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)),13);

*dVz= -delta*Ur*diffpbcz*exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))/(Rc*sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))*pow((1+exp(delta*(sqrt(pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2))-0.8*Rh)/(Rc))),2))-24*pow(Rc,24)*diffpbcz/pow((pow(diffpbcx,2)+pow(diffpbcy,2)+pow(diffpbcz,2)),13);




if ( *dVx  > 100000) {*dVx = 100000;}
if ( *dVx  < -100000) {*dVx = -100000;}
if ( *dVy >  100000) {*dVy =  100000;}
if ( *dVy < -100000) {*dVy = -100000;}
if ( *dVz >  100000) {*dVz =  100000;}
if ( *dVz < -100000) {*dVz = -100000;}





}

*/

double PotencialCSW(double r_x, double r_y, double r_z,  double S_h, double S_c, double *dVx,  double *dVy, double *dVz)
{

*dVx=-1*S_c*Ur*r_x*pow(1/cosh(S_c*(sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2))-(S_c+S_h)/2)),2)/(2*(S_h-S_c)*sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)))-24*pow(S_c,24)*r_x/pow(pow(r_x,2)+pow(r_y,2)+pow(r_z,2),13);

*dVy=-1*S_c*Ur*r_y*pow(1/cosh(S_c*(sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2))-(S_c+S_h)/2)),2)/(2*(S_h-S_c)*sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)))-24*pow(S_c,24)*r_y/pow(pow(r_x,2)+pow(r_y,2)+pow(r_z,2),13);

*dVz=-1*S_c*Ur*r_z*pow(1/cosh(S_c*(sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2))-(S_c+S_h)/2)),2)/(2*(S_h-S_c)*sqrt(pow(r_x,2)+pow(r_y,2)+pow(r_z,2)))-24*pow(S_c,24)*r_z/pow(pow(r_x,2)+pow(r_y,2)+pow(r_z,2),13);


if ( *dVx  > 100000) {*dVx = 100000;}
if ( *dVx  < -100000) {*dVx = -100000;}
if ( *dVy >  100000) {*dVy =  100000;}
if ( *dVy < -100000) {*dVy = -100000;}
if ( *dVz >  100000) {*dVz =  100000;}
if ( *dVz < -100000) {*dVz = -100000;}




}




/*********************************************/

/*********************************************/


// double PotencialAntic()
//
// double dVi;
//
//	a=(Vo/((radi[celln]+radi[celln1])*r));
//	c=exp(-r/(radi[celln]+radi[celln1]));
//	dVi=(a*bipbc*c);
//	dVpot=dVi*R*T/Dif [celln]





/*********************************************/

void inicialitzacio_sistema(void)
{
 div_t divresult_a;


 count_ES = 0;
 count_Cd1 = 0;
 count_Cd2 = 0;
 cell_max= Num_Eini + Num_Sini + Num_Oini1 + Num_Oini2;

 // Distribute E, S and O randomly on a unit square plane
 // of dimensions "costat_box_ES" and for the number of
 // obstacles "costat_box_Obstacles"

// Per Enzime:
 for (celln=1; celln<=Num_Eini; celln++) {
 px [celln] = ran2(&idum)*costat_box_ES;
 py [celln] = ran2(&idum)*costat_box_ES;
 pz [celln] = ran2(&idum)*costat_box_ES;
 radi_h[celln] = rad_h_E;
 radi_c[celln] = rad_c_E;
 Dif[celln] = D_E;
 typ [celln] = 1; // E molecules
 parella [celln] = 0; // Only C have pairs
 }

// Per Substrate:
 for (celln=Num_Eini+1; celln<=Num_Eini+Num_Sini; celln++) {
 px [celln] = ran2(&idum)*costat_box_ES;
 py [celln] = ran2(&idum)*costat_box_ES;
 pz [celln] = ran2(&idum)*costat_box_ES;
 radi_h[celln] = rad_h_S;
 radi_c[celln] = rad_c_S;
 Dif[celln] = D_S;
 typ [celln] = 2;  // S molecules
 parella [celln] = 0;
 }
///

// Per Obstacles tipus 1
 for (celln=Num_Eini+Num_Sini+1; celln<=Num_Eini+Num_Sini+Num_Oini1; celln++) {
 px [celln] = ran2(&idum)*costat_box_Obstacles;
 py [celln] = ran2(&idum)*costat_box_Obstacles;
 pz [celln] = ran2(&idum)*costat_box_Obstacles;
 radi_h[celln] = rad_h_O1;
 radi_c[celln] = rad_c_O1;
 Dif[celln] = D_O1;
 if (D_O1==0) {Dif[celln]=0.1;} // Pel temps d'equilibracio, relaxacio del sistema
 typ[celln] = 3; // Obstacles tipus 1 = index 3
 parella [celln] = 0;
 }

// Per Obstacles tipus 2
 for (celln=Num_Eini+Num_Sini+Num_Oini1+1; celln<=cell_max; celln++) {
 px [celln] = ran2(&idum)*costat_box_Obstacles;
 py [celln] = ran2(&idum)*costat_box_Obstacles;
 pz [celln] = ran2(&idum)*costat_box_Obstacles;
 radi_h[celln] = rad_h_O2;
 radi_c[celln] = rad_c_O2;
 Dif[celln] = D_O2;
 if (D_O2==0) {Dif[celln]=0.1;} // Pel temps d'equilibracio, relaxacio del sistema
 typ[celln] = 6; // Obstacles tipus 2 = index 6
 parella [celln] = 0;
 }


 divresult_a = div (sample,1);
 if (divresult_a.rem == 0) {cout<<"Sample "<<sample<<endl;}


 cE [i] = 0;
 cS [i] = 0;
 cC [i] = 0;
 cP [i] = 0;
 cES [i] = 0;
 Cd1 [i] = 0;
 Cd2 [i] = 0;



}

/*******************************************/

 void interaccions_sistema (void)
 {
 double b1,b2,b3;
 double b1pbc,b2pbc,b3pbc;

 // Molecule Counter
 count_E = 0;
 count_S = 0;
 count_C = 0;
 count_P = 0;

  for (celln=1; celln<=cell_max; celln++) {
    if (typ[celln] == 1) {count_E = count_E + 1;}
    if (typ[celln] == 2) {count_S = count_S + 1;}
    if (typ[celln] == 4) {count_C = count_C + 1;}
    if (typ[celln] == 5) {count_P = count_P + 1;}
    }

 cE [i] = cE [i] + count_E;
 cS [i] = cS [i] + count_S;
 cC [i] = cC [i] + count_C;
 cP [i] = cP [i] + count_P;
 cES [i] = cES [i] + count_ES;
 Cd1 [i] = Cd1 [i] + count_Cd1;
 Cd2 [i] = Cd2 [i] + count_Cd2;

 count_ES = 0;
 count_Cd1 = 0;
 count_Cd2 = 0;

 // Compute the gradient at the current molecule / obstacle position
 // and reaction of S and E particles

for (celln=1; celln<=cell_max; celln++) {

 sum_x = 0.0;
 sum_y = 0.0;
 sum_z = 0.0;

 for (celln1=1; celln1<=cell_max; celln1++) {  /* No es pot celln1=celln+1 pel calcul de gradients */


  if ((typ[celln] > 0)&&(typ[celln1] > 0)) {

    c_op = 0;

    b1 = px [celln] - px [celln1];
    b2 = py [celln] - py [celln1];
    b3 = pz [celln] - pz [celln1];   
 
    difpbc3D(b1,b2,b3,&b1pbc,&b2pbc,&b3pbc);

    r = sqrt( pow(b1pbc,2) + pow(b2pbc,2) + pow(b3pbc,2) );



// Interaction of S and E particles

// Nomes poden reaccionar despres de la fase d'equilibratge
if (i>i_eq && p1!=0.) {

	if (( (typ[celln] == 1)&&(typ[celln1] == 2)&&(r < (radi_h[celln]+radi_h[celln1])) )||( (typ[celln] == 2)&&(typ[celln1] == 1)&&(r < (radi_h[celln]+radi_h[celln1])))) {
		a = ran2(&idum);
		if (a < p1) {
      			if(typ[celln]==1)
			        { typ[celln1]=0;  // Desapareix el S
			          typ[celln]=4;   // Transformacio E->C
			          parella[celln]=celln1;  // Parella de C
			        }
		        else
			        { typ[celln]=0;  // Desapareix el S
			          typ[celln1]=4;   // Transformacio E->C
			          parella[celln1]=celln;  // Parella de C
			        }
        // Switch on reaction operator
	        c_op = 1;
        	count_ES += 1;
	
	      }
    	}
}

  // Collision dynamics for non-reacting particles  

  	if ((celln == celln1)||(c_op == 1)||(r>=(radi_h[celln]+radi_h[celln1]))) {	// No self-interaction
		dVy = 0.0;
		dVz = 0.0;
                dVx = 0.0;
		}
 	else {
        S_c=radi_c[celln]+radi_c[celln1];
        S_h=radi_h[celln]+radi_h[celln1];
        PotencialCSW(b1pbc,b2pbc,b3pbc,S_h,S_c,&dVx,&dVy,&dVz);
        }

	sum_x = sum_x + (dVx);
	sum_y = sum_y + (dVy);
	sum_z = sum_z + (dVz);
       
 }
 }
 
 grad_x [celln] = - (Dif [celln] * sum_x) / (R * T);
 grad_y [celln] = - (Dif [celln] * sum_y) / (R * T);
 grad_z [celln] = - (Dif [celln] * sum_z) / (R * T);



}

//FI interaccions_sistema:
}

/******************************************************/

void reaccions (void)
{
	// Decay reactions

	for (celln=1; celln<=cell_max; celln++) {

		// Decay of C particle into S and E particles 
		if (typ[celln] == 4) {

			a = ran2(&idum);

			if (a < p2) {

				typ[celln] = 1;
				cell_new=parella[celln];
				typ[cell_new]=2;
                        // 3D. Coordenades esf�riques, us dels angles tita i psi ("g" i "j")
                                hrandom=ran2(&idum); // Entre 0 i 1
                                h=hrandom*2.*M_PI;      // angle azimutal, psi

                                grandom=ran2(&idum);
                                g=acos(1.-2.*grandom);  // angle vertical, thita


                                if(px[celln]<costat_box_Obstacles/2.)   // Evitar que es generi una nova particula a la paret
                                        { px [cell_new] = px [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*cos(h);}
                                else
                                        { px [cell_new] = px [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*cos(h);}
                                if(py[celln]<costat_box_Obstacles/2.)
                                        { py [cell_new] = py [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*sin(h);}
                                else
                                        { py [cell_new] = py [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*sin(h);}
                                if(pz[celln]<costat_box_Obstacles/2.)
                                        { pz [cell_new] = pz [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g);}
                                else
                                        { pz [cell_new] = pz [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g);}


				count_Cd1 = count_Cd1 + 1;
			}
		}

		// Decay of C particle into E and P particles 
		if (typ[celln] == 4) {

			a = ran2(&idum);

			if (a < p3) {

				typ[celln] = 1;
				cell_new=parella[celln];
				typ[cell_new]=5;
				radi_h[cell_new] = rad_h_P;
				Dif[cell_new] = D_P;

                                hrandom=ran2(&idum); // Entre 0 i 1
                                h=hrandom*2.*M_PI;      // angle azimutal, psi

                                grandom=ran2(&idum);
                                g=acos(1.-2.*grandom);  // angle vertical, thita


                                if (px[celln]<costat_box_Obstacles/2.)   // Evitar que es generi una nova particula a la paret
                                        { px [cell_new] = px [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*cos(h);}
                                else
                                        { px [cell_new] = px [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*cos(h);}
                                if(py[celln]<costat_box_Obstacles/2.)
                                        { py [cell_new] = py [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*sin(h);}
                                else
                                        { py [cell_new] = py [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g)*sin(h);}
                                if(pz[celln]<costat_box_Obstacles/2.)
                                        { pz [cell_new] = pz [celln] + (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g);}
                                else
                                        { pz [cell_new] = pz [celln] - (radi_h[celln]+radi_h[cell_new]+0.1*(radi_h[celln]+radi_h[cell_new])) * sin(g);}

		      count_Cd2 = count_Cd2 + 1;

    }
  }

 }

}

/*****************************************************/

 double noves_coordenades (void)
 {
 double x_n,y_n,z_n;
 double x_pbc,y_pbc,z_pbc;
 double disp;

 disp=0;

 for (celln=1; celln<=cell_max; celln++) {
  if (typ [celln] != 0) {


 // Generate a Gaussian random number

 do { x1 = 2.0 * ran2(&idum) - 1.0;
      x2 = 2.0 * ran2(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while (( wy >= 1.0 )||( wy == 0.0));

 wy = sqrt( (-2.0 * log( wy ) ) / wy );
 r1 = x1 * wy;

 // Generate a Gaussian random number

 do { x1 = 2.0 * ran2(&idum) - 1.0;
      x2 = 2.0 * ran2(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while (( wy >= 1.0 )||( wy == 0.0 ));

 wy = sqrt( (-2.0 * log( wy ) ) / wy );
 r2 = x1 * wy;

 // Generate a Gaussian random number
 
 do { x1 = 2.0 * ran2(&idum) - 1.0;
      x2 = 2.0 * ran2(&idum) - 1.0;
      wy = x1 * x1 + x2 * x2;
    } while (( wy >= 1.0 )||( wy == 0.0 ));

 wy = sqrt( (-2.0 * log( wy ) ) / wy );
 r3 = x1 * wy;

 
 vx = (pow((2 * Dif[celln]/ts),0.5) * r1) + grad_x [celln];
 vy = (pow((2 * Dif[celln]/ts),0.5) * r2) + grad_y [celln];
 vz = (pow((2 * Dif[celln]/ts),0.5) * r3) + grad_z [celln];

 if ( (vx*ts >= (costat_box_ES*MAXDESP)) || (vy*ts >= (costat_box_ES*MAXDESP)) || (vz*ts >= (costat_box_ES*MAXDESP))) {
    cout << " WARNING in step " << i << " particle " << celln << " has a displacement greater than " << MAXDESP*costat_box_ES << endl;
    cout << " Please consider revising your time step! " << endl;

}

 
 x_n = px [celln] + (vx*ts);
 y_n = py [celln] + (vy*ts);
 z_n = pz [celln] + (vz*ts);
 disp=disp+sqrt(pow(vx*ts,2)+pow(vy*ts,2)+pow(vz*ts,2)); 
 coordpbc3D(x_n,y_n,z_n,&x_pbc,&y_pbc,&z_pbc,celln);

 px [celln] = x_pbc;
 py [celln] = y_pbc;
 pz [celln] = z_pbc;


 } /* De  if (typ [celln] != 0)  */ 

 }


 if((sample==1)&&(i>=(i_eq+1))) {if (i%freqtr==0) {posDisp << i*ts << " " <<  disp/cell_max << " "<< endl;}}
 }

/*****************************************************/
 void create_output(void)
{
 
 div_t divresult_a;
 
 for (i=(i_eq+1); i<=i_max; i++) {

 divresult_a = div (i,1);
 if (divresult_a.rem == 0) { 

 posOut<<(i-i_eq-1)<<" "<<(i-i_eq-1)*ts<<"  "<<cE[i]/sample_max<<" "<<cS[i]/sample_max<<" "<<cC[i]/sample_max<<" "<<cP[i]/sample_max<<"   "<<cES[i]/sample_max<<" "<<Cd1[i]/sample_max<<" "<<Cd2[i]/sample_max<<endl;

 }

 }

 }



/******************************************/

void Dif_coef (void)
{


 posDiff << "# Difussion coefficients obtained from the linearization of RMSD" << endl;
 posDiff << "Please consider checking these values with their respective <r2> curves!!" << endl;
 posDiff<<"D_E     D_S     D_C     D_P     D_O1     D_O2"<< endl;
 
 if (Num_Eini!=0) {D_E_RMSD = D_RegLineal (i_max,ts, r2_average_E)/6.0;}
 if (Num_Sini!=0) {D_S_RMSD = D_RegLineal (i_max,ts, r2_average_S)/6.0;}
 if (p1!=0) {D_C_RMSD = D_RegLineal (i_max, ts, r2_average_C)/6.0;}
 if (p1!=0) {D_P_RMSD = D_RegLineal (i_max, ts, r2_average_P)/6.0;}
 if (Num_Oini1!=0) {D_O1_RMSD= D_RegLineal (i_max, ts, r2_average_O1)/6.0;}
 if (Num_Oini2!=0) {D_O2_RMSD = D_RegLineal (i_max, ts, r2_average_O2)/6.0;}
 
 posDiff <<"  "<< "    "<< D_E_RMSD<<"    " << D_S_RMSD<<"    " << D_C_RMSD<<"    "  << D_P_RMSD <<"     "<< D_O1_RMSD <<"    "<< D_O2_RMSD << endl;
}
/******************************************/
 double microconstantrate (double concC, double concA, double concB, double ts)
{

 double microconstantrate;

	microconstantrate = ((concC)/ts) / (concA * concB);

 return(microconstantrate);

}

/******************************************/
 double r2_sample (double x, double y, double z, double px, double py, double pz,int xparet, int yparet, int zparet)
{
 double r2_sample;

	r2_sample = pow((x-px-(xparet*costat_box_ES)),2)+pow((y-py-(yparet*costat_box_ES)),2)+pow((z-pz-(zparet*costat_box_ES)),2);

return(r2_sample);
}

/******************************************/
 void r2_acumuladors (void)
{
	for (celln=1; celln<=cell_max; celln++) {
		r2_i = r2_sample(x_ini[celln],y_ini[celln],z_ini[celln],px[celln],py[celln],pz[celln],xparet[celln],yparet[celln],zparet[celln]);
		if (typ[celln] == 1) {acum_r2_E [i] += r2_i;}
		if (typ[celln] == 2) {acum_r2_S [i] += r2_i;}
		if (typ[celln] == 4) {acum_r2_C [i] += r2_i;}
		if (typ[celln] == 5) {acum_r2_P [i] += r2_i;}
		if (typ[celln] == 3) {acum_r2_O1 [i] += r2_i;}
		if (typ[celln] == 6) {acum_r2_O2 [i] += r2_i;}

	}

if (cE [i] != 0) {acum_i_E [i] = acum_r2_E [i] / (cE[i]/(Nfin-Nini));}
                if (cS [i] != 0) {acum_i_S [i] = acum_r2_S [i] / (cS[i]/(Nfin-Nini));}
                if (cC [i] != 0) {acum_i_C [i] = acum_r2_C [i] / (cC[i]/(Nfin-Nini));}
                if (cP [i] != 0) {acum_i_P [i] = acum_r2_P [i] / (cP[i]/(Nfin-Nini));}
                if (Num_Oini1 != 0) {acum_i_O1 [i] = acum_r2_O1 [i] / (Num_Oini1);}
                if (Num_Oini2 != 0) {acum_i_O2 [i] = acum_r2_O2 [i] / (Num_Oini2);}



}

/******************************************/
 void r2_averages (void)
{
        for (i=(i_eq+1); i<=i_max; i++) {
                r2_average_E [i] = acum_i_E [i] / sample_max;
                r2_average_S [i] = acum_i_S [i] / sample_max;
                r2_average_C [i] = acum_i_C [i] / sample_max;
                r2_average_P [i] = acum_i_P [i] / sample_max;
                r2_average_O1 [i] = acum_i_O1 [i] / sample_max;
                r2_average_O2 [i] = acum_i_O2 [i] / sample_max;
        }

}

/******************************************/
double D_RegLineal (int i_max, double ts, double r2_average [])
{
 double D_RegLineal, sum_temps_x, sum_r2_y, sum_xy, sum_temps2_x2;
 int count;

 sum_temps_x = 0.0; 
 sum_r2_y = 0.0; 
 sum_xy = 0.0; 
 sum_temps2_x2 = 0.0; 
 count = 0;

	for (i=1; i<=i_max; i++) {
		sum_temps_x += i*ts;
		sum_r2_y += r2_average [i]; 
		sum_xy += i*ts*r2_average [i];
		sum_temps2_x2 += pow(i*ts,2);
		count += 1;
	}


	D_RegLineal = (count*sum_xy - sum_temps_x*sum_r2_y) / ( count*sum_temps2_x2 - pow(sum_temps_x,2));

return(D_RegLineal);
}

/******************************************/
 void create_output_r2 (void)
{

	for (i=(i_eq+1); i<=(i_max-1); i++) {



                posRAVE << (i-i_eq)*ts << "    " <<r2_average_E[i] << "    " << r2_average_S[i] << "    " << r2_average_C[i] << "    "<< r2_average_P[i]<<"    "<< r2_average_O1[i] <<"    " <<r2_average_O2[i] << endl;



	}
}


/*********************************************************/

 #define IM1 2147483563
 #define IM2 2147483399
 #define AM (1.0/IM1)
 #define IMM1 (IM1-1)
 #define IA1 40014
 #define IA2 40692
 #define IQ1 53668
 #define IQ2 52774
 #define IR1 12211
 #define IR2 3791
 #define NTAB 32
 #define NDIV (1+IMM1/NTAB)
 #define EPS 1.2e-7
 #define RNMX (1.0-EPS)

 float ran2(long *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
       k = (*idum)/IQ1;
       *idum=IA1*(*idum-k*IQ1)-k*IR1;
       if (*idum < 0) *idum += IM1;
       if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
    }
   k=(*idum)/IQ1;
   *idum=IA1*(*idum-k*IQ1)-k*IR1;
   if (*idum < 0) *idum += IM1;
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0) idum2 += IM2;
   j = iy/NDIV;
   iy = iv[j]-idum2;
   iv[j] = *idum;
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;}
                                       
/*********************************************************************/

void intercanvi_infor(void){


// Els procesadors envien al MASTER la informacio sobre els promitjos de r2 que han realitzat per cada particula

  // Informacio del enzim

   if (Num_Eini > 0){

 if (nrank!=MASTER) {
 MPI_Isend(&acum_i_E,sizeof(acum_i_E),MPI_CHAR,MASTER,1,MPI_COMM_WORLD,&request);
    }
   // Esperamos a que todos envien sus datos al MASTER
 MPI_Barrier(MPI_COMM_WORLD);

   // El master recibe la informacion y la promedia

 if (nrank==MASTER) {

   for (icpu=1;icpu<totalcpu;icpu++) {
   MPI_Recv(&acum_aux,sizeof(acum_aux),MPI_CHAR,icpu,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  for (i=1; i<=i_max; i++) {
   acum_i_E [i]=acum_i_E [i]+acum_aux[i];
    }
   }

   }
   }

 MPI_Barrier(MPI_COMM_WORLD);
      // Dades del Sustrat
if (Num_Sini > 0){
if (nrank!=MASTER) {
 MPI_Isend(&acum_i_S,sizeof(acum_i_S),MPI_CHAR,MASTER,1,MPI_COMM_WORLD,&request);
    }
  //  Esperamos a que todos envien sus datos al MASTER
  MPI_Barrier(MPI_COMM_WORLD);

   // El master recibe la informacion y la promedia

   if (nrank==MASTER) {

   for (icpu=1;icpu<totalcpu;icpu++) {
   MPI_Recv(&acum_aux,sizeof(acum_aux),MPI_CHAR,icpu,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   for (i=1; i<=i_max; i++) {
    acum_i_S [i]=acum_i_S [i]+acum_aux[i];
    }
   }

   }
   }

 MPI_Barrier(MPI_COMM_WORLD);

        // Dades del Obstacle 1
 if (Num_Oini1 > 0){
 if (nrank!=MASTER) {
   MPI_Isend(&acum_i_O1,sizeof(acum_i_O1),MPI_CHAR,MASTER,1,MPI_COMM_WORLD,&request);
    }
   // Esperamos a que todos envien sus datos al MASTER
  MPI_Barrier(MPI_COMM_WORLD);

   // El master recibe la informacion y la promedia

   if (nrank==MASTER) {

   for (icpu=1;icpu<totalcpu;icpu++) {
   MPI_Recv(&acum_aux,sizeof(acum_aux),MPI_CHAR,icpu,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   for (i=1; i<=i_max; i++) {
    acum_i_O1 [i]=acum_i_O1 [i]+acum_aux[i];
    }
   }

   }
   }

 MPI_Barrier(MPI_COMM_WORLD);
     // Dades del Obstacle 2
 if (Num_Oini2 > 0){
 if (nrank!=MASTER) {
   MPI_Isend(&acum_i_O2,sizeof(acum_i_O2),MPI_CHAR,MASTER,1,MPI_COMM_WORLD,&request);
    }
   // Esperamos a que todos envien sus datos al MASTER
  MPI_Barrier(MPI_COMM_WORLD);

   // El master recibe la informacion y la promedia

   if (nrank==MASTER) {

   for (icpu=1;icpu<totalcpu;icpu++) {
   MPI_Recv(&acum_aux,sizeof(acum_aux),MPI_CHAR,icpu,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   for (i=1; i<=i_max; i++) {
    acum_i_O2 [i]=acum_i_O2 [i]+acum_aux[i];
    }
   }

   }
   }





}
