// BABA AHMED et RAKOTOARIJAONA Andrianarivo
#include <iostream>
#include <cmath>
#include "call_option.hpp"
#include <fstream>
using namespace std;

int  main(void){
	/*********************************************
	 On fait l'algorithme décrit dans la rapport ici
	 * avec les fonctions définit dans euler.hpp
	 **********************************************/
double VL[10],VLPl[10];
	double mean[10],meanPl[10];
	double s(0);
	double v(0),prix(0);
	double temp;
	 unsigned long long  n[10];
	double y_val[10];
	bool vrai=true;
	/***********************************************
	 * On initialise L=0
	 * *********************************************/
	int L=0;
	double complexite(0);
	double comp_monte_carlo(0);
	
	var_mean valeur;
	//Fichier pour enregistrer les résultats
	 ofstream mon_fichier0("multi0.dat");
	 ofstream mon_fichier1("multi1.dat");
	 ofstream mon_fichier2("meanPl.dat");
     ofstream mon_fichier3("mean.dat");
     ofstream mon_fichier4("nombre_l.dat");
     ofstream mon_fichier5("meanY_val.dat");
     
      /************************************************
     * Boucle qui ne s'arrete pas tant que 
     * la condition de convergence n'est pas verife
     * *******************************************/
    do{    
		//Calcul de VL pour N=10000
		    valeur=V_L_mean_call(L,N);
			VL[L]=valeur.var;
			
			  for(int i=0;i<=L;i++){
				n[i]=N_l_call(i,L,VL);

			 }
			
		// ICi on verifie si la condition N_l>N après on recalibre si c'est le cas

					   for(int i=0;i<=L;i++){
							if(n[i]>N){
								valeur=V_L_mean_call(i,n[i]);
								VL[i]=valeur.var;
								VLPl[i]=valeur.varPl;
								mean[i]=valeur.mean;
								meanPl[i]=valeur.meanPl;
								y_val[i]=Y_l_call(i,n[i]);
							}
						}
				//Si la condition d'avant n'est pas verifier on fait une réinitialisation
					if(n[L]<=N){
					mean[L]=valeur.mean;
					VLPl[L]=valeur.varPl;
					meanPl[L]=valeur.meanPl;
					y_val[L]=Y_l_call(L,N);	
				    }
				   //Ici on a la condition 5 et 6 du pseudo code 
		if(L>=2){
					if(max(abs(y_val[L-1])/M,abs(y_val[L]))<(1.0/sqrt(2))*(M-1)*eps ){
						vrai=false;
					}else{
					L++;
					}
		}else{
			L++;
		 }
		
	}while(vrai);
 
	//Calcul des complexité
	
	complexite=n[0];
     for(int i=1;i<=L;i++){
		 complexite+=n[i]*(pow(M,i)-pow(M,(i-1)));
	 }
    
    for(int i=0;i<=L;i++){
		 comp_monte_carlo+=2*(1/eps)*(1/eps)*VLPl[i]*pow(M,i);
	 }
    
    for(int i=0;i<=L;i++){
		prix+=mean[i];
		mon_fichier0<<i<<" "<<log4(VLPl[i])<<endl;
		mon_fichier1<<i<<" "<<log4(VL[i])<<endl;
		mon_fichier2<<i<<" "<<log4(abs(meanPl[i]))<<endl; 
		mon_fichier3<<i<<" "<<log4(abs(mean[i]))<<endl; 
		mon_fichier4<<i<<" "<<n[i]<<endl;
	}
	
	for(int i=1;i<=L;i++){
		mon_fichier5<<i<<" "<<log4(abs(y_val[i]-y_val[i-1]/M))<<endl; 
	}
	mon_fichier0.close();
	mon_fichier1.close();
	mon_fichier2.close();
	mon_fichier3.close();
	mon_fichier4.close();
	mon_fichier5.close();

  /********************************************************************
    * On affiche le coupe de calcul pour MT et MLMT ainsi que le prix
    * ****************************************************************/
  
  cout<<"eps="<<eps<<";"<<" le cout de calcul pour le multilevel MC : eps^2*cost="<<(eps)*(eps)*complexite<<endl;
  cout<<"eps="<<eps<<";"<<"le cout de calcul pour le MC : eps^2*cost="<<(eps)*(eps)*comp_monte_carlo<<endl;
  cout<<"prix= "<<prix<<endl;
	return 0;
}
