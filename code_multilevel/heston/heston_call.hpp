// BABA AHMED et RAKOTOARIJAONA Andrianarivo
#ifndef HESTON_H
#define HESTON_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include <vector>

/****************************************************
 * On a déclaré toute les constantes ici.
 * **************************************************/

double const r=0.05;
double const sig=0.2;
double const S_0=1;
double const T=1;
double const K=1;
double const eps=0.0005;
double const M=4;
double const N=10000;
double const V_0=0.04;
double const Rho=-0.5;
double const lambda=5;
double const kappa=0.25;


/*********************************************************
 * Simulation d'un uniforme sur [0,1]
 * *****************************************************/
double unif(void){
    static std::mt19937 gen(time(NULL));
	return (1+gen())/(1+(double)gen.max());
	
}
/*****************************************************
 * On s'est assuré de déclarer deux vecteurs gaussiennes
 * l'un avec cosinus l'autre avec sinus pour s'assurer 
 * de l'indépendance des gaussiennes qu'on va ulitiser
 * ****************************************************/
double gauss1(void){
	
	return sqrt(-2*log(unif()))*cos(2*M_PI*unif());
	
}

double gauss2(void){
	
	return sqrt(-2*log(unif()))*sin(2*M_PI*unif());
	
}

/************************************************************
 * Cette petite classe va nous permettre de faire coréller
 * deux variables aléatoires.
 * ************************************************************/
class Correlation{
	
	protected:
	
	double rho;
	
	public:
	Correlation(double a=Rho):rho(a){}
	
	void correlation(double &A, double &B){
		  A=gauss1();
		  B=A*rho+sqrt(1-rho*rho)*gauss2();
   }
   
  
};

/*****************************************************
 * La strucure qu'on a créé ici va nous permettre 
 * de renvoyer plusieurs valeur dans la fonction
 * S_T
 * **************************************************/

struct S_return{
	double a;
	double b;
};

double h_l(int l){
	return T*pow(M,-l);
}
/*******************************************************
 * Cette fonction va nous permettre de calculer le 
 * sous jacent à maturité T
 * *****************************************************/

S_return S_T(int l){

	Correlation cor(Rho);
	S_return S;
	double W1,W2,v_max;
	/*****************************************************************************
	Les deux tableaux ici va nous etre utile pour enregistrer les valeurs des browiniens 
	qu'on utilse dans S1
	*******************************************************************************/
	std::vector<double> tab1; //tableau pour stocker les browniens
	std::vector<double> tab2; //tableau pour stocker les browniens
	double S1(S_0),V1(V_0);
	double S2(S_0),V2(V_0);
	double temp,x;
	double t1(0),t2(0);
	int i(0);
	//Le if l=0 est utile pour initialiser S2=0et V2=0
	if(l==0){
		cor.correlation(W1,W2);
		S1=S1*(1+r*h_l(l)+sqrt(std::max(V1,0.0)*h_l(l))*W1);
		V1=sig*sig+exp(-lambda*h_l(l))*((V1-sig*sig)+kappa*sqrt(std::max(V1,0.0)*h_l(l))*W2);
		S2=0;
		V2=0;
	}else{
		while(t1<T){
			v_max=std::max(V1,0.0);
			cor.correlation(W1,W2);
			S1=S1*(1+r*h_l(l)+sqrt(v_max*h_l(l))*W1);
			V1=sig*sig+exp(-lambda*h_l(l))*((V1-sig*sig)+kappa*sqrt(v_max*h_l(l))*W2);
			tab1.push_back(W1);
			tab2.push_back(W2);
			t1+=h_l(l);
		}
		 /**************************************************************
		  * On ré utilise les browniens qu'on a utilisé avec S1 et V1
		  * ***********************************************************/
		while(t2<T){
			v_max=std::max(V2,0.0);
			S2=S2*(1+r*h_l(l-1)+sqrt(v_max*h_l(l))*(tab1[M*i+3]+tab1[M*i+2]+tab1[M*i+1]+tab1[M*i]));
			V2=sig*sig+exp(-lambda*h_l(l-1))*((V2-sig*sig)+kappa*sqrt(v_max*h_l(l))*(tab2[M*i+3]+tab2[M*i+2]+tab2[M*i+1]+tab2[M*i]));
			t2+=h_l(l-1);
			i++;
		}
	}

     S.a=S1;
     S.b=S2;
     
     return S;
}















#endif

