// BABA AHMED et RAKOTOARIJAONA Andrianarivo

#ifndef EULER_H
#define EULER_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

/****************************************************
 * On a déclaré toute les constantes ici.
 * **************************************************/
double const r=0.05;
double const sig=0.2;
double const S_0=1;
double const T=1;
double const K=1;
double const eps=0.0002;
double const M=4;
double const N=10000;

/*****************************************************
 * La strucure qu'on a créé ici va nous permettre 
 * de renvoyer plusieurs valeur dans la fonction
 * S_T_bar
 * **************************************************/

struct S_return{
	
	double a;  // retourne l-approximation d'Euler avec le pas h_l
	double b;   // retourne l-approximation d'Euler avec le pas h_{l-1}
};


/*****************************************************
 * Ici cette fonction nous a permi de directement utilisé
 * h_l sans passer à des variables intermediaire
 * *****************************************************/
double h_l(int l){
	return T*pow(M,-l);
}


/*********************************************************
 * Simulation d'un uniforme sur [0,1]
 * *****************************************************/
double unif(void){
    static std::mt19937 gen(time(NULL));
	return (1+gen())/(1+(double)gen.max());
	
}

/********************************************************
 * Simulation de variable gaussienne
 * *****************************************************/
double gauss(void){
	
	return sqrt(-2*log(unif()))*cos(2*M_PI*unif());
	
}

/*******************************************************
 * Cette fonction va nous permettre de calculer le 
 * sous jacent à maturité T
 * *****************************************************/
S_return S_T_bar(int l){
	
	S_return S; // structure pour stocker les approximations d'euler
	double S1(S_0),S1_Next(0); // valeur du sous-jacent pour le pas h_l
	double S2(S_0),S2_Next(0); // valeur du sous-jacent pour le pas h_{l-1}
	/*****************************************************************************
	Le tableau ici va nous etre utile pour enregistrer les valeurs des browiniens 
	qu'on utilse dans S1
	*******************************************************************************/
	std::vector<double> tab; //tableau pour stocker les browniens
	double x,t1(0),t2(0); // variable pour stocker les browniens et les 2 instants de discretisation
	int i(0); //compteur
	double somme1(0),somme2(0); // pour calculer S_bar1 et S_bar2
	//Le if l=0 est utile pour initialiser S2=0.
	if(l==0){
		x=gauss(); //simulation d'un gaussien
		S1_Next=S1*(1+r*h_l(l)+sig*sqrt(h_l(l))*x); //application d'Euler pour le pas h_l
		tab.push_back(x);
		somme1+=(S1_Next+S1)/2;
		S1=S1_Next;
		S2=0;

	}else{
		
		while(t1<T){
			x=gauss();
		    S1_Next=S1*(1+r*h_l(l)+sig*sqrt(h_l(l))*x);//application d'Euler pour le pas h_l
			somme1+=(S1_Next+S1)/2;
			S1=S1_Next;
			t1+=h_l(l);//on incremente le temps par le pas
			tab.push_back(x);
		}
		 /**************************************************************
		  * On ré utilise les browniens qu'on a utilisé avec S1
		  * ***********************************************************/
		  while(t2<T){
		  S2_Next=S2*(1+r*h_l(l-1)+sig*sqrt(h_l(l))*(tab[M*i+3]+tab[M*i+2]+tab[M*i+1]+tab[M*i]));//application d'Euler pour le pas h_l-1
			somme2+=(S2_Next+S2)/2;
			S2=S2_Next;
			t2+=h_l(l-1);//on incremente le temps par le pas
			i++;
		}	
	}
	
	S.a=somme1*h_l(l);
	S.b=somme2*h_l(l-1);
	
	
	return S;
	

}


#endif

