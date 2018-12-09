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
double const eps=0.002;
double const M=4;
double const N=10000;

/*****************************************************
 * La strucure qu'on a créé ici va nous permettre 
 * de renvoyer plusieurs valeur dans la fonction
 * S_T
 * **************************************************/
 
struct S_return{
	
	double a; // retourne l-approximation d'Euler avec le pas h_l
	double b;  // retourne l-approximation d'Euler avec le pas h_{l-1}
};

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

S_return S_T(int l){
	
	double h1,h2,t1(0),t2(0);
	double S1(S_0); // valeur du sous-jacent pour le pas h_l
	double S2(S_0); // valeur du sous-jacent pour le pas h_{l-1}
	/*****************************************************************************
	Le tableau ici va nous etre utile pour enregistrer les valeurs des browiniens 
	qu'on utilse dans S1
	*******************************************************************************/
	std::vector<double> tab; //tableau pour stocker les browniens
	S_return temp,S;
	double x;//variable pour stocker le brownien
	h1=T*pow(M,-l);
	int i(0);
	//Le if l=0 est utile pour initialiser S2=0.
	if(l==0){
         x=gauss();// simulatio du gassien
		S1=S1*(1+r*h1+sig*sqrt(h1)*x);//application d'Euler pour le pas h_l
		S2=0;
        tab.push_back(x);
	}else{	
		
		while(t1<T){
			x=gauss();
			S1=S1*(1+r*h1+sig*sqrt(h1)*x);//application d'Euler pour le pas h_l
			tab.push_back(x);
			t1+=h1;
			
		}
		
		h2=T*pow(M,-(l-1));
		if(l==1){
			S2=S2*(1+r*h2+sig*sqrt(h1)*(tab[0]+tab[1]+tab[2]+tab[3]));//application d'Euler pour le pas h_l-1
		}else{
		/**************************************************************
		  * On ré utilise les browniens qu'on a utilisé avec S1
		  * ***********************************************************/
			while(t2<T){
				S2=S2*(1+r*h2+sig*sqrt(h1)*(tab[M*i+3]+tab[M*i+2]+tab[M*i+1]+tab[M*i]));//application d'Euler pour le pas h_l-1
				t2+=h2;
				i++;
			}
		}
	}
	
	
	S.a=S1;
	S.b=S2;
	
	return S;
		
}


#endif

