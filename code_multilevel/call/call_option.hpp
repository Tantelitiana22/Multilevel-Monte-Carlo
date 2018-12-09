// BABA AHMED et RAKOTOARIJAONA Andrianarivo
#ifndef CALL_HPP
#define CALL_HPP

#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <vector>
#include "euler_call.hpp"
/********************************************************************
 * cette structur nous permet de retourner plusieurs variables 
 * comme la moyenne et la variance
 * ******************************************************************/

struct var_mean{
	
	double mean; // retourne l'esperance de p_l - p_{l-1}
	double meanPl; // retourne l'esperance de p_l 
	double var; // retourne la variance de p_l - p_{l-1}
	double varPl; // retourne la variance de p_l 
};
/*********************************************************************
  * Cette fonction va nous permettre de calculer la moyenne et la variance
  * à partir de la fonction S_T_bar qu'on a fait anvant.l ici correspond
  * au niveau ou on se trouve et n corréspond  à la taille de l'echantillion
  * qu'on utilise.Cette fonction retourne var_mean.
  * ******************************************************************/

var_mean V_L_mean_call(int l,int n){
	
	
	double s(0),s1(0); //variable pour stocker les sommes de monte carlo pour retourner l'esperance
	double v(0),v1(0);  //variable pour stocker les sommes^2 de monte carlo pour retourner la variance
	S_return temp; // structure pour stocker la discretisation d'Euler
	var_mean val; // structure pour stocker les differents variance et esperance
	double x; // variable pour stocker les gaussiens
	double buff(0); //variable temporaire
	double buff1(0); //variable temporaire
	
	
		for(int i=0;i<n;i++){
			temp=S_T(l);
			buff=exp(-r)*(std::max(temp.a-K,0.0)-std::max(temp.b-K,0.0));
			s+=buff;
			v+=buff*buff;
			buff1=exp(-r)*(std::max(temp.a-K,0.0));
			s1+=buff1;
			v1+=buff1*buff1;
			
		}
		//Ici on affect à la a une variable de type structure les valeurs
		//qu'on veut renoyer apres on retoure la variable de type var_mean.
		val.mean=s/n;
		val.meanPl=s1/n;
		val.var=v/n -(s/n)*(s/n);
		val.varPl=v1/n -(s1/n)*(s1/n);
		
		return val;
	
}

//fonction h_l pour calculer hl le pas

 double h_l(int l){
	 
	 return T*pow(M,-l);
 }

/***********************************************************************
 * cette fonction va nous permettre de calculer N_l pour cahque niveau
 * l et dont le nombre de niveau est L. le tableau ici est utile
 * car on a utilisé un tableau pour enregistrer les valeurs des variances
 * qui seront utiles pour le calcul de N_l
 * ********************************************************************/
 unsigned long long  N_l_call(int l,int L,double tab[] ){

	
	double temp(0);
	double vl;
	
	if(l>L){
		std::cout<<"Erreur! Verifier que vous ne vous trompez pas car l<L théoriquement\n";
		return EXIT_FAILURE;
	}
	
    if(L==0){
		
		return (int)(2*pow(eps,-2)*tab[0])+1;
	}else{
		
		for(int i=0;i<=L;i++){
			temp+=sqrt(tab[i]/h_l(i));
			
			if(i==l){
				vl=tab[i];
			}
		}
		
	 return (unsigned long long)(2*(1.0/eps)*(1.0/eps)*sqrt(vl*h_l(l))*temp)+1;
		
	}
}

/**********************************************************
 * Ici on a voulu ne pas réutiliser la moyenne pour avoir
 * la valeur de y_l. Donc on a fait cette fonciton pour calculer
 * à cahque niveau l y_l avec une echantillon de taille n
 * ********************************************************/

double Y_l_call(int l, unsigned long long n){
	
    S_return temp;
    double buff(0);
	 if(l==0){
		for(int i=0;i<n;i++){
			temp=S_T(0);
			buff+=exp(-r)*(std::max(temp.a-K,0.0));
		}
		
	
	}else{
		
		for(int i=0;i<n;i++){
			temp=S_T(l);
			buff+=exp(-r)*(std::max(temp.a-K,0.0)-std::max(temp.b-K,0.0));
			
		}
		
	}
	
	return buff/n;
	
}


/************************************************
 * La fonction log4 est utile pour afficher 
 *  le log4 de la variance et la moyenne.
 * *******************************************/

double log4(double x){
	
	return log(x)/log(4);
}

#endif

