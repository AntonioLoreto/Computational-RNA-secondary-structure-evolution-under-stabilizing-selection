#ifndef ecs_v6_hpp
#define ecs_v6_hpp

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "/Users/neo/Desktop/experimentos/Simul_RNA_chidas/babel_v6/rna.hpp"
#include "/Users/neo/Desktop/experimentos/Simul_RNA_chidas/babel_v6/basics.hpp"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;

std::vector<int > estrus_subopt;

class ecs{
public:
	double adec_LE(std::string objetivo, std::string secuencia, int rn, int ligando, int delta);

private: 
};

double ecs::adec_LE(std::string objetivo, std::string secuencia, int rn, int ligando, int delta){	

	basic 		basics;
	Tree 			*T1;
	Tree 			*T2;
	char			*st1;
	char			*st2; 
	char			*strus;
	int 			num_subopt;
	int  			inservible = 0 ; /*callbackfunction*/
	int 			l;
	int 			conta; //para la callbackfunction. definicion global en rna.hpp
	float 		a 	= 0.001;
	float 		be  = 110;
	float 		kas = 0.001;
	double 		kt;
	double   	alfa;
	double 		beta;
	double		suma;
	double		eq;
	double		dis_tree;
	double		particion = 0;
	double		error;
	double		u;
	double		epsilon; 
	double		resultado;
	double		antes;
	double		despues;
	std::vector<double> 					kme;
	std::vector<double> 					probabilidades;
	std::vector<double>  					b;
	std::vector<double>  					x;
	std::vector<double> 					x_pasado;
	std::vector<double> 					col;
	std::vector<double>  					s; 
	std::vector<vector<double> > 	Jacob;
	std::string archivo;

	l 		= secuencia.length();
	strus = strcpy(new char[objetivo.length() + 1], objetivo.c_str()); //string a char*
	st1 	= expand_Full(strus);
	T1 		= make_tree(st1);
	kt 		= (25 + 273.15) * 1.98717 / 1000.0; //kT in kcal/mol. Pa la funcion de particion
	conta = 0;
	free(st1);

	vrna_fold_compound_t *foco; 
	foco = vrna_fold_compound(secuencia.c_str(), NULL, VRNA_OPTION_DEFAULT);
	vrna_subopt_cb(foco, delta, &subopt_callback, (void *)&inservible); 
	vrna_fold_compound_free(foco);
	num_subopt = otipos.size();

 	//calculas las distancias de las suboptimas a la objetivo
	for(int i = 0; i < num_subopt; ++i){
		st2 			= expand_Full(otipos[i].estrus.c_str());
		T2 				= make_tree(st2);
		dis_tree 	= tree_edit_distance(T1, T2); 
		dis_tree 	= dis_tree/(2*l); 
		kme.push_back(a*exp(be*dis_tree)); 
		particion += exp((-1.0*otipos[i].eners)/(kt));
		free(st2);
		free_tree(T2);
	}

	//Probabilidades de cada estructura
  for(int i = 0; i < num_subopt; ++i){
      probabilidades.push_back(exp(-1.0*otipos[i].eners/kt)/particion);
  }

	for(int i = 0; i < num_subopt; ++i){
	    	x.push_back(0); //concentracion inicial complejo
	    	x_pasado.push_back(0);
	}

	suma = 0;
	for(int i = 0; i < num_subopt; ++i){
		suma += x[i]; 
	}

	for(int i = 0; i < num_subopt; ++i){
  	alfa = kas*(2*x[i] + (suma - x[i]) -(rn*probabilidades[i]) -ligando) - kme[i];
  	beta = kas*( x[i] - (rn*probabilidades[i]) ); //rt*P(r_1) = r_1t
		Jacob.push_back(std::vector<double>() );	
		for(int k = 0; k < num_subopt; ++k){
			if(i == k){
				Jacob[i].push_back(alfa);
				continue;
			}
			Jacob[i].push_back(beta);
		}
	}
	
  for(int i = 0; i < num_subopt; ++i){
  	eq = kas*(rn*probabilidades[i] - x[i])*(ligando - suma) - kme[i]*x[i];
  	b.push_back(-eq);
  }
  basics.gauss(Jacob, b, col, s, num_subopt);
  antes = 0;

  /****** Newton-Raphson iteraciones******/
  for(int i = 0; i < 100000; ++i){

  	for(int j = 0; j < num_subopt; ++j){
  		x_pasado[j] = x[j];
  		x[j] 				= x[j] + s[j]; // Calcular x en t+1
  	}
  	//Calcular el error
  	error = 0;
  	for(int j = 0; j < num_subopt; ++j){
  		u 		 = x[j] - x_pasado[j];
  		error += pow(u,2);
  	}
  	despues = error;
  	epsilon = sqrt(error);

		if(fabs(epsilon) <= 0.00000001){
			break;
		}
		if((despues - antes) == 0){
			++conta;
		}
		if(conta == 100){
			cout<<"No Convergencia en secuencia:"<<'\n'<<secuencia<<'\n'<<"estructura objetivo: "<<'\n'<<objetivo<<'\n'<<"parametros: "<<'\n'<<"RNA Total = "<<rn<<'\n'<<"Ligando Total = "<<ligando<<'\n'<<"Delta = "<<delta<<endl;
			return 0;
		}

		antes = despues; 
  	Jacob.clear();
  	b.clear();
  	s.clear();
  	col.clear();

  	//Obtener parciales del jacobiano y llenar la matriz 
  	suma = 0;
		for(int j = 0; j < num_subopt; ++j){
			suma += x[j]; 
		}
		for(int j = 0; j < num_subopt; ++j){
	  	alfa = kas*(2*x[j] + (suma - x[j]) -(rn*probabilidades[j]) -ligando) -kme[j];
	  	beta = kas*(x[j]-(rn*probabilidades[j])); //rt*P(r_1) = r_1t
			Jacob.push_back(vector<double>() );	
			for(int k = 0; k < num_subopt; ++k){
				if(j == k){
					Jacob[j].push_back(alfa);
					continue;
				}
				Jacob[j].push_back(beta);
			}	
		}

  	for(int j = 0; j < num_subopt; ++j){
    	eq = kas*(rn*probabilidades[j] - x[j])*(ligando - suma) - kme[j]*x[j];
    	b.push_back(-eq);
  	}

		basics.gauss(Jacob, b, col, s, num_subopt);
		if(i == 99999){
			cout<<"No convergencia, limite de iteraciones en secuencia:"<<'\n'<<secuencia<<'\n'<<"estructura objetivo: "<<'\n'<<objetivo<<'\n'<<"parametros: "<<'\n'<<"RNA Total = "<<rn<<'\n'<<"Ligando Total = "<<ligando<<'\n'<<"Delta = "<<delta<<endl;
			return 0;
		}
	}
	/*suma = 0;
	for(int i = 0; i < num_subopt; ++i){
		suma += x[i]; 
	}*/
	for(int i = 0; i < num_subopt; ++i){
		//resultado = kas*(rn*probabilidades[i]-x[i])*(ligando-suma)-kme[i]*x[i];
		x[i] = x[i]/ligando;
	}
	suma = 0;
	for(int i = 0; i < num_subopt; ++i){
		suma += x[i]; 
	}

	delete[] strus; 
	otipos.clear();
	b.clear();
	Jacob.clear();
	kme.clear();
	probabilidades.clear();
	x.clear();
	x_pasado.clear();
	col.clear();
	s.clear();
	free_tree(T1);

	return suma; 
}

#endif