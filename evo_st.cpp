#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <numeric>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <unistd.h>
#include <sstream>

#include "/home2/antonio/redneu_v7.2.5/babel_v6/rna.hpp"
#include "/home2/antonio/redneu_v7.2.5/babel_v6/basics.hpp"

using namespace std;

extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/treedist.h>
}

aleatorio alea;

int main(int argc, const char* argv[]){
    rna     rna;
    basic   basic;
    char *st1;
    char *st2;
    char *strus;
    Tree *T_obje;
    Tree *T2;
    int generaciones;
    int n_mol = 1000;
    int l;
    int posicion;
    int plegar[n_mol];
    int muto[n_mol];
    int conta   = 0;
    int a       = 0;
    int max1;
    int ss;
    int dis;
    double max2; 
    double dis_tree;
    double A; 
    double mfe;
    double dis_suma;
    double c;
    double u;
    double r;
    double p;
    double delta = 250;
    double suma;
    double disdos[n_mol];
    std::string pob_inicial;
    std::string objetivo;
    std::string archivo;
    std::string txt;
    std::string crear;
    std::string b;
    std::string STRING;
    std::string path;
    std::string aa;
    std::string str;
    std::string folder;
    std::vector<string > seqs;
    std::vector<string > bolsa_mientras;
    std::vector<int >    dis_seqs;
    std::vector<double > adec;
    std::vector<double > adec_prom;
    std::vector<double > adec_dos;
    std::vector<double > adec_max;
    std::vector<double > robus_prom;
    std::vector<double > robus_seqs;
    std::vector<double > dis_prom;
    std::vector<double > dis_max_seqs;
    std::vector<double > dis_prom_seqs;
    std::vector<double > robus_max;
    std::vector<double > distancias;
    std::vector<double > v;
    std::vector<float  > prob_mfe;
    std::vector<double >::iterator low1;
    std::vector<std::vector<float > > prob_mfe_sd;
    
    long long       semilla = atoi(argv[1]);
    stringstream    sem;

    alea.sid(semilla);
    sem << semilla;
    str = sem.str();

    /**** Parametros ****/
    generaciones = 2000;

    objetivo = "((((.((..((((........)))).(((((.......))))).....(((((.......))))))).))))....";

    /**Poblacion inicial**/
    path = "salida_seqs_caminada_canon-3.txt";
    ifstream infile;
    infile.open(path);
    if(infile.is_open()){
      for(int i = 0; i < 1000; ++i){
          getline(infile,STRING);
          seqs.push_back(STRING);
      }
      infile.close();
    }
    
    pob_inicial = seqs[semilla];
    seqs.clear();

    for(int i = 0; i < n_mol; ++i){
        seqs.push_back(pob_inicial);
        bolsa_mientras.push_back("N");
        distancias.push_back(0);
    }

    /**plegar[] te dice si hubo una mutacion en alguna secuencia. plegar = 1 -> mutacion -> calcular adecuacion; plegar = 0 -> No mutacion -> usar la misma adecuacion que en generacion (i - 1)**/
    for(int i = 0; i < n_mol; ++i){
        plegar[i]   = 1;
        muto[i]     = 0;
    }

    l       = seqs[0].size();
    strus   = new char[l + 1];

    st1         = expand_Full(objetivo.c_str());
    T_obje      = make_tree(st1); //tao
    mfe         = vrna_fold(seqs[0].c_str(), strus);
    st2         = expand_Full(strus);
    T2          = make_tree(st2);
    dis_tree    = tree_edit_distance(T_obje, T2);
    dis_tree    = dis_tree/(2*l); 
    A           = rna.funcion_seleccion_simple_tree(dis_tree);
    adec_prom.push_back(A); //
    free(st1);
    free(st2);
    free_tree(T2);

    for(int i = 0; i < n_mol; ++i){
        adec.push_back(A);
        adec_dos.push_back(A);
    }
    /********* Distancia promedio de la primera generacion *********/
    dis_prom.push_back(0);
    dis_max_seqs.push_back(0);
    dis_prom_seqs.push_back(0);
    /********* robustez promedio de la primera generacion *********/ 
    r = rna.robustez_mutacional(pob_inicial, objetivo);
    robus_prom.push_back(r); 
    robus_max.push_back(r); 
    prob_mfe_sd.push_back(vector<float>() );
    prob_mfe_sd[0].push_back(rna.probab_reperplastic_mfe(pob_inicial, objetivo, delta));
    prob_mfe_sd[0].push_back(0);
    
    /********* Mutar la primera generacion, todas las seqs tienen la misma probabilidad de mutar. *********/ 
    for(int i = 0; i < n_mol; ++i){
        seqs[i] = rna.mutar_seq(seqs[i]);
    }

    /******* REACTOR *******/
    for(int t = 1; t < generaciones; ++t){
        ++conta;
        for(int j = 0; j < n_mol; ++j){
            if(plegar[j] == 1){
                mfe = vrna_fold(seqs[j].c_str(), strus);
                st2 = expand_Full(strus);
                T2  = make_tree(st2);
                dis_tree = tree_edit_distance(T_obje, T2);
                distancias[j] = dis_tree/(2*l); 
                adec[j] = rna.funcion_seleccion_simple_tree(distancias[j]);
                free(st2);
                free_tree(T2); 
            }
        }
        if(conta == 40 || t == (generaciones - 1) ){
            suma  = 0;
            conta = 0;
            A     = 0;
            dis_prom_seqs.push_back(rna.distprom(seqs, dis_seqs));
            max1 = *max_element(dis_seqs.begin(), dis_seqs.end());
            dis_max_seqs.push_back(max1);
            dis_seqs.clear();
            max2 = *max_element(adec.begin(), adec.end());
            adec_max.push_back(max2);

            for(int i = 0; i < n_mol;++i){
                A += adec[i];
                p = rna.probab_reperplastic_mfe(seqs[i], objetivo, delta);
                suma += p;
                prob_mfe.push_back(p);
            }

            adec_prom.push_back(A/n_mol);
            basic.ceros(plegar, n_mol); //llenar de ceros el vector plegar
            dis_suma = 0;
            dis_suma = std::accumulate(distancias.begin(), distancias.end(), 0.0);
            dis_prom.push_back(dis_suma/n_mol);
            prob_mfe_sd.push_back(vector<float>());
            prob_mfe_sd[prob_mfe_sd.size() - 1].push_back(suma/n_mol);
            prob_mfe_sd[prob_mfe_sd.size() - 1].push_back(basic.SD(prob_mfe));
            prob_mfe.clear();
        }
        /* Calcular distancia promedio de las seqs, distancia maxima, robustez promedio y maxima*/
         if(t == (generaciones - 1)){
            for(int i = 0; i < n_mol; ++i){
                robus_seqs.push_back(rna.robustez_mutacional(seqs[i], objetivo));
            }
            max2 = *max_element(robus_seqs.begin(), robus_seqs.end());
            robus_max.push_back(max2);
            r = std::accumulate(robus_seqs.begin(), robus_seqs.end(), 0.0);
            robus_prom.push_back(r/n_mol);
            robus_seqs.clear();
        }
       /******* SEGMENTO NAZI ******/
        p   = 0;
        ss  = 0;
        while(ss < seqs.size()){
            mfe  = vrna_fold(seqs[ss].c_str(), strus);
            dis  = vrna_hamming_distance(strus, objetivo.c_str());
            if(dis != 0){
                seqs.erase(seqs.begin() + ss);
                adec.erase(adec.begin() + ss);
                continue;
            }
            ss += 1;
        }

        v.push_back(adec[0]);
        for(int i = 1; i < adec.size(); ++i){
            v.push_back(v[i - 1]+adec[i]);
        }
        c = v[adec.size() - 1]; //valor de la ultima casilla
        for(int i = 0; i < n_mol; ++i){
            u                   = c*alea.uno();
            low1                = std::lower_bound(v.begin(), v.end(), u);
            posicion            = (low1 - v.begin());
            rna.generador_mutaciones(seqs[posicion], l, &aa, &a);
            bolsa_mientras[i]   = aa;
            plegar[i]           = a;
            disdos[i]           = distancias[posicion];
            adec_dos[i]         = adec[posicion];
            a                   = 0; 
        }
        v.clear();
        seqs.clear();
        adec.clear();
        distancias.clear();
        for(int i = 0; i < n_mol;++i){
            seqs.push_back(bolsa_mientras[i]);
            adec.push_back(adec_dos[i]);
            distancias.push_back(disdos[i]);
        }
    }

    path = "home2/antonio/redneu_v7.2.5/res/st/";

    folder = "mkdir " + path + str;
    system(folder.c_str());
    std::ofstream file;

    folder = path + str;

    archivo = "adec_prom.txt";
    crear   = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < adec_prom.size(); ++i){
        file << adec_prom[i] << endl;
    }
    file.close();

    archivo = "adec_max.txt";
    crear   = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < adec_max.size(); ++i){
        file << adec_max[i] << endl;
    }
    file.close();

    archivo ="dis_prom_seq.txt";
    crear = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < dis_prom_seqs.size(); ++i){
        file << dis_prom_seqs[i]<< endl;
    }
    file.close();


    archivo ="dis_prom_ss.txt";
    crear = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < dis_prom.size(); ++i){
        file << dis_prom[i]<< endl;
    }

    archivo ="dis_max_seqs.txt";
    crear = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < dis_max_seqs.size(); ++i){
        file << dis_max_seqs[i]<< endl;
    }
    file.close();

    archivo ="robus_prom.txt";
    crear = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < robus_prom.size(); ++i){
        file << robus_prom[i]<< endl;
    }
    file.close();

    archivo ="robus_max.txt";
    crear = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < robus_max.size(); ++i){
        file << robus_max[i]<< endl;
    }
    file.close();

    archivo = "prob_mfess.txt";
    crear   = folder + "/"+ archivo;
    file.open(crear.c_str());
    for(int i = 0; i < prob_mfe_sd.size(); ++i){
        file << prob_mfe_sd[i][0] <<" +/- "<< prob_mfe_sd[i][1] << endl;
    }
    file.close();
    free_tree(T_obje);
    delete[] strus;
    return 0;
}
