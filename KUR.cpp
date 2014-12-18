/* 
 * File:   KUR.cpp
 * Author: ylinieml
 */

/// Evoluationary Algorithm Parameters
/// <PARAM>
#define POPULATION 100 /// EA: mu + lambda
#define ELIMINATE 50 /// EA: lambda
#define GENERATIONS 1000
#define STAT_RUNS 1

/// Method indicators
#define DO_LC 0
#define DO_PACCET 1

#include <cstdlib>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <iostream>
using namespace std;

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector>
#include <list>
#endif

#define PI 3.1415

/// Small Functions/Macros
#define LYRAND (double)rand()/RAND_MAX
double LYrand_norm(double a){
    double theta=LYRAND*2*PI;
    double rsq=-1*a*log(LYRAND);
    double x=rsq*cos(theta);   
    return(x);
}
#define SMALL 0.0001

bool pretty_print = true;

using namespace std;

double vector_median(vector<double> fit){
    double median;
    /// sort vector
    sort(fit.begin(),fit.end());
    /// even or odd size
    int even = (fit.size()+1)%2;
    int odd = fit.size()%2;
    if(even){ /// even case
        int mid = fit.size()/2;
        double m1 = fit.at(mid);
        double m2 = fit.at(mid+1);
        median = (m1+m2)/2;
    }
    if(odd){ /// odd case
        int mid = fit.size()/2;
        median = fit.at(mid);
    }
    return median;
}
double vector_mean(vector<double> fit){
    double sum=accumulate(fit.begin(), fit.end(), 0.0);
    double mean = sum / fit.size();
    return mean;
}

void report(FILE* pFILE, double value) { /// report to text file
    fprintf(pFILE, "%.5f\t", value);
}

void newline(FILE* pFILE) { /// report to text file
    fprintf(pFILE, "\b \b\n");
}


#include "Evo_Agent_KUR.h"
#include "PaCcET.h"

class KURclass{
public:
    double f1(double x1,double x2,double x3);
    double f2(double x1,double x2,double x3);
};

double KURclass::f1(double x1, double x2, double x3){
    double val=0;
    val += -10*exp(-0.2 * sqrt(pow(x1,2)+pow(x2,2)));
    val += -10*exp(-0.2 * sqrt(pow(x2,2)+pow(x3,2)));
    return val;
}
double KURclass::f2(double x1, double x2, double x3){
    double val=0;
    val+=pow(fabs(x1),0.8) + 5*sin(pow(x1,3));
    val+=pow(fabs(x2),0.8) + 5*sin(pow(x2,3));
    val+=pow(fabs(x3),0.8) + 5*sin(pow(x3,3));
    return val;
}

int main(){
    srand(time(NULL));
    
    FILE* pFILE_fit;
    FILE* pFILE_f1;
    FILE* pFILE_f2;
    pFILE_fit=fopen("fitness.txt","w");
    pFILE_f1=fopen("f1.txt","w");
    pFILE_f2=fopen("f2.txt","w");
    
    PaCcET T;
    
    for(int stat_run=0; stat_run < STAT_RUNS; stat_run++) {
        T.Pareto_Reset();
        
        /// Begin Anchor Block
        /// Comment this out for an unseeded P_I^*
        vector<double> f1_undermin;
        vector<double> f2_undermin;
        f1_undermin.push_back(-21);
        f1_undermin.push_back(1);
        f2_undermin.push_back(-13);
        f2_undermin.push_back(-12);
        T.Pareto_Check(f1_undermin);
        T.Pareto_Check(f2_undermin);
        /// End Anchor Block
            
            KURclass environment;
            KURclass* pE = &environment; /// pointer to Environment

            vector<Evo_Agent_KUR> Agents;
            vector<Evo_Agent_KUR>* pVA = &Agents; /// pointer to Vector of Agents

            for (int i = 0; i < POPULATION; i++) {
                Evo_Agent_KUR EA;
                EA.start();
                pVA->push_back(EA);
            }

            for (int gen = 0; gen < GENERATIONS; gen++){
                if (gen % (GENERATIONS / 10) == 0) {
                    cout << "Run No." << stat_run << " is " << (double) gen / GENERATIONS * 100 << " % Complete!" << endl;
                }
                
                /// For each population member, we execute KUR
                for (int mem=0; mem<POPULATION; mem++) {
                    Evo_Agent_KUR* pA = &pVA->at(mem);
                    pA->reset();
                    double f1hold = pE->f1(pA->get_action(0),pA->get_action(1),pA->get_action(2));
                    double f2hold = pE->f2(pA->get_action(0),pA->get_action(1),pA->get_action(2));
                    pA->set_f1(f1hold);
                    pA->set_f2(f2hold);
               }
                
                /// we set up the pareto front:
                for (int a = 0; a < pVA->size(); a++) {
                    vector<double> MO;
                    MO.push_back(pVA->at(a).get_f1());
                    MO.push_back(pVA->at(a).get_f2());
                    T.Pareto_Check(MO);
                }
                
                /// and now we do comparisons.
                for (int a = 0; a < pVA->size(); a++) {
                    /// <PARAM>
                    /// <Linear Combination of Objectives>
                    if(DO_LC){
                    double tr=0.5;
                    pVA->at(a).fitness = (pVA->at(a).get_f1())*tr + pVA->at(a).get_f2()*(1-tr);
                    }
                    /// <PaCcET>
                    vector<double> MO;
                    vector<double>* pMO;
                    vector<double> OMO;
                    pMO = &MO;
                    MO.push_back(pVA->at(a).get_f1());
                    MO.push_back(pVA->at(a).get_f2());
                    OMO=MO; /// Copy for pareto check after transforming

                    if(DO_PACCET){
                    T.execute_N_transform(pMO);
                    //T.Pareto_Check(OMO);
                    pVA->at(a).transformed_fitness = fabs(MO.at(0)-0.2) + fabs(MO.at(1)-0.8);
                    pVA->at(a).fitness = pVA->at(a).transformed_fitness;
                    }
                
                    
            }
                
                vector<double> fit;
                vector<double> f_ones;
                vector<double> f_twos;
                
                for (int a = 0; a < pVA->size(); a++) {
                    fit.push_back(-pVA->at(a).get_fitness()); /// needed for minimization (we always maximize fitness)
                    f_ones.push_back(pVA->at(a).get_f1());
                    f_twos.push_back(pVA->at(a).get_f2());
                   report(pFILE_fit,fit.back()); /// Report every result
                   report(pFILE_f1,f_ones.back()); // Report every result
                   report(pFILE_f2,f_twos.back()); // Report every result
                }
                
                /// always eliminate worst-performing solutions.
                for(int e=0; e<ELIMINATE; e++){
                    double minfit = *min_element(fit.begin(), fit.end());
                    int spot = min_element(fit.begin(), fit.end()) - fit.begin();
                    // kill this fitness;
                    fit.erase(fit.begin() + spot);
                    // kill this agent
                    pVA->erase(pVA->begin() + spot);
                }
                
                /// duplicate best-performing agents
                for(int d=0; d<ELIMINATE; d++){
                    double maxfit = *max_element(fit.begin(), fit.end());
                    int spot = max_element(fit.begin(), fit.end()) - fit.begin();
                    if(rand()%10)
                    {
                        spot=rand()%pVA->size();
                    }
              
                    /// create exact copy
                    pVA->push_back(pVA->at(spot));
                    /// mutate
                    pVA->back().mutate();
                }
    }
            //Start a new line in output file for next run
            fprintf(pFILE_fit,"\n");
            fprintf(pFILE_f1,"\n");
            fprintf(pFILE_f2,"\n");
}      
    fclose(pFILE_f1);
    fclose(pFILE_f2);
    fclose(pFILE_fit);
    T.cout_pareto();
}