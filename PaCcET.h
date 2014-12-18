/* 
 * File:   Procedural_Transformation.h
 * Author: ylinieml
 *
 * Created on October 22, 2013, 11:58 AM
 * Sealed as a black box Feb 18, 2014.
 */

#define PaCcET_VERBOSE 0


#ifndef PROCEDURAL_TRANSFORMATION_H
#define	PROCEDURAL_TRANSFORMATION_H


#include "PaCcET_I_Module.h"
#include "PaCcET_CC_Module.h"

#define PFRONT_THRESHOLD 250
#define PFRONT_BUFFER 1
#define K_FOR_KNN 20

#define OBJECTIVES 2
#define PI 3.1415
#define TIME_DEX 0
#define TREASURE_DEX 1

#define DISCOVERY_THRESHOLD 200

/// <Black Box Line> Last touched August 29, 2014

#ifndef VECTOR_INCLUDE
#define VECTOR_INCLUDE
#include <vector.h>
#endif

double vector_dist(vector<double> a, vector<double> b){
    int A=(int)a.size();
    int B=(int)b.size();
    if (A!=B){
        cout << endl;
        cout << "vector_dist sizes don't match" << endl;
    }
    double dx,dxsq,sum=0,d;
    for(int i=0; i<A; i++){
        dx=a.at(i)-b.at(i);
        dxsq=dx*dx;
        sum+=dxsq;
    }
    d=sqrt(sum);
    return d;
}

class Procedural_Transformation{
    friend class PaCcET_I_Module;
    
    private:
    /// COMPONENTS
    double calc_d(vector<double>);
    double calc_dprime(double, double, double);
    double calc_D(vector<double>);
    double calc_Dprime(vector<double>);
    void eliminate_not_dominating(list<vector<double> >& scPFront_temp, vector<double> td);
    
    vector< vector<double> > Surrogates;
    
    vector< vector<double> > exhaustive_PFront;

    vector<double> input;
    void take_input(vector<double>* coords);
    
    void scale();
    void N_Pro_transform();
    void N_Dummy_transform();
    void calculate_scaled_pareto();
    
    vector<double> output;
    void give_output(vector<double>* coords);
    
    bool COMMAND_I_MODULE;
    bool COMMAND_CC_MODULE;
    
public:
    vector<double> utopia;
    vector<double> nadir;
    vector< vector<double> > PFront;
    vector< vector<double> > scPFront;

    void exhaustive_to_file();
    void PFront_to_file();
    void nad_ut();
    void train();
    
    /// PARETO UTILITIES
    bool Pareto_Check(vector<double>);
    bool Dominated_Check(vector<double> coords);
    int is_in_PFront(vector<double> coords);
    bool does_v1_dominate_v2(vector<double> v1, vector<double> v2);
    void calc_utopia();
    void calc_nadir();
    
    /// PaCcET FUNCTIONALITY
    void take_objective_optimal(vector<double> coords);
    void Pareto_Reset();
    void execute_N_transform(vector<double>* pinputs);
    int get_pareto_size();
    void thresh_PFront();
    void rand_thresh();
    void knn_thresh();
    void anchor_readd();
    
    /// I/O
    void cout_pareto();
    void print_pareto(FILE*);
    void cout_scaled_pareto();
    vector<double> get_ith_pareto_approximate_member(int i);
    
    /// REVERSE PROCESS
    void execute_N_reverse_transform(vector<double>* pinputs);
    void N_Pro_reverse_transform();
    void scale_reverse();
    
    /// Interactive Module
    PaCcET_I_Module IA;
    
    /// Complete Coverage Module
    PaCcET_CC_Module CC;
};

void Procedural_Transformation::exhaustive_to_file(){
    FILE* PFILE;
    cout << "exhaustive in" << endl;
    PFILE=fopen("exhaustive_pareto.txt","w");
    for(int i=0; i<exhaustive_PFront.size(); i++){
        for(int j=0; j<exhaustive_PFront.at(i).size(); j++){
            report(PFILE,exhaustive_PFront.at(i).at(j));
        }
        newline(PFILE);
    }
    fclose(PFILE);
    cout << "exhaustive out" << endl;
}

void Procedural_Transformation::PFront_to_file(){
    FILE* PFILE;
    cout << "Pfront to file in" << endl;
    PFILE=fopen("T_final_front.txt","w");
    for(int i=0; i<PFront.size(); i++){
        for(int j=0; j<PFront.at(i).size(); j++){
            report(PFILE,PFront.at(i).at(j));
        }
        newline(PFILE);
    }
    cout << "Pfront to file out" << endl;
}

vector<double> Procedural_Transformation::get_ith_pareto_approximate_member(int i){
    return PFront.at(i);
}

int Procedural_Transformation::get_pareto_size(){
    return PFront.size();
}

void Procedural_Transformation::Pareto_Reset(){
    PFront.clear();
    scPFront.clear();
    utopia.clear();
    nadir.clear();
    input.clear();
    output.clear();
}

void Procedural_Transformation::cout_pareto(){
    cout << "Current Non-Dominated Set:" << endl;
    int P = PFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            cout << PFront.at(p).at(q) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Procedural_Transformation::print_pareto(FILE* pFILE){
    cout << "Printing pareto to file" << endl;
    pFILE=fopen("pareto_front.txt","a");
    
    /// Identifier line for split between statistical Pareto front runs.
    for(int q=0; q<PFront.at(0).size(); q++){
        fprintf(pFILE, "%.5f\t", -1.98700);
    }
    fprintf(pFILE, "\n");
    
    int P = PFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            fprintf(pFILE, "%.5f\t", PFront.at(p).at(q));
        }
        fprintf(pFILE, "\n");
    }
    fclose(pFILE);
    cout << "done printing pareto to file" << endl;        
}

void Procedural_Transformation::cout_scaled_pareto(){
    cout << endl;
    cout << "Current Non-Dominated Set (NORM):" << endl;
    int P = scPFront.size();
    for(int p=0; p<P; p++){
        for(int q=0; q<PFront.at(p).size(); q++){
            cout << scPFront.at(p).at(q) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

bool Procedural_Transformation::Dominated_Check(vector<double> unscaled_coords){
    /// <Is it dominated by any point in the Pareto front??>
    /// For each Pareto point
    for(int pt=0; pt<PFront.size(); pt++){
        int counter=0;
    for(int obj=0; obj<OBJECTIVES; obj++){
        /// If the pareto point scores higher on a criteria, increment counter
        if(unscaled_coords.at(obj) <= PFront.at(pt).at(obj)){
            counter++;
        }
    }
        if(counter==(OBJECTIVES)){
            /// If the pareto point scored higher or equal on all criteria... 
            /// coords is dominated, and is not a Pareto point.
            /// It is dominated. Return True.
            return true;
        }
    }
    /// It is not dominated. Return False.
    return false;
}

int Procedural_Transformation::is_in_PFront(vector<double> coords){
    for(int i=0; i<PFront.size(); i++){
        int yes=1;
        for(int j=0; j<PFront.at(i).size(); j++){
            if(PFront.at(i).at(j)!=coords.at(j)){yes=0;}
        }
        if(yes){return i;}
    }
    return -1;
}

bool Procedural_Transformation::does_v1_dominate_v2(vector<double> v1, vector<double> v2){
    int counter=0;
    for(int obj=0; obj<OBJECTIVES; obj++){
        /// If v1 scores better on a criteria, increment counter
        if(v1.at(obj) <= v2.at(obj)){
            counter++;
        }
    }
        if(counter==(OBJECTIVES)){
            /// If v1 scored higher or equal on all criteria...
            /// v2 is dominated.
            /// Return True.
            return true;
    }
    /// It is not dominated. Return False.
    return false;
}


bool Procedural_Transformation::Pareto_Check(vector<double> unscaled_coords){
    /// <Display the point in question>
    //cout << "Pareto Checking Point: ";
    //for(int i=0; i<coords.size(); i++){
    //    cout << coords.at(i) << "\t";
    //}
    //cout << endl;
    
    /// <Is it dominated by any point in the Pareto front??>
    /// For each Pareto point
    for(int pt=0; pt<PFront.size(); pt++){
        // determine if coords is dominated:
        if(does_v1_dominate_v2(PFront.at(pt),unscaled_coords)){
            return false;
        }
    }
    
    /// <Does it dominate any points on the Pareto front?>
    vector<int> eliminate;
    for(int pt=0; pt<PFront.size(); pt++){
        if(does_v1_dominate_v2(unscaled_coords,PFront.at(pt))){
            /// If the new point scored higher or equal on all criteria
            /// The "Pareto" point is dominated, and should be eliminated.
            eliminate.push_back(pt); 
        }
    }
    
    /// <Eliminate dominated points on the Pareto Front>
    for(int e=eliminate.size()-1; e>=0; e--){
        /// We eliminate from end -> beginning so that the indices we calculated remain valid.
        int spot = eliminate.at(e);
        PFront.erase(PFront.begin()+spot);
    }
    
    /// <Add new point in correct spot of Pareto Front>
    PFront.push_back(unscaled_coords);
    // also add it to master list.
    exhaustive_PFront.push_back(unscaled_coords);
    
    thresh_PFront();
    
    /// Since Pareto Front has changed, we recalculate the dominated hyperspace.
    nad_ut();
    calculate_scaled_pareto();
    return true;
}

void Procedural_Transformation::thresh_PFront(){
    /// Function: Makes sure that the PFront is maintained at below a threshold size.
    
    /// If PFront is over the threshold
    if(PFront.size()>=PFRONT_THRESHOLD+PFRONT_BUFFER){
        //knn_thresh();
        rand_thresh();
        anchor_readd();
    }
}

void Procedural_Transformation::rand_thresh(){
    while(PFront.size()>=PFRONT_THRESHOLD){
        PFront.erase(PFront.begin()+rand()%PFront.size());
    }
}

void Procedural_Transformation::knn_thresh(){
    /// We find the KNN distance of each.
    vector<double> PFront_KNN_distance;
    int P = PFront.size();
    PFront_KNN_distance.resize(PFront.size());
    for(int piq=0; piq<P; piq++){
        /// find distances from piq to every PFront point.
        vector<double> all_distances;
        for(int i=0; i<P; i++){
            all_distances.push_back(vector_dist(PFront.at(piq),PFront.at(i)));
        }
        /// eliminate K smallest distances
        int K=K_FOR_KNN;
        for(int i=0; i<K; i++){
            int redux = min_element(all_distances.begin(),all_distances.end()) - all_distances.begin();
            all_distances.erase (all_distances.begin()+redux);
        }
        /// remaining minimum is KNN dist.
        int spot = min_element(all_distances.begin(),all_distances.end()) - all_distances.begin();
        PFront_KNN_distance.at(piq)=all_distances.at(spot);
    }
    
    /// PFront_KNN_distance.at(i) now is the ith point on PFront's KNN distance.
    /// Let's eliminate PFRONT_BUFFER amount of PFront and of PFront_KNN_distance.
    /// We eliminate the ones with the lowest KNN distance.
    for(int i=0; i<PFRONT_BUFFER; i++){
        int redux = min_element(PFront_KNN_distance.begin(),PFront_KNN_distance.end()) - PFront_KNN_distance.begin();
        PFront_KNN_distance.erase(PFront_KNN_distance.begin()+redux);
        PFront.erase(PFront.begin()+redux);
    }
    
}

void Procedural_Transformation::anchor_readd(){
    /// We ensure that our anchor points still exist.
    bool match=true;
    bool add=true;
    for(int i=0; i<Anchors.size(); i++){
        add=true;
        for(int p=0; p<PFront.size(); p++){
            match=true;
            for(int obj=0; obj<OBJECTIVES; obj++){
                if(PFront.at(p).at(obj) != Anchors.at(i).at(obj)){
                    match=false;
                    break;
                }
            }
            if(match){
                add=false;
                break;
                /// do not add
            }
        }
        if(add){
            PFront.push_back(Anchors.at(i));
        }
    }
    /// And so it is.
}

void Procedural_Transformation::calculate_scaled_pareto(){
    /// update vector< vector<double> > scPFront;
    //cout << "Calculating Scaled Pareto\t";
    static int counter;
    if(counter%100==0){
    cout << PFront.size() << "\t" ;
    }
    counter++;
    //if(PFront.size() <= 1){return;}
    scPFront.clear();
    for(int i=0; i < PFront.size(); i++){
        vector<double> dual;
        for(int j=0; j < PFront.at(i).size(); j++){
            //cout << "Pfront at i size: " << PFront.at(i).size() << endl;
            //cout << "Nadir size: " << nadir.size() << endl;
            double val = PFront.at(i).at(j);
            double min = utopia.at(j);
            double range = nadir.at(j)-utopia.at(j);
            double scval = (val - min) / range;
            dual.push_back(scval);
        }
        scPFront.push_back(dual);
    }
    //cout_scaled_pareto();
}

void Procedural_Transformation::nad_ut(){
    /// 1) calculate utopia from PFront
    /// 2) calculate nadir from PFront
    
    // 1)
    utopia.clear();
    double min;
    double mindex;
    
    for(int o=0; o<OBJECTIVES; o++){
        min=99999999999;
        mindex=-1;
        for(int p=0; p<PFront.size(); p++){
        if(PFront.at(p).at(o) < min){
            min=PFront.at(p).at(o);
            mindex=p;
        }
        }
        utopia.push_back(PFront.at(mindex).at(o)-0.001); 
    }
    
    // 2)
    nadir.clear();
    double max;
    double maxdex;
    for(int o=0; o<OBJECTIVES; o++){
        max=-999999999999;
        maxdex=-1;
        for(int p=0; p<PFront.size(); p++){
        if(PFront.at(p).at(o)>max){
            max=PFront.at(p).at(o);
            maxdex=p;
        }
        }
        nadir.push_back(PFront.at(maxdex).at(o));   
    }
    /// Added constant avoids atan2 problems at the edges of the PFront.
    
    //cout << "Nadir: " << nadir.at(0) << "\t" << nadir.at(1) << endl;
    //cout << "Utopia: " << utopia.at(0) << "\t" << utopia.at(1) << endl;
}

void Procedural_Transformation::N_Dummy_transform() {
    output.clear();
    for(int i=0; i<input.size(); i++){
        output.push_back(input.at(i));
    }
}

void Procedural_Transformation::N_Pro_transform() {
    
    vector<double> deltas;
    int I=input.size();
    
    //input = scPFront.at(rand()%scPFront.size());
    
    for(int i=0; i<I; i++){
        deltas.push_back(input.at(i));
        //cout << "DELTAS: " << i << " " << deltas.at(i) << endl;
    }
    
    vector<double> directional_ratios;
    
    double sumdeltassq=0.0;
    for(int i=0; i<I; i++){
        sumdeltassq+=deltas.at(i)*deltas.at(i);
    }
    sumdeltassq=sqrt(sumdeltassq);
    for(int i=0; i<I; i++){
        /// More truly: directional cosines.
    directional_ratios.push_back(deltas.at(i)/sumdeltassq);
    //cout << "DIRECTIONAL COSINES " << i << " :: " << directional_ratios.back()<< endl;
    }
    
    /// <FIND d>
    double d;
    d = calc_d(deltas);
    if(PaCcET_VERBOSE > 0){
    cout << "d: " << d << endl;
    }
    
    /// <FIND D>
    double D;
    D = calc_D(directional_ratios);
    if(PaCcET_VERBOSE > 0){
    cout << "D: " << D << endl;
    }
    
    /// <FIND D'>
    double Dprime;
    Dprime = calc_Dprime(directional_ratios);
    if(PaCcET_VERBOSE > 0){
    cout << "Dprime: " << Dprime << endl;
    }
    
    /// <CALCULATE d'>
    double dprime = d*Dprime/D;
    if(PaCcET_VERBOSE > 0){
    cout << "dprime: " << dprime << endl;
    }
    
    output.clear();
    for(int i=0; i<I; i++){
        output.push_back(dprime*directional_ratios.at(i));
    }
    
    
}

void Procedural_Transformation::N_Pro_reverse_transform(){
    /// NOTATION DISAGREEMENT DUE TO BEING REVERSE PROCESS. (intended)
    
     vector<double> deltas;
    int I=input.size();
    for(int i=0; i<I; i++){
        deltas.push_back(input.at(i));
    }
    
    vector<double> directional_ratios;
    double sumdeltassq=0.0;
    for(int i=0; i<I; i++){
        sumdeltassq+=deltas.at(i)*deltas.at(i);
    }
    sumdeltassq=sqrt(sumdeltassq);
    for(int i=0; i<I; i++){
        /// More truly: directional cosines.
    directional_ratios.push_back(deltas.at(i)/sumdeltassq);
    //cout << "DIRECTIONAL COSINES " << i << " :: " << directional_ratios.back()<< endl;
    }
    
    /// <FIND d>
    double dprime;
    dprime = calc_d(deltas);
    
    /// <FIND D>
    double D;
    D = calc_D(directional_ratios);
    
    /// <FIND D'>
    double Dprime;
    Dprime = calc_Dprime(deltas);
    
    /// <CALCULATE d'>
    double d = dprime*Dprime/D;
    
    output.clear();
    for(int i=0; i<I; i++){
        output.push_back(d*directional_ratios.at(i));
    }
}

double Procedural_Transformation::calc_d(vector<double> deltas){
    int I=input.size();
    
    double d=0;
    for(int i=0; i<I; i++){
        d+=deltas.at(i)*deltas.at(i);
    }
    d=sqrt(d);
    return d;
}

void Procedural_Transformation::eliminate_not_dominating(list<vector<double> >& scPFront_temp, vector<double> td){
    
    if(scPFront_temp.size() ==1)
    {
        //cout << "SHORTCUT!" << endl;
        /// only one is dominating, no need to reduce any more.
        return;
    }
    
    bool dominated;
    
    //cout << "ELIMINATE IN: " << scPFront_temp.size() << endl;
    //for(int i=0; i< td.size(); i++){
        //cout << td.at(i) << "\t";
    //}
    //cout << endl;
    //for (std::list< vector<double> >::iterator it=scPFront_temp.begin(); it != scPFront_temp.end(); ++it){
        //vector<double> aaa = *it;
        //for(int j=0; j<aaa.size(); j++){
            //cout << "XXX " << j << " " << aaa.at(j) << endl;
        //}
    //}
    
    for (std::list< vector<double> >::iterator it=scPFront_temp.begin(); it != scPFront_temp.end(); ++it){
        dominated = does_v1_dominate_v2(*it,td);
        if(dominated==true){
            // nothing happens
            //cout << "STILL ALIVE" << endl;
        }
        if(dominated==false){
            // we don't need to use this pfront point for further calc_D calculations on this iteration:
            scPFront_temp.erase(it);
            //cout << "ELIMINATED!" << endl;
        }
    }
    
    //cout << "ELIMINATE OUT: " <<scPFront_temp.size() << endl;
    //for(int p=0; p<scPFront_temp.size(); p++){
    //    dominated = does_v1_dominate_v2(scPFront_temp.at(p),td);
    //    if(dominated==true){
            // nothing happens
     //   }
     //   if(dominated==false){
            // we don't need to use this pfront point for further calc_D calculations on this iteration:
    //        scPFront_temp.erase(scPFront_temp.begin()+p);
     //   }
//}
}

double Procedural_Transformation::calc_D(vector<double> directional_ratios){
    static int iii;
    iii++;
    
    /// normalize deltas:
    /*
    double sumdirr=0.0;
    for( int i=0; i< directional_ratios.size(); i++){
        sumdirr += directional_ratios.at(i);
    }
    for( int i=0; i< directional_ratios.size(); i++){
        directional_ratios.at(i)/=sumdirr;
    }
     */
    
    //return 0.1; 6-> 18
    int I=input.size();
    /// <FIND D>
    double D;
    /// find distance along utopia->input vector which is first dominated.
    double lowerbound = 0;
    double upperbound = 0;
    for(int i=0; i<I; i++){
        upperbound+=1;
    }
    upperbound = 2*sqrt(upperbound)+0.1;
    double candidate;
    double margin=upperbound-lowerbound;
    bool dominated;
    vector<double> td;
    td.resize(I);
    //vector< vector<double> > scPFront_temp = scPFront;
    list< vector<double> > scPFront_temp;
    std::copy( scPFront.begin(), scPFront.end(), std::back_inserter( scPFront_temp ) );

    if(IA.is_active==true){
        IA.scale_surrogates(utopia,nadir);
        std::copy( IA.scsurrogates.begin(), IA.scsurrogates.end(), std::back_inserter( scPFront_temp ) );
    }
    if(CC.is_active==true){
        CC.scale_surrogates(utopia,nadir);
        std::copy( CC.scsurrogates.begin(), CC.scsurrogates.end(), std::back_inserter( scPFront_temp) );
    }
    
    
    static int jjj;

    //return 0.1; 0->6
    while(margin>0.0001){
                jjj++;
        dominated=false;
        candidate=(upperbound+lowerbound)/2;
        for(int i=0; i<I; i++){
        td.at(i) = candidate*directional_ratios.at(i);
        }
        
        for (std::list< vector<double> >::iterator it=scPFront_temp.begin(); it != scPFront_temp.end(); ++it){
            dominated = does_v1_dominate_v2(*it,td);
            if(dominated==true){break;}
        }
     
        //for(int p=0; p<scPFront_temp.size(); p++){
        //dominated = does_v1_dominate_v2(scPFront_temp.at(p),td);
        //if(dominated==true){break;} // once we know it is dominated, we don't need to continue calculating.
        //}
        
        if(dominated==true){upperbound = candidate; eliminate_not_dominating(scPFront_temp,td);}
        if(dominated==false){lowerbound = candidate;}
        margin=upperbound-lowerbound;
        //cout << "SCPFRONT TEMP SIZE: " << scPFront_temp.size() << endl;
    }
    D=(upperbound+lowerbound)/2;
    
    //cout << "INSIDE CALC_D: " << iii << " , " << jjj << endl;
    return D;
    
}

double Procedural_Transformation::calc_Dprime(vector<double> directional_ratios){
    int I=input.size();
    /// <FIND D'>
    double Dprime;
    
    vector<double> td;
    td.resize(I);
    
    double lowerbound = 0;
    double upperbound = 0;
    for(int i=0; i<I; i++){
        upperbound+=1;
    }
    upperbound = 2*sqrt(upperbound)+0.1;
    double margin=upperbound-lowerbound;
    while(margin>0.0001){
        bool dominated=false;
        double candidate=(upperbound+lowerbound)/2;
        for(int i=0; i<I; i++){
            td.at(i) = candidate*directional_ratios.at(i);
        }
        if(accumulate(td.begin(),td.end(),0.0) <= 1){dominated = true;}
        if(dominated==true){lowerbound = candidate;}
        if(dominated==false){upperbound = candidate;}
        margin=upperbound-lowerbound;
    }
    Dprime = (upperbound+lowerbound)/2;
    return Dprime;
}

double Procedural_Transformation::calc_dprime(double d, double D, double Dprime){
    double dprime = d*Dprime/D;
    return dprime;
}

void Procedural_Transformation::take_input(vector<double>* pcoords){
    /// Assign input to class variable.
    // cout << "inside take_input\t" <<  pcoords->size() << endl;
    
    if(pcoords->size() != OBJECTIVES){cout << "Are we doing a " << OBJECTIVES << " objective problem or not?" << endl;}
    input.clear();
    for(int obj=0; obj<OBJECTIVES; obj++){
        input.push_back(pcoords->at(obj));
    }
    // cout << "inside take_input2\t" <<  input.size() << endl;
}

void Procedural_Transformation::give_output(vector<double>* pcoords){
    /// Assign output to external variable.
    for(int obj=0; obj<OBJECTIVES; obj++){
        pcoords->at(obj) = output.at(obj);
    }
}

void Procedural_Transformation::scale(){
    /// Scale values of objectives to be less than one, with the nadir point taking on (1,1).
    if(PaCcET_VERBOSE > 0){
    cout << "SCALING BEFORE!\t";
    cout << input.at(0) << "," << input.at(1) << endl;
    }
    for(int obj=0; obj<OBJECTIVES; obj++){
        double val = input.at(obj);
        double range = nadir.at(obj)-utopia.at(obj);
        double scval = (val - utopia.at(obj)) / range;
        input.at(obj) = scval;
    }
    if(PaCcET_VERBOSE > 0){
    cout << "SCALING AFTER!\t";
    cout << input.at(0) << "," << input.at(1) << endl;
    }
}

void Procedural_Transformation::scale_reverse(){
    /// Reverse scaling, for the reverse transformation
    for(int obj=0; obj<OBJECTIVES; obj++){
        //cout << "utopia " << obj << " = " << utopia.at(obj) << endl;
        //cout << "\t nadir " << obj << " = " << nadir.at(obj) << endl;
        double range = nadir.at(obj)-utopia.at(obj);
        //cout << "range : " << range << endl;
        //cout << "output before: " << output.at(obj) << endl;
        output.at(obj) = output.at(obj) * range + utopia.at(obj);
        
        //output.at(obj) = (output.at(obj) - nadir.at(obj)) / range;
        //cout << "output after: " << output.at(obj) << endl;
    }
}

void Procedural_Transformation::execute_N_transform(vector<double>* pinputs){
    take_input(pinputs);
    scale();
    N_Pro_transform();
    give_output(pinputs);
}

void Procedural_Transformation::execute_N_reverse_transform(vector<double>* pinputs){
    take_input(pinputs);
    N_Pro_reverse_transform();
    scale_reverse();
    give_output(pinputs);
}

void Pro_Pareto_Filter_Testing(){
    Procedural_Transformation T;
    vector<double> coords;
    for(int i=0; i<1000; i++){
    cout << "Trial " << i << endl;
    coords.push_back(rand()%100000);
    coords.push_back(rand()%100000);
    cout << "coords 0\t" << coords.at(0) << endl;
    cout << "coords 1\t" << coords.at(1) << endl;
    T.Pareto_Check(coords);
    coords.clear();
    T.cout_pareto();
    }
    cout << "Placement" << endl;
    coords.push_back(-100000);
    coords.push_back(-100000);
    T.Pareto_Check(coords);
    T.cout_pareto();
    coords.clear();
}

void Procedural_Testing(){
    Procedural_Transformation T;
    vector<double> coords;
    vector<double>* pcoords = &coords;
    for(int i=0; i<100; i++){
    cout << "Trial " << i << endl;
    coords.push_back(rand()%1000);
    coords.push_back(rand()%1000);
    cout << "coords 0:\t" << coords.at(0) << endl;
    cout << "coords 1:\t" << coords.at(1) << endl;
    T.Pareto_Check(coords);
    T.execute_N_transform(pcoords);
    cout << "coords after 0\t" << coords.at(0) << endl;
    cout << "coords after 1\t" << coords.at(1) << endl;
    coords.clear();
    T.cout_pareto();
    T.cout_scaled_pareto();
    }
    cout << "Placement" << endl;
    coords.push_back(-1);
    coords.push_back(-1);
    T.Pareto_Check(coords);
    T.cout_pareto();
    T.cout_scaled_pareto();
    coords.clear();
}

#endif	/* PROCEDURAL_TRANSFORMATION_H */

