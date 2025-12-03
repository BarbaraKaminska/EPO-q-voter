#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <iterator>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

struct Node {
    int pub;
    int priv;
};

double c_priv (Node* network, int N) {
    double c = 0;
    for (int i = 0; i < N; i++){
        c += (1 + network[i].priv);
    }
    return c/2/(double)N;
}

double c_pub (Node* network, int N) {
    double c = 0;
    for (int i = 0; i < N; i++){
        c += (1 + network[i].pub);
    }
    return c/2/(double)N;      
}

double dissonance (Node* network, int N) {
    double d = 0;
    for (int i = 0; i < N; i++){
        d += (1 - network[i].pub*network[i].priv);
    }
    return d/(double)(2*N);   
}

void exp_results (int N, int q, double alpha, double p, double pub, double priv, double diss){
  string fileName = "alpha_EPO_withSA_N" + to_string(N)+ "_q" + to_string(q)+"_alpha"+to_string((int)(alpha*100))+".txt";

  ofstream resFile;
  resFile.open(fileName, ios_base::app);
  resFile<<fixed<<setprecision(3)<<p<<"\t"<<fixed<<setprecision(3)<<pub<<"\t"<<fixed<<setprecision(3)<<priv<<"\t"<<fixed<<setprecision(3)<<diss<<endl;
  resFile.close();
};


int main()
{
    auto start = chrono::high_resolution_clock::now();
    int N = 100000;
    int q = 3;
    int MCS = 5000;
    int MCS_term = 3000;
    int step = 1;
    int number_of_networks = 1;
    double f = 0.5;
    double alpha = 0.1;
    double p_min = 0.0;
    double p_max = 1.;
    double p_step = 0.01;
    double c_0 = 1; // initial condition
    vector <int> nodes;
    vector <int> q_voters;   
    
    
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0);
    uniform_int_distribution<int> int_distribution(0, 1);
    uniform_int_distribution<int> nodes_distribution(0, N - 1);
    random_device rd;
    mt19937 g(rd());

    for (int i = 0; i < N; i++){
        nodes.push_back(i);
    }
    



    for (double p = p_min; p <= p_max + 0.001; p += p_step){
        double c_p_priv = 0;
        double c_p_pub = 0;
        double diss_p = 0; 
        // double d_A_p = 0;
        // double d_T_p = 0;
        double number_of_mes = 0;
        // double number_of_mes_A = 0;        
        // double number_of_mes_T = 0;
        for (int z = 0; z < number_of_networks; z++) {
            Node* network = new Node[N];  

        

            for (int i = 0; i < N; i++) {
                if(distribution(generator) <c_0 ){
                    network[i] = Node{1, 1};
                } else {
                    network[i] = Node{-1, -1};
                }
                
            }
            //cout<<"init"<<dissonance(network, N)<<endl;
            //cout<<"network done"<<endl;

            // Monte Carlo steps            
            for (int mcs = 0; mcs < MCS; mcs++) {
                //cout<<mcs<<endl;
                int node_index;
                double r;

                //single step
                for (int i = 0; i < N; i++) {
                    node_index = nodes_distribution(generator);
                    if (distribution(generator) > alpha){
                        //public update
                        
                        r = distribution(generator);
                        if (r < p) {
                            network[node_index].pub = network[node_index].priv;
                        } else {
                        
                            int j = 0;
                            int neigh; 
                            while (j < q){
                                neigh = nodes_distribution(generator);
                                if (neigh != node_index && find(q_voters.begin(), q_voters.end(), neigh) == q_voters.end()){
                                    q_voters.push_back(neigh);
                                    j++;
                                }
                            }    
                            
                        //disinhibitory contagion
                            if (network[node_index].priv != network[node_index].pub){
                                for (int i = 0; i < q; i++){
                                    if (network[q_voters[i]].pub == network[node_index].priv){
                                        network[node_index].pub = network[node_index].priv;
                                        break;
                                    }
                                }
                            } else {
                            //compliance
                                int opinion = network[q_voters[0]].pub;
                                bool unamity = true; 
                                for (int i = 0; i < q; i++){
                                    if (network[q_voters[i]].pub != opinion){
                                        unamity = false; 
                                        break;
                                    } 
                                }
                                if (unamity){
                                    network[node_index].pub = opinion;
                                }
                            }      

                            q_voters.clear();          
                        }

                    } else {
                        
                        //private update
                        
                        r = distribution(generator);
                        if (r < p) {
                            if (distribution(generator) < f ){network[node_index].priv *= (-1);}                     
                        } else {
                            
                            int j = 0;
                            int neigh; 
                            while (j < q){
                                neigh = nodes_distribution(generator);
                                if (neigh != node_index && find(q_voters.begin(), q_voters.end(), neigh) == q_voters.end()){
                                    q_voters.push_back(neigh);
                                    j++;
                                }
                            }    
                            int opinion = network[q_voters[0]].pub;
                            bool unamity = true; 
                                for (int i = 1; i < q; i++){
                                    if (network[q_voters[i]].pub != opinion){
                                        unamity = false; 
                                        break;
                                    } 
                                }
                                if (unamity){
                                    network[node_index].priv = opinion;
                                }
                            
                            q_voters.clear();
                                        
                        }                    
                    
                    }

                    
                }

                if (mcs % step == 0 && mcs > MCS_term) {
                    double c_priv_ = c_priv(network,  N);
                    double c_pub_ = c_pub(network,  N);
                    double d_ = dissonance(network,  N);
                    c_p_priv += c_priv_;
                    c_p_pub += c_pub_;
                    diss_p += d_;
                    number_of_mes++;
                }

                

            }

            
            
            delete[] network;
        }
        c_p_priv = c_p_priv / number_of_mes;
        c_p_pub = c_p_pub / number_of_mes;
        diss_p = diss_p / number_of_mes;

        std::cout<<p<<"\t"<<c_p_pub<<"\t"<<c_p_priv<<"\t"<<diss_p<< endl;
        exp_results(N, q, alpha, p, c_p_pub, c_p_priv, diss_p);
    }    
    
         
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}
