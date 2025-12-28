#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;

double MAX_DOUBLE = std::numeric_limits<double>::max();

struct node{
    int city_id;
    int x_coord;
    int y_coord;
};


void distance_array_calculation(vector<node> &nodeset,double **distance_array);
void pheromones_array_initialization(int array_size,double **pheromones_array);
double probability_generator();
vector<int> construct_ant_solution(int start_node, int num_of_nodes,double alpha, double beta, double **distance_array, double **pheromones_array);
double evaluation(vector<int> &nodeset, double **distance_array);
void update_pheromones(vector<vector<int>> &all_nodeset,vector<double> &all_nodeset_total_distance, int num_of_node, int num_of_ant, double **pheromones, double rho, double Q);
string print_to_file(vector<int> &min_nodeset, double average_min_distance, double global_min_distance);
void plot(vector<int> average_min_nodeset,vector<node> &nodeset, string &output_file_address);


int main(){
    srand(time(NULL));

    int city_id, x_coord, y_coord;
    ifstream datafile;
    string file_address;
    cout << "Input the input file address or filename:"<<endl;
    cin >> file_address;
    datafile.open(file_address);
    vector<node> nodeset;
    while(datafile >> city_id >> x_coord >> y_coord){
        node n;
        n.city_id = city_id;
        n.x_coord = x_coord;
        n.y_coord = y_coord;
        nodeset.push_back(n);
    } 

    //----------------------------------------------------------intitializtion----------------------------------------------------------
    int run_time = 30, num_of_iteration = 30, max_evaluation = 10000 * nodeset.size(), num_of_ant = 30;
    double alpha =2.0, beta = 3.0,  rho = 0.7, Q = 1.0;
    double average_min_distance=0; // average of each iteration's min distance
    double global_min_distance = MAX_DOUBLE; // Min distance of all iteration
    vector<int> global_min_nodeset; // best route of all iteration
    double **distance_array = new double *[nodeset.size()];
    for(int i=0;i<nodeset.size();i++){
        distance_array[i] = new double[nodeset.size()];
    }
    distance_array_calculation(nodeset,distance_array);

    double **pheromones_array = new double *[nodeset.size()];
    for(int i=0;i<nodeset.size();i++){
        pheromones_array[i] = new double[nodeset.size()];
    }
    //----------------------------------------------------------intitializtion----------------------------------------------------------

    //----------------------------------------------------------execution----------------------------------------------------------
    for(int i=0;i<run_time;i++){
        //Initialize pheromone of all route to 1
        pheromones_array_initialization(nodeset.size(),pheromones_array);
        vector<int> run_min_nodeset;
        double run_min_distance = MAX_DOUBLE;
        for(int j=0;j<num_of_iteration;j++){
            vector<vector<int>> local_nodeset_list; // Every ant's route
            vector<int> iteration_min_nodeset; // Minimum ant's route for that iteration
            vector<double> local_min_distance_list;// Distance for every ant's route
            double iteration_min_distance = MAX_DOUBLE;// Minimum ant's route distance for that iteration
            for(int k=0;k<num_of_ant;k++){
                int start_node = rand() % nodeset.size();// starting city is randomly generated
                local_nodeset_list.push_back(construct_ant_solution(start_node, nodeset.size(), alpha, beta, distance_array, pheromones_array));
                double current_ant_distance = evaluation(local_nodeset_list[k],distance_array);
                local_min_distance_list.push_back(current_ant_distance);
                if(current_ant_distance < iteration_min_distance){
                    iteration_min_distance = current_ant_distance;
                    iteration_min_nodeset = local_nodeset_list[k];
                }
            }
            update_pheromones(local_nodeset_list,local_min_distance_list,nodeset.size(),num_of_ant,pheromones_array,rho,Q);
            if(iteration_min_distance < run_min_distance){
                run_min_distance = iteration_min_distance;
                run_min_nodeset = iteration_min_nodeset;
            }
        }
        if(run_min_distance < global_min_distance){
            global_min_distance = run_min_distance;
            global_min_nodeset = run_min_nodeset;
        }
        average_min_distance += run_min_distance;
    }
    average_min_distance /= run_time;
    //----------------------------------------------------------execution----------------------------------------------------------
    string output_file_address;
    output_file_address = print_to_file(global_min_nodeset, average_min_distance, global_min_distance);
    plot(global_min_nodeset, nodeset, output_file_address);
    datafile.close();
}

void distance_array_calculation(vector<node> &nodeset,double **distance_array){
    for(int i=0;i<nodeset.size();i++){
        distance_array[i][i] = 0;
        for(int j=i+1;j<nodeset.size();j++){
            distance_array[i][j] = sqrt(pow((nodeset[i].x_coord - nodeset[j].x_coord),2) + pow((nodeset[i].y_coord - nodeset[j].y_coord),2));
            distance_array[j][i] = distance_array[i][j]; 
        }
    }
    return;
}

void pheromones_array_initialization(int array_size, double **pheromones_array){
    for(int i=0;i<array_size;i++){
        pheromones_array[i][i] = 0;
        for(int j=i+1;j<array_size;j++){
            pheromones_array[i][j] = 1;
            pheromones_array[j][i] = 1;
        }
    }
    return;
}

double probability_generator(){
    double probability = (double)rand() / RAND_MAX;
    return probability;
}

vector<int> construct_ant_solution(int start_node, int num_of_nodes,double alpha, double beta, double **distance_array, double **pheromones_array){
    vector<int> local_nodeset;
    local_nodeset.push_back(start_node); // start city is pushed into route
    vector<bool> visited_node(num_of_nodes,false); 
    visited_node[start_node] = true;// start city is visited
    int current_node = start_node;// Ant's current locate city

    // To see if the route of ant is constructed
    while(local_nodeset.size() < num_of_nodes){
        vector<int> node_not_visited;
        vector<double> numerator_list;
        double denominator=0;

        // calculation of probabilities information
        for(int next_node=0;next_node<num_of_nodes;next_node++){
            if(!visited_node[next_node]){
                double numerator;
                double heuristic_information = 1/distance_array[current_node][next_node]; //n(ij) = 1 / d(ij)
                numerator = pow(pheromones_array[current_node][next_node],alpha) * pow(heuristic_information,beta); // tau(ij)^alpha * n(ij)^beta
                numerator_list.push_back(numerator);
                node_not_visited.push_back(next_node);
                denominator += numerator;
            }
        }

        int next_node_selection = -1;
        if(denominator > 0){
            double random_prob = probability_generator(); // probability = 0~1
            double cumulative_prob = 0; // 累積機率
            for(int next_node_index=0;next_node_index<node_not_visited.size();next_node_index++){
                cumulative_prob += numerator_list[next_node_index] / denominator;
                // 累積機率超果 random probability => choose that city as next city
                if(random_prob <= cumulative_prob){
                    next_node_selection = node_not_visited[next_node_index];
                    break;
                }
            }
            
            //If selected city is valid
            if(next_node_selection != -1){
                visited_node[next_node_selection] = true; // change the selected node to visited
                local_nodeset.push_back(next_node_selection); //update ant's route
                current_node = next_node_selection; // change current city to next city for next iteration 
            }
        }
    }
    return local_nodeset;
}

double evaluation(vector<int> &nodeset, double **distance_array){
    double total_distance=0;
    for(int i=0;i<nodeset.size();i++){
        total_distance += distance_array[nodeset[i]][nodeset[(i+1) % nodeset.size()]];
    }
    return total_distance;
}

void update_pheromones(vector<vector<int>> &all_nodeset,vector<double> &all_nodeset_total_distance, int num_of_node, int num_of_ant, double **pheromones, double rho, double Q){
    //tau(ij) = (1-rho) * tau(ij)
    for(int i=0;i<num_of_node;i++){
        for(int j=i+1;j<num_of_node;j++){
            pheromones[i][j] *= (1 - rho);
            pheromones[j][i] = pheromones[i][j]; 
        }
    }
    
    for(int k=0;k<num_of_ant;k++){
        double delta = Q / all_nodeset_total_distance[k]; //delta tau(ij)
        for(int i=0;i<all_nodeset[k].size()-1;i++){
            pheromones[all_nodeset[k][i]][all_nodeset[k][i+1]] += delta;
            pheromones[all_nodeset[k][i+1]][all_nodeset[k][i]] += delta;
        }
        pheromones[all_nodeset[k][all_nodeset[k].size()-1]][all_nodeset[k][0]] += delta;
        pheromones[all_nodeset[k][0]][all_nodeset[k][all_nodeset[k].size()-1]] += delta;
    }
}

string print_to_file(vector<int> &average_min_nodeset, double average_min_distance, double global_min_distance){
    ofstream output_datafile;
    string output_file_address;
    cout << "Input the output file address or filename:" << endl;
    cin >> output_file_address;
    output_datafile.open(output_file_address);
    output_datafile << "mean distance: " << average_min_distance << endl;
    output_datafile << "distance: " << global_min_distance << endl; 
    for(int i=0;i<average_min_nodeset.size();i++){
        output_datafile << average_min_nodeset[i] + 1 << endl; 
    }
    output_datafile.close();
    return output_file_address;
}

void plot(vector<int> min_nodeset, vector<node> &nodeset,string &output_file_address){
    min_nodeset.push_back(min_nodeset[0]);
    FILE *gnuplotpipe = _popen("gnuplot -persist","w");
    if(!gnuplotpipe){
        return;
    }
    string png_name = output_file_address.substr(0, output_file_address.find_last_of(".") );

    fprintf(gnuplotpipe, "set terminal wxt\n");
    fprintf(gnuplotpipe, "set title '%s'\n",png_name.c_str());
    fprintf(gnuplotpipe, "set xrange [0:100]\n");
    fprintf(gnuplotpipe, "set yrange [0:100]\n");
    fprintf(gnuplotpipe, "unset key\n");
    fprintf(gnuplotpipe, "plot '-' lt 1 lc 1 w lp\n");
    for(int i=0;i<min_nodeset.size();i++){
        fprintf(gnuplotpipe, "%d %d\n", nodeset[min_nodeset[i]].x_coord,nodeset[min_nodeset[i]].y_coord);
    }
    fprintf(gnuplotpipe, "e\n");

    fprintf(gnuplotpipe, "set terminal pngcairo\n");
    fprintf(gnuplotpipe, "set title '%s'\n",png_name.c_str());
    fprintf(gnuplotpipe, "set output '%s.png'\n",png_name.c_str());
    fprintf(gnuplotpipe, "set xrange [0:100]\n");
    fprintf(gnuplotpipe, "set yrange [0:100]\n");
    fprintf(gnuplotpipe, "unset key\n");
    fprintf(gnuplotpipe, "plot '-' lt 1 lc 1 w lp\n");
    for(int i=0;i<min_nodeset.size();i++){
        fprintf(gnuplotpipe, "%d %d\n", nodeset[min_nodeset[i]].x_coord,nodeset[min_nodeset[i]].y_coord);
    }
    fprintf(gnuplotpipe, "%d %d\n", nodeset[min_nodeset[0]].x_coord,nodeset[min_nodeset[0]].y_coord);
    fprintf(gnuplotpipe, "e\n");

    _pclose(gnuplotpipe);
}

