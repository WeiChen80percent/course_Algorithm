#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <map>
#include <iomanip>
using namespace std;

#define PI 3.14159
double MAX_DOUBLE = std::numeric_limits<double>::max();
double MIN_DOUBLE = std::numeric_limits<double>::min();

struct node{
    double x_coord;
    double y_coord;
};

struct city{
    int city_id;
    node coord;
};

void distance_array_calculation(vector<city> &nodeset,double **distance_array);
double distance_square(node &a, node &b);
node centroid(vector<city> &cityset);
double gaussian_calculation(double distance_square,double K);
void elastic_net(vector<city> &cityset, vector<node> &path_nodeset, double alpha, double beta, double K_initial, double K_final, int cool_iter_interval, double cool_rate, ofstream &path_evolve_file, double min_x, double min_y, double scale);
void get_tsp_result(vector<city> &cityset, vector<node> &path_nodeset, vector<vector<int>> &route_list, vector<double> &distance_list, double **distance_array);
string print_to_file(vector<int> &min_nodeset, double average_min_distance, double best_min_distance);
void plot_gif(string output_file_address);

int main(){
    // srand(time(NULL));

    int city_id, x_coord, y_coord;
    ifstream datafile;
    string file_address;
    cout << "Input the input file address or filename:"<<endl;
    cin >> file_address;
    datafile.open(file_address);
    vector<city> cityset;
    ofstream city_file("city.dat");
    while(datafile >> city_id >> x_coord >> y_coord){
        city n;
        n.city_id = city_id;
        n.coord.x_coord = x_coord;
        n.coord.y_coord = y_coord;
        city_file << x_coord << " " << y_coord << endl;
        cityset.push_back(n);
    } 
    city_file.close();

    double **distance_array = new double *[cityset.size()];
    for(int i=0;i<cityset.size();i++){
        distance_array[i] = new double[cityset.size()];
    }
    distance_array_calculation(cityset,distance_array);

    double min_x = MAX_DOUBLE, min_y = MAX_DOUBLE;
    double max_x = MIN_DOUBLE, max_y = MIN_DOUBLE;
    for(int i=0;i<cityset.size();i++){
        if(cityset[i].coord.x_coord < min_x){
            min_x = cityset[i].coord.x_coord;
        }
        if(cityset[i].coord.y_coord < min_y){
            min_y = cityset[i].coord.y_coord;
        }
        if(cityset[i].coord.x_coord > max_x){
            max_x = cityset[i].coord.x_coord;
        }
        if(cityset[i].coord.y_coord > max_y){
            max_y = cityset[i].coord.y_coord;
        }
    }

    double scale = max(max_x - min_x, max_y - min_y);
    if(scale == 0){
        scale = 1;
    }

    //正規化
    for(int i=0;i<cityset.size();i++){
        cityset[i].coord.x_coord = (cityset[i].coord.x_coord - min_x) / scale;
        cityset[i].coord.y_coord = (cityset[i].coord.y_coord - min_y) / scale;
    }


    //----------------------------------------------------------intitializtion----------------------------------------------------------
    node centroid_node = centroid(cityset); //calculate centroid of all cities
    int M = round(2.5 * cityset.size());
    double radius = 0.5;
    vector<node> path_nodeset;
    // Initialize net 的點
    for(int i=0;i<M;i++){
        node temp_node;
        temp_node.x_coord = centroid_node.x_coord + (radius * cos(2 * PI * i / M));
        temp_node.y_coord = centroid_node.y_coord + (radius * sin(2 * PI * i / M));
        path_nodeset.push_back(temp_node);
    }
    double alpha = 0.2, beta = 2, K_initial = 0.2, K_final = 0.05, cool_rate = 0.99;
    int cool_iter_interval = 25, run_time = 30;
    vector<vector<int>> route_list;
    vector<double> distance_list;
    double average_min_distance = 0;
    double best_min_distance = MAX_DOUBLE;
    vector<int> best_min_distance_route;

    ofstream path_evolve_file("evolution.dat");
    path_evolve_file << fixed << setprecision(6);

    //----------------------------------------------------------intitializtion----------------------------------------------------------
    //----------------------------------------------------------execution----------------------------------------------------------
    for(int i=0;i<run_time;i++){
        vector<node> path_nodeset_copy;
        path_nodeset_copy = path_nodeset;
        elastic_net(cityset, path_nodeset_copy, alpha, beta, K_initial, K_final, cool_iter_interval, cool_rate, path_evolve_file, min_x, min_y, scale);
        get_tsp_result(cityset, path_nodeset_copy,route_list, distance_list, distance_array);
        average_min_distance += distance_list[i];
        if(distance_list[i] < best_min_distance){
            best_min_distance_route = route_list[i];
            best_min_distance = distance_list[i];
        }
    }
    path_evolve_file.close();
    average_min_distance /= run_time;

    string output_file_address;
    output_file_address = print_to_file(best_min_distance_route, average_min_distance, best_min_distance);
    plot_gif(output_file_address);
    //----------------------------------------------------------execution----------------------------------------------------------
}

void distance_array_calculation(vector<city> &nodeset,double **distance_array){
    for(int i=0;i<nodeset.size();i++){
        distance_array[i][i] = 0;
        for(int j=i+1;j<nodeset.size();j++){
            double dx = nodeset[i].coord.x_coord - nodeset[j].coord.x_coord;
            double dy = nodeset[i].coord.y_coord - nodeset[j].coord.y_coord;
            distance_array[i][j] = sqrt(dx *dx + dy*dy);
            distance_array[j][i] = distance_array[i][j]; 
        }
    }
    return;
}

double distance_square(node &a, node &b){
    double dx = a.x_coord - b.x_coord;
    double dy = a.y_coord - b.y_coord;
    return (dx * dx + dy * dy);
}

node centroid(vector<city> &cityset){
    int num_of_city = cityset.size();
    double mean_x_coord=0;
    double mean_y_coord=0;
    for(int i=0;i<num_of_city;i++){
        mean_x_coord += cityset[i].coord.x_coord;
        mean_y_coord += cityset[i].coord.y_coord;
    }
    mean_x_coord /= num_of_city;
    mean_y_coord /= num_of_city;
    node n;
    n.x_coord = mean_x_coord;
    n.y_coord = mean_y_coord;
    return n;
}

double gaussian_calculation(double distance_square,double minus_two_k_square_reciprocal){
    double result = exp(distance_square * minus_two_k_square_reciprocal);
    return result;
}

void elastic_net(vector<city> &cityset, vector<node> &path_nodeset, double alpha, double beta, double K_initial, double K_final, int cool_iter_interval, double cool_rate, ofstream &path_evolve_file, double min_x, double min_y, double scale){
    double K = K_initial; // Use to degrade later
    long long int iter = 0; // calculate iteration times
    vector<node> delta_path_nodeset(path_nodeset.size(),{0,0}); 
    vector<double> gaussian_values(path_nodeset.size(),{0});
    while(iter <= 10000){
        fill(delta_path_nodeset.begin(), delta_path_nodeset.end(), node{0,0});
        // fill(gaussian_values.begin(), gaussian_values.end(), 0);
        double minus_two_k_square_reciprocal = -1 / (2 * K * K);
        for(int city=0;city<cityset.size();city++){
            double sum_of_gaussian = 0;
            for(int path_node=0;path_node<path_nodeset.size();path_node++){
                gaussian_values[path_node] = gaussian_calculation(distance_square(cityset[city].coord,path_nodeset[path_node]),minus_two_k_square_reciprocal);
                sum_of_gaussian += gaussian_values[path_node];
            }
            double sum_of_gaussian_reciprocal = 1 / sum_of_gaussian;
            for(int path_node=0;path_node<path_nodeset.size();path_node++){
                double gaussian_city_path_node = gaussian_values[path_node];
                double weight_city_path_node = gaussian_city_path_node * sum_of_gaussian_reciprocal;
                delta_path_nodeset[path_node].x_coord += weight_city_path_node * (cityset[city].coord.x_coord - path_nodeset[path_node].x_coord) * alpha;
                delta_path_nodeset[path_node].y_coord += weight_city_path_node * (cityset[city].coord.y_coord - path_nodeset[path_node].y_coord) * alpha;
            }
        }
        double beta_K =  beta * K;
        for(int path_node=0;path_node<path_nodeset.size();path_node++){
            int path_node_plus1 = (path_node + 1) % path_nodeset.size();
            int path_node_minus1 = (path_node - 1 + path_nodeset.size()) % path_nodeset.size();
            delta_path_nodeset[path_node].x_coord += beta_K * (path_nodeset[path_node_plus1].x_coord - 2 * path_nodeset[path_node].x_coord + path_nodeset[path_node_minus1].x_coord); 
            delta_path_nodeset[path_node].y_coord += beta_K * (path_nodeset[path_node_plus1].y_coord - 2 * path_nodeset[path_node].y_coord + path_nodeset[path_node_minus1].y_coord);
        }

        for(int path_node=0;path_node<path_nodeset.size();path_node++){
            path_nodeset[path_node].x_coord += delta_path_nodeset[path_node].x_coord;
            path_nodeset[path_node].y_coord += delta_path_nodeset[path_node].y_coord;
        }

        iter++;
        if((iter % cool_iter_interval) == 0){
            K *= cool_rate;
        }

        if(iter % 100 ==0){
            for(int path_node=0;path_node<path_nodeset.size();path_node++){
                double original_x = path_nodeset[path_node].x_coord * scale + min_x;
                double original_y = path_nodeset[path_node].y_coord * scale + min_y;
                path_evolve_file << original_x << " " << original_y << endl;
            }
            double original_x = path_nodeset[0].x_coord * scale + min_x;
            double original_y = path_nodeset[0].y_coord * scale + min_y;
            path_evolve_file << original_x << " " << original_y <<endl;
            path_evolve_file << endl << endl;
        }
    }
}

void get_tsp_result(vector<city> &cityset, vector<node> &path_nodeset, vector<vector<int>> &route_list, vector<double> &distance_list, double **distance_array){
    int N = cityset.size();
    int M = path_nodeset.size();

    map<int,int> path_node_to_city;
    bool *path_node_used = new bool[M];
    for(int i=0;i<M;i++){
        path_node_used[i] = false;
    }

    for(int i=0;i<N;i++){
        double min_distance = MAX_DOUBLE;
        int closest_path_node_index = -1;
        for(int j=0;j<M;j++){
            double distance = distance_square(cityset[i].coord,path_nodeset[j]);
            if(distance < min_distance && path_node_used[j] == false){
                min_distance = distance;
                closest_path_node_index = j;
            }
        }
        path_node_to_city[closest_path_node_index] = i;
        path_node_used[closest_path_node_index] = true;
    }
    vector<int> city_route;
    for(auto &pair:path_node_to_city){
        int city = pair.second;
        city_route.push_back(city);
    }

    double distance = 0;
    for(int i=0;i<N;i++){
        distance += distance_array[city_route[i]][city_route[(i+1) % N]];
    }
    
    route_list.push_back(city_route);
    distance_list.push_back(distance);
}

string print_to_file(vector<int> &min_nodeset, double average_min_distance, double best_min_distance){
    ofstream output_datafile;
    string output_file_address;
    cout << "Input the output file address or filename:" << endl;
    cin >> output_file_address;
    output_datafile.open(output_file_address);
    output_datafile << "mean distance: " << average_min_distance << endl;
    output_datafile << "distance: " << best_min_distance << endl;
    for(int i=0;i<min_nodeset.size();i++){
        output_datafile << min_nodeset[i] + 1 << endl; 
    }
    output_datafile.close();
    return output_file_address;
}

void plot_gif(string output_file_address){
    string gif_name = output_file_address.substr(0, output_file_address.find_last_of(".") )+ ".gif";
    FILE *gnuplotpipe = _popen("gnuplot -persist","w");
    if(!gnuplotpipe){
        return;
    }

    fprintf(gnuplotpipe, "set terminal gif animate delay 10 size 800,800\n");
    fprintf(gnuplotpipe, "set output '%s'\n",gif_name.c_str());
    fprintf(gnuplotpipe, "set title '%s'\n",gif_name.c_str());

    fprintf(gnuplotpipe, "set xrange [0:100]\n");
    fprintf(gnuplotpipe, "set yrange [0:100]\n");
    fprintf(gnuplotpipe, "unset key\n");
    fprintf(gnuplotpipe, "stats 'evolution.dat' nooutput\n");

    fprintf(gnuplotpipe, "do for [i=0:STATS_blocks-1]{\n");
    fprintf(gnuplotpipe, "      plot 'city.dat' with points pt 7 ps 1.5 lc rgb 'red', \\\n");
    fprintf(gnuplotpipe, "          'evolution.dat' index i with lines lw 2 lc rgb 'blue'\n");
    fprintf(gnuplotpipe, "}\n");
    _pclose(gnuplotpipe);
}
