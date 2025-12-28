#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstdio>
using namespace std;

double MAX_DOUBLE = std::numeric_limits<double>::max();

struct node{
    int city_id;
    int x_coord;
    int y_coord;
};

string print_to_file(vector<node> &min_nodeset, double min_total_distance);
void distance_array_calculation(vector<node> &nodeset,double **distance_array);
void plot(vector<node> min_nodeset,string &output_file_address);
double dp(int current_city, int record, double **distance_array,double **dp_table,int num_of_city,int **path_table);
void path_reconstruct(int **path_table,int num_of_city,vector<node> &nodeset,vector<node> &min_nodeset);

int main(){
    int city_id, x_coord, y_coord;
    ifstream datafile;
    string file_address;
    cout << "Input the input file address or filename:"<<endl;
    cin >> file_address;
    datafile.open(file_address);
    vector<node> nodeset,min_nodeset;
    while(datafile >> city_id >> x_coord >> y_coord){
        node n;
        n.city_id = city_id;
        n.x_coord = x_coord;
        n.y_coord = y_coord;
        nodeset.push_back(n);
    } 

    //calculate distance between two cities
    double **distance_array = new double *[nodeset.size()];
    for(int i=0;i<nodeset.size();i++){
        distance_array[i] = new double[nodeset.size()];
    }
    distance_array_calculation(nodeset,distance_array);
    
    //create table for dp dp_table[start_node => nodeset.size()][record => 2^(nodeset.size()-1) not including first point]
    double **dp_table = new double *[nodeset.size()];
    int max_record_value = pow(2,nodeset.size()-1);
    for(int i=0;i<nodeset.size();i++){
        dp_table[i] = new double[max_record_value];
    }
    for(int i=0;i<nodeset.size();i++){
        for(int j=0;j<max_record_value;j++){
            dp_table[i][j] = -1;
        }
    }
   
    //path table use to save next city for each path_table[current city][record] 
    int **path_table = new int *[nodeset.size()];
    for(int i=0;i<nodeset.size();i++){
        path_table[i] = new int[max_record_value];
    }
    // start from city1(index 0) and record is 0
    double min_distance = dp(0,0,distance_array,dp_table,nodeset.size(),path_table);
    path_reconstruct(path_table,nodeset.size(),nodeset,min_nodeset);

    string output_file_address = print_to_file(min_nodeset, min_distance);
    plot(min_nodeset, output_file_address);
    datafile.close();
}

string print_to_file(vector<node> &min_nodeset, double min_total_distance){
    ofstream output_datafile;
    string output_file_address;
    cout << "Input the output file address or filename:" << endl;
    cin >> output_file_address;
    output_datafile.open(output_file_address);
    output_datafile << "Distance: " << min_total_distance <<endl;
    for(int i=0;i<min_nodeset.size();i++){
        output_datafile << min_nodeset[i].city_id << endl; 
    }
    output_datafile.close();
    return output_file_address;
}

void distance_array_calculation(vector<node> &nodeset,double **distance_array){
    for(int i=0;i<nodeset.size();i++){
        for(int j=i+1;j<nodeset.size();j++){
            distance_array[i][j] = sqrt(pow((nodeset[i].x_coord - nodeset[j].x_coord),2) + pow((nodeset[i].y_coord - nodeset[j].y_coord),2));
            distance_array[j][i] = distance_array[i][j]; 
        }
    }
    return;
}

void plot(vector<node> min_nodeset,string &output_file_address){
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
        fprintf(gnuplotpipe, "%d %d\n", min_nodeset[i].x_coord,min_nodeset[i].y_coord);
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
        fprintf(gnuplotpipe, "%d %d\n", min_nodeset[i].x_coord,min_nodeset[i].y_coord);
    }
    fprintf(gnuplotpipe, "e\n");

    _pclose(gnuplotpipe);
}

double dp(int current_city, int record, double **distance_array, double **dp_table,int num_of_city,int **path_table){
    
    //If all cities are visited
    int max_record_value = pow(2,num_of_city-1)-1;

    // Already raveled to all cities so directly return current city to start city
    if(record == max_record_value){
        return distance_array[0][current_city];
    }

    // If it is not empty so return
    if(dp_table[current_city][record] != -1){
        return dp_table[current_city][record];
    }

    double min_distance=MAX_DOUBLE;
    int minnode_index = -1;
    //city1 is skip, i is next city
    for(int i=1;i<num_of_city;i++){
        //If bit i equal to 0 => city i+1 not visited
        if(!(record & (1<<(i-1)))){
            int new_record = record | (1<<(i-1)); //change not visited city to visited  
            double distance = distance_array[current_city][i] + dp(i,new_record,distance_array,dp_table,num_of_city,path_table);  
            if(distance<min_distance){
                min_distance = distance;
                minnode_index = i;
            }
        }
    }
    dp_table[current_city][record] = min_distance;
    path_table[current_city][record] = minnode_index; //route next city
    return min_distance;
}

void path_reconstruct(int **path_table,int num_of_city,vector<node> &nodeset,vector<node> &min_nodeset){
    int max_record = pow(2,num_of_city-1)-1;
    
    int current_city=0,record=0;
    
    min_nodeset.push_back(nodeset[0]);

    // To see if answer is countructed(all city visited)
    while(record < max_record){
        int next_city = path_table[current_city][record];
        min_nodeset.push_back(nodeset[next_city]);

        //update for next iteration
        current_city = next_city;
        record = record | (1 << (next_city - 1)); //bit for the city is set to 1(visited)
    }
}