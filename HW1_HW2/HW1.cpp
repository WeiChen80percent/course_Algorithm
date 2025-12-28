#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
using namespace std;

struct node{
    int city_id;
    int x_coord;
    int y_coord;
};

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

double calculate_distance(vector<node> &nodeset,double **distance_array){
    double total_distance = 0;
    for(int i=0;i<nodeset.size()-1;i++){
        total_distance += distance_array[nodeset[i].city_id-1][nodeset[i+1].city_id-1];
    }
    total_distance += distance_array[nodeset[nodeset.size()-1].city_id-1][nodeset[0].city_id-1];
    return total_distance;
}

void permutation(vector<node> &nodeset,int start_node, int len, double **distance_array,vector<node> &min_nodeset, double &min_total_distance){
    if(start_node == len){
        double total_distance = 0;
        total_distance = calculate_distance(nodeset,distance_array);
        if(total_distance<min_total_distance){
            min_total_distance = total_distance;
            min_nodeset = nodeset;
        }
        return;
    }
    for(int i=start_node;i<len;i++){
        swap(nodeset[i],nodeset[start_node]);
        permutation(nodeset,start_node+1,len,distance_array,min_nodeset,min_total_distance);
        swap(nodeset[i],nodeset[start_node]);
    }
    return;
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
    
    fprintf(gnuplotpipe, "set terminal wxt\n");
    fprintf(gnuplotpipe, "set title '%s'\n",output_file_address.c_str());
    fprintf(gnuplotpipe, "set xrange [0:100]\n");
    fprintf(gnuplotpipe, "set yrange [0:100]\n");
    fprintf(gnuplotpipe, "unset key\n");
    fprintf(gnuplotpipe, "plot '-' lt 1 lc 1 w lp\n");
    for(int i=0;i<min_nodeset.size();i++){
        fprintf(gnuplotpipe, "%d %d\n", min_nodeset[i].x_coord,min_nodeset[i].y_coord);
    }
    fprintf(gnuplotpipe, "e\n");

    fprintf(gnuplotpipe, "set terminal pngcairo\n");
    fprintf(gnuplotpipe, "set title '%s'\n",output_file_address.c_str());
    fprintf(gnuplotpipe, "set output '%s.png'\n",output_file_address.c_str());
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

int main(){
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
    double **distance_array = new double *[nodeset.size()];
    for(int i=0;i<nodeset.size();i++){
        distance_array[i] = new double[nodeset.size()];
    }
    distance_array_calculation(nodeset,distance_array);
    vector<node> min_nodeset;
    double min_total_distance= numeric_limits<double>::max();
    permutation(nodeset,0,nodeset.size(),distance_array,min_nodeset,min_total_distance);
    string output_file_address;
    output_file_address = print_to_file(min_nodeset,min_total_distance);
    plot(min_nodeset,output_file_address);
    datafile.close();
}