#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace std;

struct node{
    int city_id;
    int x_coord;
    int y_coord;
};

void print_out(vector<node> &nodeset, double total_distance){
    double min_total_distance=INT32_MAX;
    for(int i=0;i<nodeset.size();i++){
        cout<<nodeset[i].city_id<<" ";
    }
    cout<<"distance: "<<total_distance<<endl;
}

void print_to_file(vector<node> &min_nodeset, double min_total_distance){
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
        //print_out(nodeset, total_distance); // print it on terminal
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
    print_to_file(min_nodeset,min_total_distance);
    datafile.close();
}