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

double calculate_distance(vector<node> &nodeset){
    double total_distance = 0;
    for(int i=0;i<nodeset.size();i++){
        total_distance += sqrt(pow((nodeset[i%nodeset.size()].x_coord - nodeset[(i+1)%nodeset.size()].x_coord),2) + pow((nodeset[i%nodeset.size()].y_coord - nodeset[(i+1)%nodeset.size()].y_coord),2));
    }
    return total_distance;
}

void greedy(vector<node> &nodeset,int start_node, int len, vector<node> &global_min_nodeset, double &min_total_distance){
   for(int i=0;i<nodeset.size();i++){
        vector<node> copy_nodeset = nodeset;
        vector<node> local_min_nodeset;
        local_min_nodeset.push_back(copy_nodeset[i]);
        copy_nodeset.erase(copy_nodeset.begin()+i);
        while(!copy_nodeset.empty()){
            int min_index;
            double min_value = numeric_limits<double>::max();
            for(int j=0;j<copy_nodeset.size();j++){
                if(sqrt(pow((local_min_nodeset.back().x_coord - copy_nodeset[j].x_coord),2) + pow((local_min_nodeset.back().y_coord - copy_nodeset[j].y_coord),2)) < min_value){
                    min_value = sqrt(pow((local_min_nodeset.back().x_coord - copy_nodeset[j].x_coord),2) + pow((local_min_nodeset.back().y_coord - copy_nodeset[j].y_coord),2));
                    min_index = j;
                }
            }
            local_min_nodeset.push_back(copy_nodeset[min_index]);
            copy_nodeset.erase(copy_nodeset.begin() + min_index);
        }
        double cal_total_distance = calculate_distance(local_min_nodeset);
        if(cal_total_distance < min_total_distance){
            min_total_distance = cal_total_distance;
            global_min_nodeset = local_min_nodeset;
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
    // fprintf(gnuplotpipe, "set output '%s.png'\n",output_file_address.c_str());
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
    vector<node> min_nodeset;
    double min_total_distance= numeric_limits<double>::max();
    greedy(nodeset,0,nodeset.size(),min_nodeset,min_total_distance);
    string output_file_address;
    output_file_address = print_to_file(min_nodeset,min_total_distance);
    plot(min_nodeset,output_file_address);
    datafile.close();
}