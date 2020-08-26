#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

// for string delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}


int main() {
    string flow_id;
    int time = 0;

    fstream ddos_file, output_file;
    // ddos_file.open("unbalanced_20_80_dataset.csv", ios::in);
    // output_file.open("ddos.dat", ios::out);
    ddos_file.open("final_dataset.csv", ios::in);
    output_file.open("ddos_2.dat", ios::out);
    int max_length = 0;
    bool first = true;
    string line;
    vector<string> labels, values;
    map<string, bool> if_ddos;
    if (ddos_file.is_open()){   //checking whether the file is open
        while (getline(ddos_file, line)) { //read data from file object and put it into string.
            if (first) {
                first = false;
                continue;
            }
            labels = split(line, ",");
            flow_id = labels[2];
            if (flow_id.length() > max_length)
                max_length = flow_id.length();
            bool if_ddos_ = (labels[84] == "ddos");
            output_file << flow_id << ' ' << time << ' ' << if_ddos_ << endl;
            
            /*
            if (if_ddos.find(flow_id) != if_ddos.end()) 
                if (if_ddos[flow_id] != if_ddos_)
                    cout << "not the same"  << endl;
            if_ddos[flow_id] = if_ddos_;
            */
            time += 1;
            if (time % 100000 == 0)
                cout << flow_id << ' ' << max_length << endl;
        }
        cout << time << endl;
        ddos_file.close(); //close the file object.
        output_file.close();
    }
}