#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <random>
#include <ctime>
#include <map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace std;

    class Interval {
    public:
        int net_index;
        pair<int, int> boundary;
        bool update;
    };

    class Net {
    public:
        pair<char, int> track;
        set<int> loc;
        pair<int, int> bounary;
        int preCount;
        int sucCount;
        vector<int> successor;
        vector<int> predecessor;
    };

    class Routing { 
    public:
        Routing(){};
        void read_file(fstream &);
        void write_file(fstream &);
        void left_edge();
        void init();

    private:
        void build_graph();
        void get_interval();
        map<int, map<int, int>> get_side_track(map<int, map<int, int>>);

        vector<int> BOT;
        vector<int> TOP;
        map<int, map<int, int>> TOP_boundary;
        map<int, map<int, int>> BOT_boundary;
        map<int, map<int, int>> TOP_track;
        map<int, map<int, int>> BOT_track;

        map<int, Net> Net_list;
        vector<Interval> I;
        int added_track;
        
    };

    bool mycompare(Interval, Interval);

int main(int argc, char *argv[]) {

    fstream myFile;
    Routing R;
    string intputFile = argv[1]; 
    string outputFile = argv[2];

    myFile.open(intputFile, ios::in);
    // myFile.open("case3.in", ios::in);
    R.read_file(myFile);
    myFile.close();

    R.init();

    R.left_edge();

    myFile.open(outputFile, ios::out);
    // myFile.open("case3.out", ios::out);
    R.write_file(myFile);
    myFile.close();
    
    return 0;
}

void Routing::write_file(fstream & myFile) {

    myFile << "Channel density: " << added_track;
    for(auto it = Net_list.begin(); it != Net_list.end(); it++) {
        myFile << endl;
        myFile << "Net " << (*it).first << endl;
        myFile << (*it).second.track.first << (*it).second.track.second << " ";
        myFile << (*it).second.bounary.first << " " << (*it).second.bounary.second;
    }
}

void Routing::left_edge() {

    sort(I.begin(), I.begin()+I.size(), mycompare);

    int ICount = I.size();

    for(auto it = (TOP_track.rbegin()); it != TOP_track.rend(); it++) {
        vector<int> rm_list;

        for(auto inner = (*it).second.begin(); inner != (*it).second.end(); inner++) {
            int left = (*inner).first;
            int right = (*inner).second;

            int watermark = left - 1;
            for(int j = 0; j < I.size(); j++) {
                if(!Net_list[I[j].net_index].preCount && !I[j].update && I[j].boundary.first > watermark && I[j].boundary.second <= right) {
                    I[j].update = true;
                    ICount--;
                    watermark = I[j].boundary.second;

                    Net_list[I[j].net_index].track.first = 'T';
                    Net_list[I[j].net_index].track.second = (*it).first;

                    rm_list.push_back(I[j].net_index);

                }
            }
        }


        for(int n = 0; n < rm_list.size(); n++) {
            for(int k = 0; k < Net_list[rm_list[n]].successor.size(); k++) {
                int sucNum = Net_list[rm_list[n]].successor[k];
                Net_list[sucNum].preCount--;
            }
        }
    }

    for(auto it = (BOT_track.rbegin()); it != BOT_track.rend(); it++) {
        vector<int> rm_list;

        for(auto inner = (*it).second.begin(); inner != (*it).second.end(); inner++) {
            int left = (*inner).first;
            int right = (*inner).second;

            int watermark = left - 1;
            for(int j = 0; j < I.size(); j++) {
                if(!Net_list[I[j].net_index].sucCount && !I[j].update && I[j].boundary.first > watermark && I[j].boundary.second <= right) {
                    I[j].update = true;
                    ICount--;
                    watermark = I[j].boundary.second;

                    Net_list[I[j].net_index].track.first = 'B';
                    Net_list[I[j].net_index].track.second = (*it).first;

                    rm_list.push_back(I[j].net_index);

                }
            }
        }


        for(int n = 0; n < rm_list.size(); n++) {
            for(int k = 0; k < Net_list[rm_list[n]].predecessor.size(); k++) {
                int sucNum = Net_list[rm_list[n]].predecessor[k];
                Net_list[sucNum].sucCount--;
            }
        }
    }

    
    while(ICount > 0) {

        int watermark = -1;
        vector<int> rm_list;
        added_track++;

        for(int i = 0; i < I.size(); i++) {
            
            if(!I[i].update && !Net_list[I[i].net_index].sucCount && I[i].boundary.first > watermark) {
                I[i].update = true;
                ICount--;
                watermark = I[i].boundary.second;

                Net_list[I[i].net_index].track.first = 'C';
                Net_list[I[i].net_index].track.second = added_track;

                rm_list.push_back(I[i].net_index);
            }

        }
        for(int n = 0; n < rm_list.size(); n++) {
            for(int k = 0; k < Net_list[rm_list[n]].predecessor.size(); k++) {
                int sucNum = Net_list[rm_list[n]].predecessor[k];
                Net_list[sucNum].sucCount--;
            }
        }
    }
    


}

void Routing::build_graph() {
    for(int i = 0; i < TOP.size(); i++) {
        if(TOP[i]) {
            Net_list[TOP[i]].loc.insert(i);
        }
        if(BOT[i]) {
            Net_list[BOT[i]].loc.insert(i);
        }
        if(TOP[i] && BOT[i] && TOP[i] != BOT[i]) {
            Net_list[TOP[i]].successor.push_back(BOT[i]);
            Net_list[TOP[i]].sucCount++;
            Net_list[BOT[i]].predecessor.push_back(TOP[i]);
            Net_list[BOT[i]].preCount++;
        }
    }
}
void Routing::get_interval() {
    for(auto it = Net_list.begin(); it != Net_list.end(); it++) {
        Interval temp_I;

        (*it).second.bounary.first = *((*it).second.loc.begin());
        (*it).second.bounary.second = *(--(*it).second.loc.end());

        temp_I.boundary = (*it).second.bounary;
        temp_I.net_index = (*it).first;
        temp_I.update = false;

        I.push_back(temp_I);
    }
}

void Routing::init() {
    build_graph();
    get_interval();
    TOP_track = get_side_track(TOP_boundary);
    BOT_track = get_side_track(BOT_boundary);
    added_track = 0;

}

map<int, map<int, int>> Routing::get_side_track(map<int, map<int, int>> boundary) {

    map<int, map<int, int>> side_track;
    int left, right;

    for(auto it = boundary.begin(); it != --boundary.end(); it++) {
        map<int, int> temp_map;
        for(auto inner = boundary.begin(); inner != next(it); inner++) {
            for(auto in = (*inner).second.begin(); in != (*inner).second.end(); in++) {
                temp_map[(*in).first] = (*in).second;
            }
        }

        for(auto inner = temp_map.begin(); inner != --temp_map.end(); inner++) {
            left = (*inner).second + 1;
            right = (*next(inner)).first - 1;
            if(right - left < 1)
                continue;

            side_track[(*it).first][left] = right;
        }

        auto tmp = temp_map.begin();
        if((*tmp).first != 0 && (*tmp).first != 1) {
            left = 0;
            right = (*tmp).first - 1;
            side_track[(*it).first][left] = right;
        }
        tmp = --temp_map.end();
        if((*tmp).second != TOP.size()-1 && (*tmp).second != TOP.size()-2) {
            left = (*tmp).second + 1;
            right = TOP.size()-1;
            side_track[(*it).first][left] = right;
        }

    }

    return side_track;
}

void Routing::read_file(fstream& myFile) {
    
    string line;
    string temp_str;
    char temp_char;
    int temp_int, temp_int1, temp_int2;
    map<int, int> temp_map;

    while(getline(myFile, line)) {
        stringstream ss(line);

        if(line[0] == 'T') {
            ss >> temp_char >> temp_int >> temp_int1 >> temp_int2;
            TOP_boundary[temp_int][temp_int1] = temp_int2;
        }
        else if(line[0] == 'B') {
            ss >> temp_char >> temp_int >> temp_int1 >> temp_int2;
            BOT_boundary[temp_int][temp_int1] = temp_int2;
        }
        else {
            vector<int> temp_vec;
            while(ss >> temp_str) {
                if(temp_str[0] == '/')
                    break;

                temp_int = stoi(temp_str);
                temp_vec.push_back(temp_int);

                if(temp_int)
                    Net_list[temp_int];
            }
            if(!TOP.size()) 
                TOP = temp_vec;
            else
                BOT = temp_vec;
        }
    }
}

bool mycompare(Interval a, Interval b){
   return a.boundary.first < b.boundary.first;
}
