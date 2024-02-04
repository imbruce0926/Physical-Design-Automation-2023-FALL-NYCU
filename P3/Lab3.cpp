#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <random>
#include <ctime>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cmath>


using namespace std;

// time(NULL)
#define seed         time(NULL)
#define iter_T       5
#define T0           30000
#define freezing     0.1

class MOS {
public:
    string name;
    string gate;
    string drain;
    string source;
    bool type;
    bool D_connected;
    bool S_connected;
    bool dir; // D -> S: 0 // S -> D: 1
};

class placement {
public:
    placement();
    void read_file(fstream &, clock_t);
    void write_file(fstream &);
    void init();
    void SA();

    double best_cost;

private:
    vector<MOS> MOS_buger;
    vector<MOS> PMOS_buger;
    vector<MOS> NMOS_buger;
    vector<int> sequence;
    vector<int> sequence_dummy;
    int MOS_Num;
    int Net_Num;

    vector<int> PMOS_dummy_pos;
    vector<int> NMOS_dummy_pos;
    set<int> dummy_pos;

    unordered_map<string, set<int>> pin_loc;
    unordered_map<string, pair<bool, bool>> pin_state; // pair<bool, bool> = pair<inPMOS, inNMOS>
    vector<string> PMOS_pin;
    vector<string> NMOS_pin;

    double PMOS_width;
    double NMOS_width;
    double cost;


    vector<int> best_sequence_dummy;
    vector<string> best_PMOS_pin;
    vector<string> best_NMOS_pin;
    

    clock_t start;

    void seperate_MOS();
    void get_sequence_with_dummy();
    vector<int> get_PN_dummy_pos(vector<MOS> &);
    set<int> merge_dummy_pos();
    vector<int> to_sequence_with_dummy();
    void SWAP(int *, int *, bool);
    void get_MOS_location();
    void get_pin();
    double HPWL();
    void clear();
    void update_cooling_factor(double *, double);
    bool is_movable(double, double);
    double calculate_cost();

};

int main(int argc, char *argv[]) {
    fstream myFile;
    string intputFile = argv[1]; 
    string outputFile = argv[2];

    placement P;

    clock_t start = clock();

    myFile.open(intputFile, ios::in);
    // myFile.open("testcase.sp", ios::in);
    // L.read_file(myFile, outputFile);
    P.read_file(myFile, start);
    myFile.close();

    P.init();
    P.SA();

    // myFile.open("testcase_output.txt", ios::out);
    myFile.open(outputFile, ios::out);
    P.write_file(myFile);
    myFile.close();

    // cout << endl << P.best_cost << endl;
    // cout << "time = " << double(clock()-start)/CLOCKS_PER_SEC << " s" << endl;
    return 0;

}

placement::placement() {
    srand(seed);
    best_cost = 1.79769e+308;
}

void placement::init() {
    for(int i = 0; i < MOS_Num; i++) 
        sequence.push_back(i+1);
    seperate_MOS();
}

void placement::clear() {

    // sequence_dummy.clear();
    // dummy_pos.clear();
    PMOS_pin.clear();
    NMOS_pin.clear();
    // for(auto it = pin_loc.begin(); it != pin_loc.end(); it++) 
    //     (*it).second.clear();
    pin_loc.clear();
    // PMOS_dummy_pos.clear();
    // NMOS_dummy_pos.clear();
    for(int i = 1; i < PMOS_buger.size(); i++) {
        PMOS_buger[i].D_connected = PMOS_buger[i].S_connected = false;
        NMOS_buger[i].D_connected = NMOS_buger[i].S_connected = false;
        PMOS_buger[i].dir = NMOS_buger[i].dir = 0;
    }
}

double placement::calculate_cost() {
    double cost;
    get_sequence_with_dummy();
    get_pin();
    cost = HPWL();
    return cost;
}

void placement::SA() {


    // clear();

    double T, w, curr_cost, new_cost, delta_cost;
    int MoveCount = 0, uphill = 0, N;
    // int reject = 0;
    int i, j;

    // vector<int> best_sequence;


    best_cost = calculate_cost();
    curr_cost = best_cost;
    new_cost = best_cost;

    // best_sequence = sequence;
    best_sequence_dummy = sequence_dummy;
    best_PMOS_pin = PMOS_pin;
    best_NMOS_pin = NMOS_pin;

    T = T0;
    // N = iter_T * sqrt(100 * 200);
    N = iter_T * sqrt(MOS_Num * Net_Num);

    // cout << MOS_Num << " " <<  Net_Num << endl;

    clear();

    vector<MOS> temp_NMOS, temp_PMOS;
    temp_NMOS = NMOS_buger;
    temp_PMOS = PMOS_buger;
    bool is_moved = true;

    do {
        
        MoveCount = 0; uphill = 0; 
        // reject = 0;
        update_cooling_factor(&w, T);
        
        do {
            // cout << "MoveCount " << MoveCount << endl;
            // curr_cost = calculate_cost();
            // clear();

            if(is_moved) {
                curr_cost = new_cost;
            }
            // else {
            //     PMOS_buger = temp_PMOS;
            //     NMOS_buger = temp_NMOS;
            // }

            // temp_PMOS = PMOS_buger;
            // temp_NMOS = NMOS_buger;

            // cout << curr_cost << endl;;


            SWAP(&i, &j, true);
            new_cost = calculate_cost();

            delta_cost = new_cost - curr_cost;
            delta_cost /= (Net_Num * MOS_Num);

            MoveCount++;
            
            is_moved = is_movable(delta_cost, T);

            if(is_moved) {
                if(delta_cost > 0)
                    uphill++;

                if(new_cost < best_cost) {
                    // cout << new_cost << " " << best_cost << endl;

                    // best_sequence = sequence;
                    best_cost = new_cost;
                    best_sequence_dummy = sequence_dummy;
                    best_PMOS_pin = PMOS_pin;
                    best_NMOS_pin = NMOS_pin;

                    
                }
            }
            else {
                SWAP(&i, &j, false);
                // reject++;
            }
            clear();

            // cout << is_moved << endl << endl;
        }
        while((uphill < N) && (MoveCount < 2 * N) /*&& (double(clock()-start)/CLOCKS_PER_SEC) < 7080*/);
        // cout << "T = " << T  << " time = " << double(clock()-start)/CLOCKS_PER_SEC << " s"<<endl;
        // cout << "MoveCount: " << MoveCount << " 2N: " << 2 * N << endl;
        // cout << "uphill: " << uphill << " N: " << N << endl << endl;

        T = w * T;
        // sequence = best_sequence;
    }
    while(/*((double)reject/MoveCount) < 0.95 && */(T > freezing) /*&& (double(clock()-start)/CLOCKS_PER_SEC) < 7080*/);
    // cout << "reject: " << reject << " MoveCount: " << MoveCount << endl;

}

void placement::get_sequence_with_dummy() {

    PMOS_dummy_pos = get_PN_dummy_pos(PMOS_buger);
    NMOS_dummy_pos = get_PN_dummy_pos(NMOS_buger);
    // get_PN_dummy_pos();
    dummy_pos = merge_dummy_pos();
    sequence_dummy = to_sequence_with_dummy();

}

vector<int> placement::get_PN_dummy_pos(vector<MOS>& buger) {
    vector<int> dummy;
    bool broken = false;

    for(int i = 0; i < sequence.size()-1; i++) {

        if(buger[sequence[i]].drain == buger[sequence[i+1]].drain && !buger[sequence[i]].D_connected) {
            buger[sequence[i]].D_connected = 1;
            buger[sequence[i]].dir = 1;
            buger[sequence[i+1]].D_connected = 1;
            buger[sequence[i+1]].dir = 0;
        }
        else if(buger[sequence[i]].drain == buger[sequence[i+1]].source && !buger[sequence[i]].D_connected) {
            buger[sequence[i]].D_connected = 1;
            buger[sequence[i]].dir = 1;
            buger[sequence[i+1]].S_connected = 1;
            buger[sequence[i+1]].dir = 1;
        }
        else if(buger[sequence[i]].source == buger[sequence[i+1]].source && !buger[sequence[i]].S_connected) {
            buger[sequence[i]].S_connected = 1;
            buger[sequence[i]].dir = 0;
            buger[sequence[i+1]].S_connected = 1;
            buger[sequence[i+1]].dir = 1;
        }
        else if(buger[sequence[i]].source == buger[sequence[i+1]].drain && !buger[sequence[i]].S_connected) {
            buger[sequence[i]].S_connected = 1;
            buger[sequence[i]].dir = 0;
            buger[sequence[i+1]].D_connected = 1;
            buger[sequence[i+1]].dir = 0;
        }
        else {
            dummy.push_back(i);
        }
    }


    return dummy;
}

set<int> placement::merge_dummy_pos() {

    set<int> pos;

    for(int i = 0; i < PMOS_dummy_pos.size(); i++)
        pos.insert(PMOS_dummy_pos[i]);
    for(int i = 0; i < NMOS_dummy_pos.size(); i++)
        pos.insert(NMOS_dummy_pos[i]);

    return pos;
}
  
vector<int> placement::to_sequence_with_dummy() {
    vector<int> seq_dummy;
    int i = 0;
    for(auto it = dummy_pos.begin(); it != dummy_pos.end(); it++) {
        for(; i <= (*it); i++) 
            seq_dummy.push_back(sequence[i]);
        seq_dummy.push_back(0);
    }
    for(; i < sequence.size(); i++) 
        seq_dummy.push_back(sequence[i]);
    
    return seq_dummy;

}

void placement::get_pin() {

    sequence_dummy.push_back(0);
    for(int i = 0; i < sequence_dummy.size()-1; i++) {
        if(sequence_dummy[i] == 0) {
            for(int i = 0; i < 3; i++) {
                PMOS_pin.push_back("Dummy");
                NMOS_pin.push_back("Dummy");
            }
            continue;
        }

        if(PMOS_buger[sequence_dummy[i]].dir == 0) {
            PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].drain);
            pin_loc[PMOS_buger[sequence_dummy[i]].drain].insert(PMOS_pin.size()-1);

            PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].gate);
        }
        else {
            PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].source);
            pin_loc[PMOS_buger[sequence_dummy[i]].source].insert(PMOS_pin.size()-1);
            
            PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].gate);
        }

        if(NMOS_buger[sequence_dummy[i]].dir == 0) {
            NMOS_pin.push_back(NMOS_buger[sequence_dummy[i]].drain);
            pin_loc[NMOS_buger[sequence_dummy[i]].drain].insert(NMOS_pin.size()-1);
            
            NMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].gate);
        }
        else {
            NMOS_pin.push_back(NMOS_buger[sequence_dummy[i]].source);
            pin_loc[NMOS_buger[sequence_dummy[i]].source].insert(NMOS_pin.size()-1);
            
            NMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].gate);
        }

        if(sequence_dummy[i+1] == 0) {
            if(PMOS_buger[sequence_dummy[i]].dir == 0) {
                PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].source);
                pin_loc[PMOS_buger[sequence_dummy[i]].source].insert(PMOS_pin.size()-1);
            }
            else {
                PMOS_pin.push_back(PMOS_buger[sequence_dummy[i]].drain);
                pin_loc[PMOS_buger[sequence_dummy[i]].drain].insert(PMOS_pin.size()-1);
            }

            if(NMOS_buger[sequence_dummy[i]].dir == 0) {
                NMOS_pin.push_back(NMOS_buger[sequence_dummy[i]].source);
                pin_loc[NMOS_buger[sequence_dummy[i]].source].insert(NMOS_pin.size()-1);
            }
            else {
                NMOS_pin.push_back(NMOS_buger[sequence_dummy[i]].drain);
                pin_loc[NMOS_buger[sequence_dummy[i]].drain].insert(NMOS_pin.size()-1);
            }
        }
    }
    sequence_dummy.pop_back();

    Net_Num = pin_loc.size();

}

double placement::HPWL() {
    double cost = 0;
    int left, right;
    const double verti = PMOS_width * 0.5 + NMOS_width * 0.5 + 27;
    int pinNum = NMOS_pin.size();

    for (auto it = pin_state.begin(); it != pin_state.end(); it++) {
        double x = 0;
        double y = 0;
        left = *(pin_loc[(*it).first].begin());
        right = *(--(pin_loc[(*it).first].end()));
        
        if((*it).second.first == true && (*it).second.second == true) { // NMOS/PMOS
            y = verti;
        }
        else {
            if(left == right)
                continue;
        }

        if(left != 0 && right != pinNum-1) 
            x = ((right-left)/2)*20 + ((right-left)/2-1)*34 + 0.5*34 + 0.5*34;
        else if(left == 0 && right == pinNum-1) 
            x = ((right-left)/2)*20 + ((right-left)/2-1)*34 + 0.5*25 + 0.5*25;
        else 
            x = ((right-left)/2)*20 + ((right-left)/2-1)*34 + 0.5*34 + 0.5*25;
        

        cost += x + y;
    }

    return cost;
}

void placement::write_file(fstream& myFile) {
    myFile << best_cost << endl;
    for(int i = 0; i < best_sequence_dummy.size(); i++) {
        myFile << PMOS_buger[best_sequence_dummy[i]].name << " ";
    }
    myFile << endl;
    for(int i = 0; i < best_PMOS_pin.size(); i++) {
        if(best_PMOS_pin[i] == "Dummy")
            i+=2;
        myFile << best_PMOS_pin[i] << " ";
    }
    myFile << endl;

    for(int i = 0; i < best_sequence_dummy.size(); i++) {
        myFile << NMOS_buger[best_sequence_dummy[i]].name << " ";
    }
    myFile << endl;
    for(int i = 0; i < best_NMOS_pin.size(); i++) {
        if(best_NMOS_pin[i] == "Dummy")
            i+=2;
        myFile << best_NMOS_pin[i] << " ";
    }
    // myFile << endl;

    // cout << endl << best_cost << endl;
    // cout << "time = " << double(clock()-start)/CLOCKS_PER_SEC << " s" << endl;

}


void placement::SWAP(int *i, int *j, bool ctrl) {

    if(ctrl){
        (*i) = rand() % sequence.size();
        do {
            (*j) = rand() % sequence.size();
        }
        while((*i) == (*j));
    }

    int temp;
    temp = sequence[(*i)];
    sequence[(*i)] = sequence[(*j)];
    sequence[(*j)] = temp;
}


void placement::read_file(fstream& myFile, clock_t start_T) {

    start = start_T;

    MOS temp_mos;
    char temp_char;
    string temp_str;
    double temp_double;

    getline(myFile, temp_str);

    while(myFile >> temp_char) {

        myFile >> temp_str;
        if(temp_str == "ENDS")
            break;

        temp_mos.name = temp_str;

        myFile >> temp_str;
        temp_mos.drain = temp_str;

        myFile >> temp_str;
        temp_mos.gate = temp_str;

        myFile >> temp_str;
        temp_mos.source = temp_str;

        myFile >> temp_str >> temp_str;
        if(temp_str == "nmos_rvt") {
            temp_mos.type = 0;
            myFile >> temp_char >> temp_char >> temp_double;
            NMOS_width = temp_double;
            pin_state[temp_mos.source].second = true; // pair<bool, bool> = pair<inPMOS, inNMOS>
            pin_state[temp_mos.drain].second = true;
        }
        else {
            temp_mos.type = 1;
            myFile >> temp_char >> temp_char >> temp_double;
            PMOS_width = temp_double;
            pin_state[temp_mos.source].first = true;
            pin_state[temp_mos.drain].first = true;
        }

        myFile >> temp_char >> temp_str >> temp_str;

        temp_mos.D_connected = temp_mos.S_connected = temp_mos.dir = 0;        

        MOS_buger.push_back(temp_mos);
    }

    MOS_Num = MOS_buger.size() / 2;
}

void placement::seperate_MOS() {
    MOS temp_mos;

    temp_mos.name = "Dummy";
    temp_mos.drain = temp_mos.gate = temp_mos.source = "";
    temp_mos.type = 0;
    temp_mos.dir = 0;
    temp_mos.D_connected = temp_mos.S_connected = 1;
    PMOS_buger.push_back(temp_mos); // index 0 is dummy
    NMOS_buger.push_back(temp_mos); // index 0 is dummy

    vector<bool> is_visited(MOS_buger.size(), false);
    for(int i = 0; i < MOS_buger.size(); i++) {
        if(is_visited[i])
            continue;

        int type = MOS_buger[i].type;
        string gate = MOS_buger[i].gate;

        if(type == 0) 
            NMOS_buger.push_back(MOS_buger[i]);
        else
            PMOS_buger.push_back(MOS_buger[i]);

        is_visited[i] = true;

        for(int j = i+1; j < MOS_buger.size(); j++) {
            if(is_visited[j])
                continue;
            
            if(MOS_buger[j].gate == gate) {
                if(type == 0 && MOS_buger[j].type == 1) {
                    PMOS_buger.push_back(MOS_buger[j]);
                    is_visited[j] = true;
                    break;
                }
                if(type == 1 && MOS_buger[j].type == 0) {
                    NMOS_buger.push_back(MOS_buger[j]);
                    is_visited[j] = true;
                    break;
                }
                
            }

        }   
    }
}

void placement::update_cooling_factor(double *w, double T) {

    if(T <= T0 && T > 1000)
        *w = 0.6;
    else if(T <= 1000 && T > 100)
        *w = 0.7; 
    else if(T <= 100 && T > 50)
        *w = 0.85; 
    else if(T <= 100 && T > 10)
        *w = 0.9; 
    else if(T <= 10 && T > 1)
        *w = 0.95;
    else if(T <= 1 && T > freezing)
        *w = 0.98;

}

bool placement::is_movable(double delta_cost, double T) {
    if(delta_cost <= 0)
        return true;
    
    double x = (double) rand() / (RAND_MAX + 1.0);

    // cout << exp(-(delta_cost/T)) << endl;
    if(x < exp(-(delta_cost/T)))
        return true;

    return false;
}