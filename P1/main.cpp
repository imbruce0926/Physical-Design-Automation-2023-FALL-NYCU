#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <time.h> 
#include <ctime>

using namespace std;



class Cell_table {
public:
    vector<int> net_loc;
    bool valid;
};

class Net_table {
public:
    vector<int> cell;
};

class Cell_info {
public:
    int gain;
    bool A_B;
};

class Net_info {
public:
    int A_count;
    int B_count;
};

class node {
public:
    // int gain;
    int index;
    bool lock;
    node *last, *next;
};

class FM {
public:
    FM();
    void parser(fstream &, string);
    void partitioning(fstream &, string);
private:
    double ratio;
    vector<Cell_table> cell_table;
    vector<Net_table> net_table;
    vector<Cell_info> cell_info;
    vector<Net_info> net_info;
    vector<node *> node_array;
    vector<node> A_bucket, B_bucket;
    
    int A_size, B_size;
    double min_size, max_size;
    int cellCount; 
    int totalPinCount;
    int Pmax;
    int zero_gain;

    vector<int> moved_cell;
    vector<int> moved_gain;
    vector<int> accu_gain;
    int step;

    vector<Cell_info> cell_info_op;
    vector<Net_info> net_info_op;
    int A_size_op;
    int B_size_op;
    int last_accuGain;
    clock_t start;
    int cnt;


    void initial();
    void initial_partition();
    void initial_gain();
    void construct_gain_bucket();
    int cell_moving();
    void free_cell_A_to_B(int, int);
    void free_cell_B_to_A(int, int);
    void update_cell_gain(int);
    void update_neighbor_A(int);
    void update_neighbor_B(int);
    void print_bucket();
    void output_file(fstream &, string);

};

int main(int argc, char *argv[]) {

    fstream myFile;
    string intputFile = argv[1]; 
    FM F;

    myFile.open(intputFile, ios::in);
    F.parser(myFile, intputFile);
    myFile.close();

    string outputFile = argv[2]; 
    F.partitioning(myFile, outputFile);

    // cout << "End main";
    return 0;
}

FM::FM() {
    ratio = 0;
    cellCount = 0; 
    totalPinCount = 0;
    Pmax = 0;
    last_accuGain = -1;
    start = clock();
    cnt = 0;
}

void FM::parser(fstream& myFile, string intputFile) {

    
    myFile >> ratio;

    Cell_table temp_c;
    temp_c.valid = 0;

    Net_table temp_n;

    net_table.push_back(temp_n);

    char in;
    int cell_index, net_index, max_cell_index = 0;
    set<int> cell_in_net;

    while(myFile >> in) {
        
        if(in == 'N') {
            net_table.push_back(temp_n);
            myFile >> in >> in >> in >> net_index;
        }
        else if(in == 'c') {
            myFile >> cell_index;
            cell_in_net.insert(cell_index);
        }
        else if(in == ';') {
            for(auto iter = cell_in_net.begin(); iter != cell_in_net.end(); iter++) {
                net_table.back().cell.push_back(*iter);
                if(*iter > max_cell_index) 
                    max_cell_index = *iter;
        
                while(max_cell_index+1 > cell_table.size()) 
                    cell_table.push_back(temp_c);

                cell_table[*iter].valid = 1;
                cell_table[*iter].net_loc.push_back(net_table.size()-1);
            }
            cell_in_net.clear();
        }
    }


    for(int i = 1; i < cell_table.size(); i++) {
        if(cell_table[i].valid) {
            if(cell_table[i].net_loc.size() > Pmax) 
                Pmax = cell_table[i].net_loc.size();
        }
    }
    for(int i = 1; i < cell_table.size(); i++) {
        if(cell_table[i].valid)
            cellCount++;
    }

    // cout << cellCount << " " << net_table.size()-1 << endl;

}



void FM::partitioning(fstream& myFile, string outputFile) {

    initial();
    initial_partition();
    // print_bucket();
    
    int cnt = 0;

    while(1){

        if(cell_moving() == -1) {
            output_file(myFile, outputFile);
            
            break;
        }
        output_file(myFile, outputFile);
    }


}

void FM::initial() {
    Cell_info temp_c;
    temp_c.A_B = 0;
    temp_c.gain = 0;

    Net_info temp_n;
    temp_n.A_count = 0;
    temp_n.B_count = 0;
    
    cell_info.resize(cell_table.size(), temp_c);
    net_info.resize(net_table.size(), temp_n);

    A_size = cellCount;
    B_size = 0;

    // cout << "cellCount is " << cellCount << endl;
    min_size = (1-ratio)*cellCount/2;
    max_size = (1+ratio)*cellCount/2;

    // cout << min_size << " " << max_size << endl;

    node_array.resize(cell_table.size());
    for(int i = 1; i < node_array.size(); i++) {
        node_array[i] = new node();
        node_array[i]->index = i;
        node_array[i]->lock = 0;
        node_array[i]->last = nullptr;
        node_array[i]->next = nullptr;
    }

    node temp_node;
    temp_node.index = 0;
    temp_node.last = nullptr;
    temp_node.next = nullptr;
    temp_node.lock = 1;

    A_bucket.resize(2*Pmax+1, temp_node);
    B_bucket.resize(2*Pmax+1, temp_node);

    zero_gain = (2*Pmax+1) / 2;


}

void FM::initial_partition() {

    for(int i = 1; i < cell_table.size() && A_size > min_size+1 ; i++) {
        if(cell_table[i].valid) {
            cell_info[i].A_B = 1;
            A_size--;
            B_size++;
        }
    }

    for(int i = 1; i < net_table.size(); i++) {
        for(int j = 0; j < net_table[i].cell.size(); j++) {
            int cell_index = net_table[i].cell[j];
            if(cell_info[cell_index].A_B == 0)
                net_info[i].A_count++;
            else
                net_info[i].B_count++;
        }
    }

    int cut_size = 0;
    for(int i = 1; i < net_info.size(); i++) {
        if(net_info[i].A_count != 0 && net_info[i].B_count != 0)
            cut_size++;
    }


    // cout << "Cutsize = " << cut_size << endl;


}

int FM::cell_moving() {
    double START = clock();
    initial_gain();
    construct_gain_bucket();
    // print_bucket();

    moved_cell.clear();
    moved_gain.clear();
    accu_gain.clear();
    accu_gain.push_back(0);

    /* unlock all node */
    for(int i = 1; i < node_array.size(); i++) 
        node_array[i]->lock = 0;

    cell_info_op = cell_info;
    net_info_op = net_info;
    A_size_op = A_size;
    B_size_op = B_size;

    cnt++;

    
    for(int i = 0; i < cellCount; i++) {
        int A_max, B_max;
        bool A_found = 0, B_found = 0;
        int selected;

        if((double(clock()-start)/CLOCKS_PER_SEC) > 300)
            return -1;

        for(int i = 2*Pmax; i >= 0; i--) {
            if(A_bucket[i].next == nullptr)
                continue;
            else {
                A_found = 1;
                A_max = i;
                break;
            }
        }

        for(int i = 2*Pmax; i >= 0; i--) {
            if(B_bucket[i].next == nullptr)
                continue;
            else {
                B_found = 1;
                B_max = i;
                break;
            }
        }
        
        if(A_found && B_found) {
            if(A_max >= B_max) {
                if((A_size_op-1) > min_size){
                    selected = A_bucket[A_max].next->index;
                    free_cell_A_to_B(selected, A_max);
                    moved_gain.push_back(A_max-zero_gain);
                }
                else {
                    selected = B_bucket[B_max].next->index;
                    free_cell_B_to_A(selected, B_max);
                    moved_gain.push_back(B_max-zero_gain);
                }
            }
            else {
                if((B_size_op-1) > min_size){
                    selected = B_bucket[B_max].next->index;
                    free_cell_B_to_A(selected, B_max);
                    moved_gain.push_back(B_max-zero_gain);
                }
                else {
                    selected = A_bucket[A_max].next->index;
                    free_cell_A_to_B(selected, A_max);
                    moved_gain.push_back(A_max-zero_gain);
                }
            }
        }
        else if(A_found) {
            selected = A_bucket[A_max].next->index;
            free_cell_A_to_B(selected, A_max);
            moved_gain.push_back(A_max-zero_gain);
        }
        else if(B_found) {
            selected = B_bucket[B_max].next->index;
            free_cell_B_to_A(selected, B_max);
            moved_gain.push_back(B_max-zero_gain);
        }
        else {
            // cout << "HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
            break;
        }
        

        moved_cell.push_back(selected);

        node_array[selected]->lock = 1;

        update_cell_gain(selected);
        // print_bucket();
        
        

    }



    for(int i = 0; i < moved_gain.size(); i++) 
        accu_gain.push_back(moved_gain[i] + accu_gain[i]);
    
    int max_accu_gain = -2147483647;

    for(int i = 1; i < accu_gain.size(); i++) {
        if(accu_gain[i] >= max_accu_gain) {
            max_accu_gain = accu_gain[i];
            step = i - 1;
        }
    }

    // cout << "max_accu_gain is " << max_accu_gain << endl; 

    long long int complexity = cellCount * (net_table.size() -1);
    if(last_accuGain == 0 && max_accu_gain == 0 ) 
        return -1;

    last_accuGain = max_accu_gain;
    // if(max_accu_gain <= 0) 
    //     return -1;

    /* update to new state */

    for(int i = 0; i <= step; i++) {

        cell_info[moved_cell[i]].A_B = !cell_info[moved_cell[i]].A_B;

        for(int j = 0; j < cell_table[moved_cell[i]].net_loc.size(); j++) {
            int net_loc = cell_table[moved_cell[i]].net_loc[j];
            if(cell_info[moved_cell[i]].A_B == 0) {
                net_info[net_loc].A_count++;
                net_info[net_loc].B_count--;
            }
            else {
                net_info[net_loc].A_count--;
                net_info[net_loc].B_count++;
            }
        }

        if(cell_info[moved_cell[i]].A_B == 0) {
            A_size++;
            B_size--;
        }
        else {
            A_size--;
            B_size++;
        }
    }

    // cout << "iter " << cnt;
    if(complexity > 150000000000 && cnt == 4) {
        // cout << complexity << endl;
        return -1;
    }
    
    // cout << "cell_moving takes " << clock() - START << endl;
    return 0;

}

void FM::free_cell_A_to_B(int selected, int A_max) {
    if(node_array[selected]->next)
        node_array[selected]->next->last = node_array[selected]->last;
    A_bucket[A_max].next = node_array[selected]->next;
    node_array[selected]->last = nullptr;
    node_array[selected]->next = nullptr;
}

void FM::free_cell_B_to_A(int selected, int B_max) {
    if(node_array[selected]->next)
        node_array[selected]->next->last = node_array[selected]->last;
    B_bucket[B_max].next = node_array[selected]->next;
    node_array[selected]->last = nullptr;
    node_array[selected]->next = nullptr;
}


void FM::initial_gain() {


    for(int i = 1; i < cell_table.size(); i++) {
        if(cell_table[i].valid) {
            cell_info[i].gain = 0;
        }
    }

    for(int i = 1; i < cell_table.size(); i++) {
        if(cell_table[i].valid) {
            for(int j = 0; j < cell_table[i].net_loc.size(); j++) {
                int net_loc = cell_table[i].net_loc[j];
                /* From set == 1*/
                if((cell_info[i].A_B == 0 && net_info[net_loc].A_count == 1) || (cell_info[i].A_B == 1 && net_info[net_loc].B_count == 1))
                    cell_info[i].gain++;
                /* To set == 0 */
                if((cell_info[i].A_B == 0 && net_info[net_loc].B_count == 0) || (cell_info[i].A_B == 1 && net_info[net_loc].A_count == 0))
                    cell_info[i].gain--;
            }
        }
    }

}

void FM::construct_gain_bucket() {
    for(int i = 0; i < cell_table.size(); i++) {
        if(cell_table[i].valid) {
            
            int gain = cell_info[i].gain;
            int set = cell_info[i].A_B;
            if(set == 0) {
                node_array[i]->next = A_bucket[gain + zero_gain].next;
                node_array[i]->last = &A_bucket[gain + zero_gain];
            }
            else {
                node_array[i]->next = B_bucket[gain + zero_gain].next;
                node_array[i]->last = &B_bucket[gain + zero_gain];
            }
            node_array[i]->last->next = node_array[i];
            if(node_array[i]->next)
                node_array[i]->next->last = node_array[i];
            
        }
    }
}

void FM::update_cell_gain(int selected) {
    // cout << "*************************HERE IS UPDATING GAIN *************************" << endl;
    for(int i = 0; i < cell_table[selected].net_loc.size(); i++) {
        int from, to;
        int net_loc = cell_table[selected].net_loc[i];
        
        from = cell_info_op[selected].A_B ? net_info_op[net_loc].B_count : net_info_op[net_loc].A_count;
        to = cell_info_op[selected].A_B ? net_info_op[net_loc].A_count : net_info_op[net_loc].B_count;

        if((to == 0 || to == 1)) {
            for(int j = 0; j < net_table[net_loc].cell.size(); j++) {
                int cell_index = net_table[net_loc].cell[j];
                if(node_array[cell_index]->lock == 1)
                    continue;

                if(to == 0) 
                    cell_info_op[cell_index].gain++;
                else if(to == 1) {
                    if((cell_info_op[selected].A_B == !cell_info[cell_index].A_B)) {
                        cell_info_op[cell_index].gain--;
                        break;
                    }
                }

            }
        }

        from = from - 1;
        to = to + 1;
        
        if((from == 0 || from == 1)) {
            for(int j = 0; j < net_table[net_loc].cell.size(); j++) {
                int cell_index = net_table[net_loc].cell[j];

                if(node_array[cell_index]->lock == 1)
                    continue;

                if(from == 0) 
                    cell_info_op[cell_index].gain--;
                else if(from == 1) {
                    if((cell_info_op[selected].A_B == cell_info[cell_index].A_B)) {
                        cell_info_op[cell_index].gain++;
                        break;
                    }
                }

            }
        }
        
        vector<bool> is_check(cell_table.size(), false);

        for(int j = 0; j < net_table[net_loc].cell.size(); j++) {
            int cell_index = net_table[net_loc].cell[j];

            if(node_array[cell_index]->lock == 1 || is_check[cell_index])
                continue;

            if(cell_info_op[cell_index].A_B == 0) 
                update_neighbor_A(cell_index);
            else 
                update_neighbor_B(cell_index); 
            
            is_check[cell_index] = true;    
             
        }

        net_info_op[net_loc].A_count = cell_info_op[selected].A_B ? to : from;
        net_info_op[net_loc].B_count = cell_info_op[selected].A_B ? from : to;

    }

    A_size_op = cell_info_op[selected].A_B ? A_size_op + 1 : A_size_op - 1;
    B_size_op = cell_info_op[selected].A_B ? B_size_op - 1 : B_size_op + 1;

    cell_info_op[selected].A_B = !cell_info_op[selected].A_B;

    // print_bucket();

}

void FM::update_neighbor_A(int cell_index) {
    if(node_array[cell_index]->next)
        node_array[cell_index]->next->last = node_array[cell_index]->last;
    node_array[cell_index]->last->next = node_array[cell_index]->next;

    int gain = cell_info_op[cell_index].gain;
    node_array[cell_index]->next = A_bucket[gain+zero_gain].next;
    node_array[cell_index]->last = &A_bucket[gain+zero_gain];
    if(A_bucket[gain+zero_gain].next)
        A_bucket[gain+zero_gain].next->last = node_array[cell_index];
    A_bucket[gain+zero_gain].next = node_array[cell_index];
}

void FM::update_neighbor_B(int cell_index) {
    if(node_array[cell_index]->next)
        node_array[cell_index]->next->last = node_array[cell_index]->last;
    node_array[cell_index]->last->next = node_array[cell_index]->next;

    int gain = cell_info_op[cell_index].gain;
    node_array[cell_index]->next = B_bucket[gain+zero_gain].next;
    node_array[cell_index]->last = &B_bucket[gain+zero_gain];

    if(B_bucket[gain+zero_gain].next) { 
        B_bucket[gain+zero_gain].next->last = node_array[cell_index];
    }
    B_bucket[gain+zero_gain].next = node_array[cell_index]; 
}

void FM::output_file(fstream& myFile, string outputFile) {
    myFile.open(outputFile, ios::out);

    int cut_size = 0;
    for(int i = 0; i < net_info.size(); i++) {
        if(net_info[i].A_count != 0 && net_info[i].B_count != 0)
            cut_size++;
    }

    myFile << "Cutsize = " << cut_size << endl;
    myFile << "G1 " << A_size << endl;
    for(int i = 0; i < cell_info.size(); i++) {
        if(cell_table[i].valid) {
            if(cell_info[i].A_B == 0)
                myFile << "c" << i << " ";
        }
    }
    myFile << ";" << endl;
    myFile << "G2 " << B_size << endl;
    int count = 0;
    for(int i = 0; i < cell_info.size(); i++) {
        if(cell_table[i].valid) {
            if(cell_info[i].A_B == 1) {
                myFile << "c" << i << " ";
                count++;
            }
        }
    }
    myFile << ";";
    // cout << "Cutsize = " << cut_size << endl;

    myFile.close();
}

void FM::print_bucket() {
    cout << "***************** Set A *****************" << endl;
    for(int i = 2*Pmax; i >= 0; i--) {
        cout << "gain is " << i - zero_gain << endl;
            for(node *temp = A_bucket[i].next; temp!=nullptr; temp = temp->next)
                cout << temp->index << " ";
        cout << endl;
    }
    cout << endl;
    cout << "***************** Set B *****************" << endl;
    for(int i = 2*Pmax; i >= 0; i--) {
        cout << "gain is " << i - zero_gain << endl;
            for(node *temp = B_bucket[i].next; temp!=nullptr; temp = temp->next)
                cout << temp->index << " ";
        cout << endl;
    }
}