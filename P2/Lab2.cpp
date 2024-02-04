#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stack>
#include <queue>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <random>

using namespace std;


// time(NULL)
#define seed         time(NULL)
#define iter_T       1.5 // 1.5
#define T0           30000
#define iter_temp    50
#define freezing     0.005
#define init_P       0.99
#define t0           -1
#define cost_f       100
#define alpha        0.7

    class Size {
    public:
        double width, height;
    };
    
    class Module {
    public:
        string name;
        int x, y;
        Size aspect; 
    };
    
    class Block {
    public:
        int Lindex, Rindex;
        Size aspect;
    };

    class Node {
    public:
        class Node *left, *right, *parent;
        string data;
        int x, y;
        bool update;
        vector<Block> blks;
    };

    class Floorplanning { 
    public:
        Floorplanning();
        void read_file(fstream &, string, clock_t);
        vector<string> SA(vector<string> &);
        vector<string> init();

        bool is_success;

    private:
        int moduleCount;
        double total_area;
        Node *root, *avail;
        vector<Size> input_blk_aspect;
        vector<Module> module;
        vector<string> SPE;
        string outputname;
        fstream myFile;
        double R_lb, R_ub;
        double R_ub_temp;
        clock_t start;
        int chip_area;
        double chip_ratio;


        void print_tree_post(Node *);
        vector<string> init_SPE(vector<string> &);
        Node * creat_node(string, bool);
        void choose_possible_size(int, Node *);
        Node * SPE_to_tree(vector<string>);
        void bottom_up(Node *);
        bool is_equal(double, double);
        void M1(bool, Node **, Node **, int *, int *);
        void M2(bool, Node **, int *, int *);
        Node * M3(int *, int *);

        bool legal_check(int);
        int find_int(int);
        int find_VH(int);
        void find_int_VH(int, int *, int *);
        void find_node_ij(Node **, Node **, int, int);
        void find_node_i(Node **, int, int);
        void swap(int, int);

        double init_temp();
        double calculate_cost(Node *);
        int select_move();
        void delete_tree(Node *);

        void update_cooling_factor(double *, double);
        bool is_movable(double, double);

        void top_down(Node *, int);
        void write_file(fstream &);

    };

    bool mycompare(Block, Block);
    bool mycompare2(Block, Block);

int main(int argc, char *argv[]) {
    fstream myFile;
    Floorplanning F;
    string intputFile = argv[1]; 
    string outputFile = argv[2];


    clock_t start = clock();

    myFile.open(intputFile, ios::in);
    F.read_file(myFile, outputFile, start);
    myFile.close();

    vector<string> spe = F.init();


    while(!F.is_success) 
        spe = F.SA(spe);

        
    // cout << "time = " << double(clock()-start)/CLOCKS_PER_SEC << " s" << endl;
}

Floorplanning::Floorplanning() {
    avail = NULL;
    root = NULL;
    total_area = 0;
    moduleCount = 0;
    is_success = false;
    srand(seed);
}

vector<string> Floorplanning::SA(vector<string> &spe) {

    vector<string> best_SPE;
    Node *temp_root = NULL, *iptr = NULL, *jptr = NULL;
    double T, w, curr_cost, new_cost, delta_cost, best_cost = 1.79769e+308;
    int Move_type = 0, MoveCount = 0, uphill = 0, N, reject = 0, i, j;

    SPE = spe;

    N = iter_T * moduleCount;
    
    // T0 = 30000;
    T = T0;

    SPE = spe;

    root = SPE_to_tree(SPE);
    bottom_up(root);

    best_cost = calculate_cost(root);
    best_SPE = SPE;

    do {

        MoveCount = 0; uphill = 0; reject = 0;
        update_cooling_factor(&w, T);

        do {
            curr_cost = calculate_cost(root);
            Move_type = select_move();    

            switch(Move_type) {
                case 1:
                    M1(true, &iptr, &jptr, &i, &j);
                    new_cost = calculate_cost(root);
                    break;
                case 2:
                    M2(true, &iptr, &i, &j);
                    new_cost = calculate_cost(root);
                    break;
                case 3:
                    temp_root = M3(&i, &j);
                    if(temp_root) 
                        new_cost = calculate_cost(temp_root);
                    else
                        new_cost = curr_cost;
                    break;
            }
            delta_cost = new_cost - curr_cost;
            MoveCount++;
            
            if(is_movable(delta_cost, T)) {
                if(Move_type == 3 && temp_root) {
                    delete_tree(root);
                    root = temp_root;
                }
                if(delta_cost > 0)
                    uphill++;

                double chip_R;
                if(root->blks[0].aspect.width > root->blks[0].aspect.height)
                    chip_R = root->blks[0].aspect.width / root->blks[0].aspect.height;
                else
                    chip_R = root->blks[0].aspect.height / root->blks[0].aspect.width;

                int chip_A = root->blks[0].aspect.width * root->blks[0].aspect.height;

                if(new_cost < best_cost && chip_R <= R_ub_temp) {

                    // cout << chip_A << " / " << chip_R << endl;
                    best_cost = new_cost;
                    best_SPE = SPE;

                    is_success = true;

                    top_down(root, 0);
                    chip_area = root->blks[0].aspect.width * root->blks[0].aspect.height;
                    chip_ratio = root->blks[0].aspect.width / root->blks[0].aspect.height;

                }
            }
            else {
                switch(Move_type) {
                    case 1:
                        M1(false, &iptr, &jptr, &i, &j);
                        break;
                    case 2:
                        M2(false, &iptr, &i, &j);
                        break;
                    case 3:
                        if(temp_root)
                            delete_tree(temp_root);                        
                        swap(i, j);
                        break;
                }
                reject++;
            }
        }
        while((uphill < N) && (MoveCount < 2 * N) /*&& (double(clock()-start)/CLOCKS_PER_SEC) < 7080*/);
        // cout << "T= " << T  << " time = " << double(clock()-start)/CLOCKS_PER_SEC << " s"<<endl;
        // cout << "MoveCount: " << MoveCount << " 2N: " << 2 * N << endl;
        // cout << "uphill: " << uphill << " N: " << N << endl << endl;

        T = w * T;
    }
    while(((double)reject/MoveCount) < 0.95 && (T > freezing));
    
    delete_tree(root);

    myFile.open(outputname, ios::out);
    write_file(myFile);
    myFile.close();

    return SPE;
}

void Floorplanning::top_down(Node *ptr, int index) {

    // calculate corrdinate

    if(!ptr)
        return;    
    if(ptr->data == "V" || ptr->data == "H") {

        if(ptr->data == "V") {
            ptr->left->x = ptr->x;
            ptr->left->y = ptr->y;
            ptr->right->x = ptr->x + ptr->left->blks[ptr->blks[index].Lindex].aspect.width; 
            ptr->right->y = ptr->y;
        }
        else {
            ptr->left->x = ptr->x;
            ptr->left->y = ptr->y;
            ptr->right->x = ptr->x;
            ptr->right->y = ptr->y + ptr->left->blks[ptr->blks[index].Lindex].aspect.height; 
        }

        if(ptr->left->data != "V" && ptr->left->data != "H") {
            module[stoi(ptr->left->data)].x = ptr->left->x;
            module[stoi(ptr->left->data)].y = ptr->left->y;

            module[stoi(ptr->left->data)].aspect.width = ptr->left->blks[ptr->blks[index].Lindex].aspect.width;
            module[stoi(ptr->left->data)].aspect.height = ptr->left->blks[ptr->blks[index].Lindex].aspect.height;
            
        }
        if(ptr->right->data != "V" && ptr->right->data != "H") {
            module[stoi(ptr->right->data)].x = ptr->right->x;
            module[stoi(ptr->right->data)].y = ptr->right->y;

            module[stoi(ptr->right->data)].aspect.width = ptr->right->blks[ptr->blks[index].Rindex].aspect.width;
            module[stoi(ptr->right->data)].aspect.height = ptr->right->blks[ptr->blks[index].Rindex].aspect.height;
        }
    }

    top_down(ptr->left, ptr->blks[index].Lindex);
    top_down(ptr->right, ptr->blks[index].Rindex);
}

void Floorplanning::write_file(fstream& myFile) {
    if(chip_ratio >= R_lb && chip_ratio <= R_ub) {
        myFile << "A = " << chip_area << endl;
        myFile << "R = " << chip_ratio << endl;
        for(int i = 0; i < module.size(); i++) {
            myFile << module[i].name << " " << module[i].x << " " << module[i].y;
            if(module[i].aspect.width == input_blk_aspect[i].height)
                myFile << " R";

            myFile << endl;
        }
    }
    else {
        myFile << "A = " << chip_area << endl;
        myFile << "R = " << 1 / chip_ratio << endl;
        for(int i = 0; i < module.size(); i++) {
            myFile << module[i].name << " " << module[i].y << " " << module[i].x;
            if(module[i].aspect.width != input_blk_aspect[i].height)
                myFile << " R";

            myFile << endl;
        }
    }

}

double Floorplanning::calculate_cost(Node *ptr) {
    double A = ptr->blks[0].aspect.width * ptr->blks[0].aspect.height;
    double A_norm = A / total_area-1;
    double R;

    if(ptr->blks[0].aspect.width > ptr->blks[0].aspect.height)
        R = ptr->blks[0].aspect.width / ptr->blks[0].aspect.height;
    else
        R = ptr->blks[0].aspect.height / ptr->blks[0].aspect.width;

    // double R_star = (R_ub_temp + 1) / 2;

    // cout << R_star << endl;

    double cost_R = (R <= R_ub_temp) ? 0 : R-1 ;

    // if(R <= R_ub_temp)
    //     cost_R = 0;
    // else 
    //     cost_R = abs(R-1);

    return cost_f * (alpha* A_norm + (1-alpha) * cost_R);

}

int Floorplanning::select_move() {
    return rand() % 3 + 1;
}

void Floorplanning::bottom_up(Node *ptr) {
    if(!ptr)
        return; 
    bottom_up(ptr->left);
    bottom_up(ptr->right);

    vector<Block> blks;
    Block temp;
    Size tmp;
    
    if((ptr->data == "V" || ptr->data == "H") && ptr->update == false) {
        ptr->update = true;
        sort(ptr->left->blks.begin(), ptr->left->blks.begin()+ptr->left->blks.size(), mycompare2);
        sort(ptr->right->blks.begin(), ptr->right->blks.begin()+ptr->right->blks.size(), mycompare2);

        if(ptr->data == "V") {
            int i = 0, j = 0;
            while(i < ptr->left->blks.size() && j < ptr->right->blks.size()){
                tmp.width = ptr->left->blks[i].aspect.width + ptr->right->blks[j].aspect.width;
                tmp.height = max(ptr->left->blks[i].aspect.height, ptr->right->blks[j].aspect.height);

                temp.aspect = tmp;
                temp.Lindex = i;
                temp.Rindex = j;
                blks.push_back(temp);

                if(is_equal(tmp.height, ptr->left->blks[i].aspect.height))
                    i++;
                else
                    j++;
            }  
        }
        else {
            int i = ptr->left->blks.size()-1, j = ptr->right->blks.size()-1;
            while(i >= 0 && j >= 0){
                tmp.height = ptr->left->blks[i].aspect.height + ptr->right->blks[j].aspect.height;
                tmp.width = max(ptr->left->blks[i].aspect.width, ptr->right->blks[j].aspect.width);

                temp.aspect = tmp;
                temp.Lindex = i;
                temp.Rindex = j;
                blks.push_back(temp);

                if(is_equal(tmp.width, ptr->left->blks[i].aspect.width))
                    i--;
                else
                    j--;
            }
        }

        ptr->blks.clear();

        sort(blks.begin(), blks.begin()+blks.size(), mycompare);

        for(int k = 0 ; k < blks.size(); k++)     // pts or blks.size(), this condition is changlable
             ptr->blks.push_back(blks[k]);
        
    }
}

void Floorplanning::read_file(fstream& myFile, string outputfile, clock_t start_T) {

    string in;

    double temp;
    Module temp_; 
    Size tmp;
    int cell_id = 0;

    start = start_T;

    temp_.x = temp_.y = 0;
    
    outputname = outputfile;

    myFile >> R_lb >> R_ub;

    R_ub_temp = R_ub > (1/R_lb) ? R_ub : (1/R_lb) ;
    // cout << "R_ub_temp: " << R_ub_temp << endl;

    while(myFile >> in) {
        temp_.name = in;

        myFile >> temp;
        tmp.width = temp;
        myFile >> temp;
        tmp.height = temp;

        input_blk_aspect.push_back(tmp);
        module.push_back(temp_);

        moduleCount += 1;
        total_area += tmp.width * tmp.height;
    }

    // for(int i = 0; i < input_blk_aspect.size(); i++) {
    //     cout << input_blk_aspect[i].width << " " << input_blk_aspect[i].height << endl;
    // }

}

vector<string> Floorplanning::init() {
    Node *tmp = NULL;
    vector<string> spe;

    for(int i = 0; i < 4 * moduleCount-2; i++) {
        tmp = creat_node("", false);
        tmp->right = avail;
        avail = tmp;
    }
    
    spe = init_SPE(spe);
    return spe;

}

void Floorplanning::print_tree_post(Node *ptr) {
    if(ptr) {
        print_tree_post(ptr->left);
        print_tree_post(ptr->right);
        cout << ptr->data << " ";
        // cout << endl;
        // for(int i = 0; i < ptr->blks.size(); i++) {
        //     cout << "width: " << ptr->blks[i].aspect.width << " height: " << ptr->blks[i].aspect.height << endl;
        // }
        // //  for(int i = 0; i < ptr->blks.size(); i++) {
        //     if(ptr->data == "H" || ptr->data == "V") {
        //         cout << "Lindex index: " << ptr->blks[0].Lindex << endl;
        //         cout << "Rindex index: " << ptr->blks[0].Rindex << endl;
        //     }
        // // }

        // if(ptr->data != "V" && ptr->data != "H") { // !
        //     cout << ptr->data << " ";
        //     cout << ptr->x << " " << ptr->y  << " width is " << module[stoi(ptr->data)].aspect.width << endl;
        // }
        // // !           FOR CHECKER           //
        // if(ptr->data == "V" || ptr->data == "H")
        //     myFile << ptr->data << " ";
        // else
        //     myFile << stoi(ptr->data)-1 << " ";
        // // !           FOR CHECKER           //
    }
}

void Floorplanning::M1(bool control, Node **iptr, Node **jptr, int *i, int *j) {

    int x;

    if(control){
        
        x = rand() % SPE.size();
        if(SPE[x] == "V" || SPE[x] == "H") 
            (*i) = find_int(x);
        else 
            (*i) = x;
        
        (*j) = find_int((*i));

        find_node_ij(&(*iptr), &(*jptr), (*i), (*j));
    }
    
    swap((*i), (*j));

    Node *parent_i = (*iptr)->parent;
    Node *parent_j = (*jptr)->parent;

    if(parent_i != parent_j) {
        if(parent_i->left == (*iptr))
            parent_i->left = (*jptr);
        else
            parent_i->right = (*jptr);
        if(parent_j->right == (*jptr))
            parent_j->right = (*iptr);
        else
            parent_j->left = (*iptr);
    }
    else {
        if(parent_i->left == (*iptr)){
            parent_i->left = (*jptr);
            parent_i->right = (*iptr);
        }
        else {
            parent_i->right = (*jptr);
            parent_i->left = (*iptr);
        }
    }

    (*iptr)->parent = parent_j;
    (*jptr)->parent = parent_i;
    
    for(Node *temp = (*iptr)->parent; temp; temp = temp->parent)
        temp->update = false;
    
    for(Node *temp = (*jptr)->parent; temp; temp = temp->parent) {
        temp->update = false;
        if(temp->parent && temp->parent->update == false)
            break;
    }

    bottom_up(root);
}

void Floorplanning::M2(bool control, Node **iptr, int *i, int *j) {

    int x;
    if(control){
        x = rand() % SPE.size();
    
        if(SPE[x] != "V" && SPE[x] != "H")
            x = find_VH(x);

        for(int k = x; k < SPE.size(); k++) {
            if(SPE[k] == "V" || SPE[k] == "H") 
                (*j) = k;
            else
                break;
        }
        for(int k = x; k >= 0; k--) {
            if(SPE[k] == "V" || SPE[k] == "H") 
                (*i) = k;
            else
                break;
        }
    }

    find_node_i(&(*iptr), (*i), (*j));

    for(int k = (*i); k <= (*j); k++) {
        if(SPE[k] == "V")
            SPE[k] = "H";
        else
            SPE[k] = "V";
    }

    for(Node *temp = (*iptr); temp; temp = temp->parent) 
        temp->update = false;

    bottom_up(root);
}

Node * Floorplanning::M3(int *i, int *j) {
    
    int x, counter = 0; 
    Node *temp = NULL;
    

    while(counter < 10) {
        x = rand() % (SPE.size() - 1);
        find_int_VH(x, &(*i), &(*j));
        swap((*i), (*j));
        if(!legal_check((*j))) 
            swap((*i), (*j));
        else{
            break;
        }
        counter++;
    }
        
    if(counter != 10) { 
        temp = SPE_to_tree(SPE);
        bottom_up(temp);
        return temp;
    }
    return NULL;
}

void Floorplanning::find_node_ij(Node **iptr, Node **jptr, int i, int j) {

    stack<Node *> s1, s2;
    Node *temp = NULL;
    int flag = 0;

    s1.push(root);
    while(!s1.empty()) {
        temp = s1.top();
        s1.pop();
        s2.push(temp);
        if(temp->left) 
            s1.push(temp->left);
        if(temp->right) 
            s1.push(temp->right);
    }
    while(!s2.empty()) {
        temp = s2.top();
        s2.pop();

        if(temp->data == SPE[i]) {
            (*iptr) = temp;
            if(!flag)
                flag = 1;
            else 
                break;
        }
        if(temp->data == SPE[j]) {
            (*jptr) = temp;
            if(!flag)
                flag = 1;
            else 
                break;
        }
    }
}

void Floorplanning::find_node_i(Node **iptr, int i, int j) {

    stack<Node *> s1, s2;
    Node *temp = NULL;
    int flag = 0;

    s1.push(root);
    while(!s1.empty()) {
        temp = s1.top();
        s1.pop();
        s2.push(temp);
        if(temp->left) 
            s1.push(temp->left);
        if(temp->right) 
            s1.push(temp->right);
    }
    while(!s2.empty()) {
        temp = s2.top();
        s2.pop();

        if(flag && i-j != 1) {
            if(temp->data == "V")
                temp->data = "H";
            else
                temp->data = "V";
            j--;
        }
        
        if(i - j == 1)
            break;

        if(temp->data == SPE[i-1]) {
            (*iptr) = temp->parent;
            flag = 1;
        }
    }

}

void Floorplanning::find_int_VH(int x, int *i, int *j) {
    
    if(SPE[x] == "V" || SPE[x] == "H") {
        for(int k = x; k < SPE.size()-1; k++) {
            if(SPE[k+1] != "V" && SPE[k+1] != "H") {
                (*i) = k;
                (*j) = k+1;
                return;
            }  
        }
        
        for(int k = x; k > 0; k--) {
            if(SPE[k-1] != "V" && SPE[k-1] != "H") {
                (*i) = k-1;
                (*j) = k;
                return;
            }  
            
        }
    }
    else {
        for(int k = x; k < SPE.size()-1; k++) {
            if(SPE[k+1] == "V" || SPE[k+1] == "H") {
                (*i) = k;
                (*j) = k+1;
                return;
            }  
        }
        for(int k = x; k > 0; k--) {
            if(SPE[k-1] == "V" || SPE[k-1] == "H") {
                (*i) = k-1;
                (*j) = k;
                return;
            }  
        }
    }
}

int Floorplanning::find_int(int x) {
    for(int k = x; k < SPE.size(); k++) {
        if(k == x) 
            continue;
        if(SPE[k] != "V" && SPE[k] != "H") {
            return k;
        }
    }
    for(int k = x; k >= 0; k--) {
        if(k == x) 
            continue;
        if(SPE[k] != "V" && SPE[k] != "H") {
            return k;
        }
    }
    return -1;
}

int Floorplanning::find_VH(int x) {
    for(int k = x; k < SPE.size(); k++) {
        if(k == x)
            continue;
        if(SPE[k] == "V" || SPE[k] == "H") {
            return k;
        }    
    }
    for(int k = x; k >= 0; k--) {
        if(k == x) 
            continue;
        if(SPE[k] == "V" || SPE[k] == "H") {
            return k;
        }
    }
    return -1;
}

bool Floorplanning::legal_check(int j) {
    int operandCount = 0, operatorCount = 0;

    if(j == SPE.size()-1) {
        if(SPE[j-2] == SPE[j-1]) 
            return false;
    }
    else if(j == 1) {
        if(SPE[j] == SPE[j+1]) 
            return false;
    }
    else {
        if(SPE[j-2] == SPE[j-1] || SPE[j] == SPE[j+1]) 
            return false;
    }

    if(SPE[j] != "V" && SPE[j] != "H") { // if right side is number, need to check balloting
        for(int i = 0; i <= j; i++) {
            if(SPE[i] == "V" || SPE[i] == "H")
                operatorCount++;
            else
                operandCount++;
            if(operatorCount >= operandCount) {
                return false;      
            }  
        }
    }
    return true;
}

vector<string> Floorplanning::init_SPE(vector<string> &spe) {
    spe.clear();
    spe.push_back("0");
    spe.push_back("1");
    spe.push_back("V");
    for(int i = 2; i < moduleCount; i++) {
        string x = to_string(i);
        spe.push_back(x);
        if(rand() % 2 == 0)
            spe.push_back("V");
        else
            spe.push_back("H");

    }

    // for(int i = 0; i < spe.size(); i++)
    //     cout << spe[i] << " ";

    // cout << endl;
    return spe;
}



Node *Floorplanning::creat_node(string data, bool is_avail) {
    Node *temp;
    if(!is_avail) 
        temp = new Node();
    else {
        temp = avail;
        avail = avail->right;
    }  
    temp->data = data;
    temp->left = NULL;
    temp->right = NULL;
    temp->update = false;
    temp->x = 0;
    temp->y = 0;
    return temp;
}

Node * Floorplanning::SPE_to_tree(vector<string> SPE) {

    stack<Node *> s;
    Node *temp = NULL;

    for(int i = 0; i < SPE.size(); i++) {
        if(SPE[i] == "V" || SPE[i] == "H") {
            if(SPE[i] == "V") 
                    temp = creat_node("V", true);
            else 
                    temp = creat_node("H", true);
            temp->right = s.top();
            temp->right->parent = temp;
            s.pop();

            temp->left = s.top();
            temp->left->parent = temp;
            s.pop();
            
            s.push(temp);
        }
        else { 
            temp = creat_node(SPE[i], true);
            choose_possible_size(stoi(SPE[i]), temp);
            s.push(temp);
        }
    }
    temp = s.top();
    temp->parent = NULL;
    s.pop();
    return temp;
}

void Floorplanning::choose_possible_size(int index, Node *ptr) { // no need to **!!!!!!!!!!!!!!!!
    
    Block temp;

    temp.aspect.width = input_blk_aspect[index].width;
    temp.aspect.height = input_blk_aspect[index].height;
    ptr->blks.push_back(temp);

    temp.aspect.width = input_blk_aspect[index].height;
    temp.aspect.height = input_blk_aspect[index].width;
    ptr->blks.push_back(temp);
}

void Floorplanning::swap(int i, int j) {
    string tmp = SPE[i];
    SPE[i] = SPE[j];
    SPE[j] = tmp;
}

bool mycompare(Block a, Block b){
   return a.aspect.width * a.aspect.height < b.aspect.width * b.aspect.height;
}

bool mycompare2(Block a, Block b){
   return a.aspect.width< b.aspect.width;
}
bool Floorplanning::is_equal(double a, double b) {
    if((a-b) < 0.00000001 && (a-b) > -0.00000001)
        return true;
    else
        return false;
}

void Floorplanning::delete_tree(Node *ptr) {
    if(ptr){
        delete_tree(ptr->left);
        delete_tree(ptr->right);
        if(ptr->parent) {
            if(ptr->parent->left == ptr)
                ptr->parent->left = NULL;
            else
                ptr->parent->right = NULL;
            ptr->parent = NULL;
        }
        ptr->blks.clear();
        ptr->right = avail;
        avail = ptr;
    }
}

bool Floorplanning::is_movable(double delta_cost, double T) {
    if(delta_cost <= 0)
        return true;
    
    double x = (double) rand() / (RAND_MAX + 1.0);

    if(x < exp(-(delta_cost/T)))
        return true;

    return false;
}

void Floorplanning::update_cooling_factor(double *w, double T) {
    // if(T <= 1000 && T > 66)
    //     *w = 0.7;
    // else if(T <= 66 && T > 33)
    //     *w = 0.85;
    // else if(T <= 33 && T > 16)
    //     *w = 0.95;
    // else if(T <= 16 && T > 6)
    //     *w = 0.96; 
    // else if(T <= 6 && T > 3)
    //     *w = 0.97;
    // else if(T <= 3 && T > 0.01)
    //     *w = 0.98;

    // if(T <= 30000 && T > 1000)
    //     *w = 0.6;
    // else if(T <= 1000 && T > 500)
    //     *w = 0.7;
    // else if(T <= 500 && T > 200)
    //     *w = 0.8;
    // else if(T <= 200 && T > 10)
    //     *w = 0.9; 
    // else if(T <= 10 && T > 0.1)
    //     *w = 0.97;
    // else if(T <= 0.1 && T > 0.001)
    //     *w = 0.98;

    if(T <= 30000 && T > 100)
        *w = 0.65;
    else if(T <= 100 && T > 10)
        *w = 0.9; 
    else if(T <= 10 && T > 0.1)
        *w = 0.95;
    else if(T <= 0.1 && T > 0.005)
        *w = 0.98;

}