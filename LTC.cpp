#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <fstream>
#include <time.h>
#include "BOBHASH32.h"
using namespace std;


//-------------------------------------------------------------------------------------
// global default settings
// string dataset = "ddos_2.dat";
string dataset = "stack";
int num_periods = 1000, co_p = 3, co_f = 10, top_k = 100, window_size, experiment_idx = 0;
int th_p = 30, th_f = 100;
bool THRESHOLD = true; bool has_auxiliary_array = false;
const int MAX_ENTRIES = 13000005;
const int MAX_ITEMS = 1000005;
const int MAX_BUCKETS = 10005; // 10K buckets will consume 640K memory
const int LONG_TAIL_BOUND = 100000;

// map <string, int> deducted_count, enter_count;
// int ratio_hist[10000];
//-------------------------------------------------------------------------------------


class LongTailClock {
    private:
    BOBHash32 *bobhash; // hash function
    int X, Y;           // X, Y is the position of the pointer in the CLOCK algorithm
    float step_size;    // the time that the pointer takes to move to the next cell
    int Last_T;         // when the last insetion happens
    const int ODD=1;
    
    public:
    struct CELL {int persistency; bool odd, even; string ID; int count;};
    CELL ** LTC;
    set<string> sig_items; // used for finding significant items with threshold
    int M;
    bool freq_replacement;
    bool two_level;
    int M1, M2;
    int cells_per_bucket;
    
    LongTailClock(int MEM, bool _freq_replacement, int _cells_per_bucket, int hash_seed, float first_layer) {
        freq_replacement = _freq_replacement;
        cells_per_bucket = _cells_per_bucket;
        two_level = !(first_layer == 1);
        bobhash = new BOBHash32(hash_seed);
        if (THRESHOLD) {
            // if using auxiliary array, there is additional memory overhead
            if (has_auxiliary_array)
                MEM = MEM * 0.8;
        }
        M = MEM * 1024 / cells_per_bucket / 8;  //the number of buckets
        if (two_level) {
            M1 = M * first_layer;
            M2 = M - M1;
        }
        // M = MEM * 1024 / CELLS_PER_BUCKET / 16;
        printf("number of buckets M: %d\n", M);
        step_size = (window_size) / (M * cells_per_bucket+0.0);
        LTC = new CELL*[cells_per_bucket];
        for (int i = 0; i < cells_per_bucket; i ++) {
            LTC[i] = new CELL[M];
        }
        for (int i=0; i<cells_per_bucket; i++)
            for (int j=0; j<M; j++) 
                LTC[i][j].persistency = LTC[i][j].odd = LTC[i][j].even = LTC[i][j].count=0;
        Last_T=0;
        X=Y=0;
    }

    ~LongTailClock() {
        for (int i = 0; i < cells_per_bucket; i ++) {
            delete []LTC[i];
        }
        delete []LTC;
    }

    int get_scanned_time() {
        return int(step_size * (X*M + Y));
    }

    void insert(string x, int T)
    {
        if (THRESHOLD && has_auxiliary_array) {
            if (sig_items.find(x) != sig_items.end())
                return;
        }
        int P_idx = (Last_T / window_size) % 2; // judge whether this period is an even-numbered period or an odd-numbered period
        if (T/window_size == Last_T / window_size) { // if the item x does not cause the increment of period
            // the pointer p moves clockwise
            while (get_scanned_time() <= T % window_size && Y!=M) {
                if (P_idx==ODD && LTC[X][Y].even) 
                    LTC[X][Y].persistency++, LTC[X][Y].even=0; 
                if (P_idx!=ODD && LTC[X][Y].odd)
                    LTC[X][Y].persistency++, LTC[X][Y].odd=0;
                X++; 
                if (X==cells_per_bucket) 
                    X=0, Y++;
            }
            if (Y==M)
                Y=0;
        } else { // if the item x casues the increment of period
            // similarly, the pointer p moves clockwise
            for (int i = Last_T / window_size; i < T / window_size; i++) {
                while (Y!=M) {
                    if (P_idx==ODD && LTC[X][Y].even)
                        LTC[X][Y].persistency++, LTC[X][Y].even=0;
                    if (P_idx!=ODD && LTC[X][Y].odd)
                        LTC[X][Y].persistency++, LTC[X][Y].odd=0;
                    X++; 
                    if (X==cells_per_bucket) 
                        X=0, Y++;
                }
                Y=0;
                P_idx ^= 1;
            }
            while (get_scanned_time() <= T % window_size && Y != M) {
                if (P_idx==ODD && LTC[X][Y].even) 
                    LTC[X][Y].persistency++, LTC[X][Y].even=0;
                if (P_idx!=ODD && LTC[X][Y].odd)
                    LTC[X][Y].persistency++, LTC[X][Y].odd=0;
                X++;
                if (X==cells_per_bucket)
                    X=0, Y++;
            }
            if (Y==M)
                Y=0;
        }

        int min_index = -1, min_value = -2147483648, second_min_value=-2147483648;
        bool has_item = false;
        unsigned int p, q;
        p = bobhash->run(x.c_str(), x.size());
        // calculate the index of the mapped bucket
        if (!two_level) {
            p = p % M;
        } else {
            q = p % M2 + M1;
            p = p % M1;
            for (int k=0; k<cells_per_bucket; k++) {
                if (LTC[k][q].ID == x) {
                    if (P_idx == ODD) 
                        LTC[k][q].odd = 1; 
                    else
                        LTC[k][q].even = 1; 
                    LTC[k][q].count++; 
                    has_item = true;
                }
            }
        }
        if (!has_item) {
            for (int k=0; k<cells_per_bucket; k++) {
                if (LTC[k][p].ID == x) {
                    if (P_idx == ODD) 
                        LTC[k][p].odd = 1; 
                    else
                        LTC[k][p].even = 1; 
                    LTC[k][p].count++; 
                    has_item = true;
                    if (THRESHOLD) {
                        if (LTC[k][p].count >= th_f && LTC[k][p].persistency >= th_p-1) {
                            if (has_auxiliary_array) {
                                LTC[k][p].count = 0;
                                LTC[k][p].persistency = 0;
                            }
                            sig_items.insert(LTC[k][p].ID);
                        }
                    }
                    break;
                } // if x is found in that bucket
                else // calculate the smallest cell
                {
                    int current_value = (LTC[k][p].persistency + 
                        LTC[k][p].odd + LTC[k][p].even)*co_p + LTC[k][p].count*co_f;
                    if (min_index == -1 || current_value < min_value) {
                        min_index = k;
                        second_min_value = min_value;
                        min_value = current_value;
                    }
                    // deducted_count[x] ++;
                }
            }
        }
        if (!has_item) { // if x is not found in that bucket
            if (two_level) {
                q = bobhash->run(LTC[min_index][p].ID.c_str(), LTC[min_index][p].ID.size());
                q = q % M2 + M1;
                bool has_item_q = false;
                int min_index_q = -1, min_value_q = -2147483648;
                for (int k=0; k<cells_per_bucket; k++) {
                    if (LTC[k][q].ID == x) {
                        has_item_q = true;
                        break;
                    }
                    int current_value = (LTC[k][q].persistency + 
                        LTC[k][q].odd + LTC[k][q].even)*co_p + LTC[k][q].count*co_f;
                    if (min_index_q == -1 || current_value < min_value_q) {
                        min_index_q = k;
                        min_value_q = current_value;
                    }
                }
                if (!has_item_q && min_value > min_value_q) {
                    LTC[min_index_q][q].persistency = LTC[min_index][p].persistency;
                    LTC[min_index_q][q].count = LTC[min_index][p].count;
                    LTC[min_index_q][q].odd = LTC[min_index][p].odd;
                    LTC[min_index_q][q].even = LTC[min_index][p].even;
                    LTC[min_index_q][q].ID = LTC[min_index][p].ID;
                    // ----------------------------------------------
                    // LTC[min_index][p].persistency = 0;
                    // LTC[min_index][p].count = 0;
                    // ----------------------------------------------
                }
            }
            // Singificance Decrementing operation
            if (!freq_replacement)
                LTC[min_index][p].persistency--; 
            LTC[min_index][p].count--;

            if (LTC[min_index][p].persistency < 0)
                LTC[min_index][p].persistency=0;
            if (LTC[min_index][p].count < 0)
                LTC[min_index][p].count=0;
            // judge whether that item could be replaced with x
            // persistency is always <= frequency so we never need to look at persistency
            if (LTC[min_index][p].count*co_f <= 0) {
                LTC[min_index][p].ID = x;
                int MINC=LONG_TAIL_BOUND, MINcount=LONG_TAIL_BOUND;
                for (int k=0; k<cells_per_bucket; k++)  // Long-tail Replacement technique
                    if (k != min_index) {
                        MINC=min(MINC, LTC[k][p].persistency);
                        MINcount=min(MINcount, LTC[k][p].count);
                    }
                LTC[min_index][p].persistency = max(0, MINC-1);
                LTC[min_index][p].count = max(1, MINcount-1);
                LTC[min_index][p].odd = LTC[min_index][p].even = 0;
                if (P_idx == ODD)
                    LTC[min_index][p].odd=1;
                else
                    LTC[min_index][p].even=1;
                // enter_count[x] = LTC[min_index][p].count;
            }
        }
        // cout << "ok" << endl;
        Last_T = T;
    }

    void clean_up_persistency()
    {
        for (int k=0; k<cells_per_bucket; k++)
            for (int p=0; p<M; p++) {
                // the estimated persistency should be incremented if one of the flags is not equal to 0
                LTC[k][p].persistency += LTC[k][p].even+LTC[k][p].odd;
                LTC[k][p].even = LTC[k][p].odd = 0;
            }
    }
};


//---------------------------------------------------------------------------------
// input data
int n_entries = 0;
char* ID[MAX_ENTRIES];
int timestamp[MAX_ENTRIES];

// data structures used by calculate_real_significance() and passed to run_ltc()
map <string, int> items_of_period, real_persistency, item_count, top_significance, all_significance;
map <string, bool> if_ddos;
set<string> real_sig_items;
struct NODE{string x; int y;} real_significance[MAX_ITEMS], LTC_significance[MAX_ITEMS];
int bin_stat[32]; float bin_total;
int ddos_total = 0;
float ddos_recall;
//---------------------------------------------------------------------------------


bool cmp(NODE i, NODE j) {return i.y > j.y;}


void read_file() {
    n_entries = 0;
    if (dataset=="stack") {
        // there are m flows in total, every flow contains an ID and its arriving time.
        freopen("stack-new.txt", "r", stdin);
        char id[12];
        for (int i=1; i<=MAX_ENTRIES; i++)
        {
            if (scanf("%s", id) == EOF)
                break;
            ID[i] = new char[12];
            strcpy(ID[i], id);
            scanf("%d", &timestamp[i]);
            n_entries ++;
        }
    } else if (dataset.substr(0, 4) == "ddos") {
        freopen(dataset.c_str(), "r", stdin);
        int temp;
        char id[16];
        for (int i=1; i<=MAX_ENTRIES; i++)
        {
            ID[i] = new char[16];
            if (scanf("%s", id) == EOF)
                break;
            scanf("%d", &timestamp[i]);
            scanf("%d", &temp);
            strcpy(ID[i], id);
            if_ddos[string(ID[i])] = bool(temp);
            n_entries ++;
        }
    }
    else {
        cout << "dataset not found" << endl;
        return;
    }
    // initialize the time.
    for (int i=n_entries; i>=1; i--)
        timestamp[i] -= timestamp[1];
    // window_size denotes the time span of one period
    window_size = timestamp[n_entries] / num_periods;
    cout << "read lines: " << n_entries << endl;
}


void calculate_real_significance() {
    // Last denotes the arriving time of the previous flow
    int Last = 0;
    items_of_period.clear();
    item_count.clear();
    real_persistency.clear();
    top_significance.clear();
    all_significance.clear();

    // prepare for the real top-k significant items
    for (int i=1; i<=n_entries; i++)
    {
        string s=ID[i]; int T=timestamp[i];
        if (T / window_size != Last / window_size)  // meet a new period
        {
            for (map<string,int>::iterator sit=items_of_period.begin(); sit!=items_of_period.end(); sit++)
                real_persistency[sit->first]++;  //real_persistency[i] denotes the persistency of i
            items_of_period.clear();
        }
        Last = T;
        // items_of_period[i] denotes whether i has appeared in this period
        items_of_period[s]++; item_count[s]++;
    }
    for (map<string,int> :: iterator sit = items_of_period.begin(); sit!=items_of_period.end(); sit++)
        real_persistency[sit->first]++;  // execute the last period
    // the process of preparing ends

    int CNT=0;
    top_significance.clear(); all_significance.clear();
    for (map<string,int>::iterator sit = item_count.begin(); sit != item_count.end(); sit ++)
    {
        if (dataset.substr(0, 4) == "ddos" && if_ddos[sit->first]) {
            ddos_total ++;
        }
        if (THRESHOLD) {
            if (sit->second >= th_f && real_persistency[sit->first] >= th_p) {
                real_sig_items.insert(sit->first);
            }
        }
        real_significance[++CNT].x = sit->first; // means the ID
        real_significance[CNT].y = real_persistency[sit->first]*co_p + sit->second*co_f; // means the significance
    }
    cout << real_sig_items.size() << endl;

    sort(real_significance+1, real_significance+CNT+1, cmp); // sort all values

    for (int i=1; i<=top_k; i++) {
        top_significance[real_significance[i].x] = real_significance[i].y;
    }
    for (int i=1; i<=CNT; i++) {
        all_significance[real_significance[i].x] = real_significance[i].y;
    }
    cout << "calculating ground truth is done" << endl;
}


void run_ltc(string mode, int MEM, bool freq_replacement, int cells_per_bucket,
             float first_layer, FILE* f_pre, FILE* f_are) {   
    LongTailClock ltc = LongTailClock(MEM, freq_replacement,
                                      cells_per_bucket, experiment_idx++, first_layer);
    for (int i=1; i<=n_entries; i++) {
        string s=ID[i];
        int T=timestamp[i];
        ltc.insert(s, T); // an item S arrives at time T
    }
    ltc.clean_up_persistency(); // execute the last period

    // get all the significances recorded in the lossy table
    int CNT = 0;
    for (int i=0; i<ltc.M; i++) {
        int topK_counter = 0;
        for (int j=0; j<cells_per_bucket; j++) {
            LTC_significance[++CNT].x=ltc.LTC[j][i].ID;
            LTC_significance[CNT].y=ltc.LTC[j][i].persistency*co_p + ltc.LTC[j][i].count*co_f;
            if (top_significance[LTC_significance[CNT].x]) {
                topK_counter ++;
            }
        }
        bin_stat[topK_counter] ++;
        bin_total += 1;
    }

    sort(LTC_significance+1, LTC_significance+CNT+1, cmp); // get the estimated top-k significant items
    
    double PRE = 0, ARE = 0.0; // calculate the precision and ARE of LTC
    int ddos_cnt = 0;
    for (int i=1; i<=top_k; i++)
    {
        if (if_ddos[LTC_significance[i].x])
            ddos_cnt += 1;
        if (top_significance[LTC_significance[i].x])
            PRE += 1.0; // means it is indeed a significant item
        ARE += fabs(all_significance[LTC_significance[i].x] - LTC_significance[i].y) 
                    / (all_significance[LTC_significance[i].x] + 0.0);
        // printf("id: %s, esti: %d, real: %d\n",
        //     LTC_significance[i].x.c_str(), LTC_significance[i].y, all_significance[LTC_significance[i].x]);
    }

    if (ddos_total != 0)
        ddos_recall = ddos_cnt / float(ddos_total);
    if (THRESHOLD) {
        float recall = 0;
        for (set<string>::iterator sit = real_sig_items.begin(); sit != real_sig_items.end(); sit ++) {
            if (ltc.sig_items.find(*sit) != ltc.sig_items.end()) {
                recall ++;
            }
        }
        recall = recall/real_sig_items.size();
        float precision = 0;
        for (set<string>::iterator sit = ltc.sig_items.begin(); sit != ltc.sig_items.end(); sit ++) {
            if (real_sig_items.find(*sit) != real_sig_items.end()) {
                precision ++;
            }
        }
        precision = precision/ltc.sig_items.size();
        fprintf(f_pre, "%s,%d,%.5f\n", mode.c_str(), MEM, 2*precision*recall/(precision+recall));
    }
    if (f_are != NULL && f_pre != NULL) {
        ARE/=top_k; PRE/=top_k;
        printf("mode:%s, precision: %.5f, ARE: %.5f\n", mode.c_str(), PRE, ARE);
        fprintf(f_pre, "%s,%d,%.5f\n", mode.c_str(), MEM, PRE);
        fprintf(f_are, "%s,%d,%.5f\n", mode.c_str(), MEM, ARE);
    }
    /* 
    for (map<string,int>::iterator sit = enter_count.begin(); sit != enter_count.end(); sit ++)
    {
        int ratio = float(sit->second) / deducted_count[sit->first] * 100;
        if (sit->second > 100)
            ratio_hist[ratio] += 1;
    }
    for (int i = 0; i <= 150; i ++) {
        printf("%d,%d\n", i, ratio_hist[i]);
        ratio_hist[i] = 0;
    }
    enter_count.clear();
    deducted_count.clear();
    */
}


// measuring: precision, are, statistics of buckets
void test_v1_v2_mem(string f1, string f2) {
    FILE * f_pre = fopen(f1.c_str(), "w");
    fprintf(f_pre, "Algo,Memory size (KB),Precision\n");
    FILE * f_are = fopen(f2.c_str(), "w");
    fprintf(f_are, "Algo,Memory size (KB),ARE\n");
    FILE * f_stat = fopen("freq_replacement/ltc_stat.txt", "w");
    fprintf(f_stat, "Memory,Number of top $k$ items in the bucket,Percentage of buckets\n");
    
    int NUM_REPEAT = 5;
    int cells_per_bucket = 8;
    for (int MEM = 10; MEM <= 50; MEM += 10) {
        for (int i = 0; i < cells_per_bucket; i ++) {
            bin_stat[i] = 0;
            bin_total = 0.0;
        }
        for (int i = 0; i < NUM_REPEAT; i ++) {   
            run_ltc("LTC\\_F", MEM, true, cells_per_bucket, 1, f_pre, f_are);
            run_ltc("LTC", MEM, false, cells_per_bucket, 1, f_pre, f_are);
        }
        for (int i = 0; i < cells_per_bucket; i ++) {
            fprintf(f_stat, "memory=%dK,%d,%.5f\n", MEM, i, bin_stat[i]/bin_total);
            fflush(f_stat);
        }
    }
    fclose(f_stat);
    fclose(f_pre);
    fclose(f_are);
}

// measuring: precision, are, speed
void test_v2_v3_mem(string f1, string f2) {
    FILE * f_pre = fopen(f1.c_str(), "w");
    fprintf(f_pre, "Algo,Memory size (KB),Precision\n");
    FILE * f_are = fopen(f2.c_str(), "w");
    fprintf(f_are, "Algo,Memory size (KB),ARE\n");
    
    int NUM_REPEAT = 10;
    int cells_per_bucket = 8;
    for (int MEM = 10; MEM <= 50; MEM += 10) {
        for (int i = 0; i < NUM_REPEAT; i ++) {   
            run_ltc("50\\\% secondary part", MEM, true, cells_per_bucket, 0.5, f_pre, f_are);
            run_ltc("40\\\% secondary part", MEM, true, cells_per_bucket, 0.6, f_pre, f_are);
            run_ltc("30\\\% secondary part", MEM, true, cells_per_bucket, 0.7, f_pre, f_are);
            run_ltc("15\\\% secondary part", MEM, true, cells_per_bucket, 0.85, f_pre, f_are);
            run_ltc("no secondary part", MEM, true, cells_per_bucket, 1, f_pre, f_are);
        }
    }
    fclose(f_pre);
    fclose(f_are);
}


void test_threshold() {
    FILE * f_pre = fopen("freq_replacement/threshold.txt", "w");
    fprintf(f_pre, "Algo,Memory size (KB),$F_1$ score\n");
    
    int NUM_REPEAT = 1;
    int cells_per_bucket = 8;
    for (int MEM = 10; MEM <= 50; MEM += 10) {
        for (int i = 0; i < NUM_REPEAT; i ++) {   
            has_auxiliary_array = true;
            run_ltc("LTC\\_Y", MEM, true, cells_per_bucket, 1, f_pre, NULL);
            has_auxiliary_array = false;
            run_ltc("LTC\\_N", MEM, true, cells_per_bucket, 1, f_pre, NULL);
        }
    }
    fclose(f_pre);
}


void test_v1_v3_time() {
    FILE * f_time = fopen("two_layers/time.txt", "w");
    fprintf(f_time, "Algo,Time\n");
    
    int NUM_REPEAT = 5;
    int cells_per_bucket = 8;
    int MEM = 30;

    clock_t t1 = clock();
    cout << t1 << endl;
    for (int i = 0; i < NUM_REPEAT; i ++) {   
        run_ltc("50\\\% secondary part", MEM, true, cells_per_bucket, 0.5, NULL, NULL);
    }
    fprintf(f_time, "ltc_two_levels,%f\n", (clock() - t1) * 1.0 / CLOCKS_PER_SEC * 1000);
    t1 = clock();
    cout << t1 << endl;
    for (int i = 0; i < NUM_REPEAT; i ++) {   
        run_ltc("no secondary part", MEM, true, cells_per_bucket, 1, NULL, NULL);
    }
    cout << clock() << endl;
    fprintf(f_time, "ltc_two_levels,%f\n", (clock() - t1) * 1.0 / CLOCKS_PER_SEC * 1000);
    fclose(f_time);
}

// measuring: precision, are
void test_v1_cells_per_bucket() {
    FILE * f_pre = fopen("cells_per_bucket/ltc_pre.txt", "w");
    fprintf(f_pre, "Algo,Memory size (KB),Precision\n");
    FILE * f_are = fopen("cells_per_bucket/ltc_are.txt", "w");
    fprintf(f_are, "Algo,Memory size (KB),ARE\n");
    
    int NUM_REPEAT = 5;
    for (int MEM = 10; MEM <= 50; MEM += 10) {
        for (int i = 0; i < NUM_REPEAT; i ++) {   
            run_ltc("$d$=4", MEM, true, 4, 1, f_pre, f_are);
            run_ltc("$d$=8", MEM, true, 8, 1, f_pre, f_are);
            run_ltc("$d$=12", MEM, true, 12, 1, f_pre, f_are);
            run_ltc("$d$=16", MEM, true, 16, 1, f_pre, f_are);
        }
    }
    fclose(f_pre);
    fclose(f_are);
}

// measuring: recall
void test_v1_ddos_detection() {
    FILE * f_ddos = fopen("ddos_detection/recall.txt", "w");
    fprintf(f_ddos, "Algo,$k$,recall rate\n");
    
    int NUM_REPEAT = 1;
    int MEM = 30;
    num_periods = 100000;
    int top_k_list[5] = {100, 500, 1000, 1500, 2000};
    for (int j = 0; j < 5; j ++) {
        top_k = top_k_list[j];
        for (int i = 0; i < NUM_REPEAT; i ++) {
            co_p = 1; co_f = 1;
            run_ltc("$\\alpha$=1 $\\beta$=1", MEM, false, 8, 1, NULL, NULL);
            fprintf(f_ddos, "$\\alpha$=1 $\\beta$=1,%d,%f\n", top_k, ddos_recall);
            co_p = 1; co_f = 0;
            run_ltc("$\\alpha$=0, $\\beta$=1", MEM, false, 8, 1, NULL, NULL);
            fprintf(f_ddos, "$\\alpha$=0 $\\beta$=1,%d,%f\n", top_k, ddos_recall);
            co_p = 0; co_f = 1;
            run_ltc("$\\alpha$=1, $\\beta$=0", MEM, false, 8, 1, NULL, NULL);
            fprintf(f_ddos, "$\\alpha$=1 $\\beta$=0,%d,%f\n", top_k, ddos_recall);
            co_p = -1; co_f = 1;
            run_ltc("$\\alpha$=1, $\\beta$=-1", MEM, false, 8, 1, NULL, NULL);
            fprintf(f_ddos, "$\\alpha$=1 $\\beta$=-1,%d,%f\n", top_k, ddos_recall);
            fflush(f_ddos);
        }
    }
    num_periods = 1000;
    top_k = 100;
    fclose(f_ddos);
}


int main() {
    cout << "running" << endl;

    // co_f corresponds to the alpha in the paper, co_p corresponds to the beta in the paper.
    // cin >> co_f >> co_p;
    // top_k means the number of significant items we want to retrieve.
    // cin >> top_k;
    // MEM means the memomy size (KB)
    // cin >> MEM;
    // num_periods means the number of periods
    // cin >> num_periods;

    read_file();
    calculate_real_significance();
    if (dataset == "stack") {
        test_threshold();
        // test_v1_v2_mem("freq_replacement/ltc_pre.txt", "freq_replacement/ltc_are.txt");
        // test_v2_v3_mem("two_layers/ltc_pre.txt", "two_layers/ltc_are.txt");
        // test_v1_v3_time();
        // test_v1_cells_per_bucket();
        // co_p = 1; co_f = 0;
        // calculate_real_significance();
        // test_v1_v2_mem("freq_replacement/ltc_pre_p.txt", "freq_replacement/ltc_are_p.txt");
        // co_p = 0; co_f = 1;
        // calculate_real_significance();
        // test_v1_v2_mem("freq_replacement/ltc_pre_f.txt", "freq_replacement/ltc_are_f.txt");
    } else if (dataset.substr(0, 4) == "ddos") {
        test_v1_ddos_detection();
    }

    printf("# individual items: %d, # all items: %d\n", item_count.size(), n_entries);
    return 0;
}
