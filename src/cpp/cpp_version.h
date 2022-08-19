#ifndef MY_CLASS_H // include guard
#define MY_CLASS_H

#include <queue>
#include <list>
#include <cstring>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "../scs/include/scs.h"

//#include "../scs/include/scs.h"
//load the shared library containing the Fortran functions (when using c++ in w/ jupyter)
//#pragma cling load("../fortran_funcs/fortran_funcs")
//#pragma cling load("../scs/out/libscsdir")

//define the struct that is put into the queue
//upper bound, free sensors, u \in {-1,0,1}, warm start info,
//and sensors placed
struct queue_item {scs_float ub;
                    scs_float* u; ScsSolution* sol; scs_int sp;};

//define a function so that we can sort the queue_inputs
struct compare_queue_item{
    bool operator()(const queue_item &item1, const queue_item &item2){
    return item1.ub < item2.ub;
    }
    //< is reversed of > in previous line bc want to process max first
};

//load the external functions from the Fortran library
//extern"C"{ void spm_greedy_(const double[], const int*, const int*, double[], double*);}
extern"C"{ void spm_greedy_(double[], int*, int*, double[], double*);}
extern"C"{ void find_best_in_subtree_(double[], double[], int*, int*, int*, double[], double*, int*);}
extern"C"{ void choosebranchingvar_greedy_good_(double[], int*, double[], int*);}
extern"C"{ void choosebranchingvar_greedy_bad_(double[], int*, double[], int*);}
extern"C"{ void choosebranchingvar_max_var_(double[], int*, double[], int*);}
extern"C"{ void choosebranchingvar_random_(int*, double[], int*);}

extern"C"{ void fast_logdet_(double[], int*, double*);}
extern"C"{ void best_from_many_outer_rot_roundings_(double[], int*, double[], int*, int*, double*, double*);}
extern"C"{ void choose_(int*, int*, long*);}

class BranchAndBound{
    public:
        //load the external functions from the Fortran library
        //extern"C"{
        //declarations of internal functions
        //Arguments for default constructor is:
        //Sigma, n, k, tol, depth, rounds
        BranchAndBound(scs_float*, scs_int, scs_int, scs_float, scs_int, scs_int, scs_int);
        void flat2mat_nodiag(const scs_int&, scs_int&, scs_int&);
        void flat2mat(const scs_int&, scs_int&, scs_int&);
        void mat2flat(scs_int&, const scs_int&, const scs_int&);
        void bigmat2flat(scs_int&, const scs_int&, const scs_int&);
        void bigflat2mat(const scs_int&, scs_int&, scs_int&);
        void DeclareProblem();
        scs_int BeginBranch();
        void TerminateBranch(scs_float[], scs_int, scs_int);
        scs_int Branch(queue_item&, queue_item&, queue_item&, scs_int);
        scs_int UpdateLbUb(queue_item&);
        void Updateb(scs_float*);
        void clear_memory();
        //declarations of branch and bound variables
        scs_int depth;
        scs_int n;
        scs_int k;
        scs_int ind1;
        scs_int ind2;
        scs_int bigi;
        scs_float *Sigma;
        scs_float best_lb;
        scs_float new_lb;
        scs_float *best_u;
        scs_float *new_u;
        scs_float tol;
        scs_float *A_vals;
        scs_int A_rows;
        scs_int A_cols;
        scs_int *p;
        scs_int *row_inds;
        scs_int rounds;
        scs_float *X;
        scs_float logdet_sigma;
        scs_int branch_var;
        std::priority_queue<queue_item, std::vector<queue_item>, compare_queue_item> Queue;
        scs_int sdps_solved = 0;
        int leaf_nodes = 0;
        bool warmstart = 0;

        //declarations of SCS variables
        // ScsMatrix *A = SCS_NULL;
        ScsMatrix A;
        ScsData *D = new ScsData;
        ScsCone *K = new ScsCone;
        ScsWork *W = SCS_NULL;
        ScsSettings *stgs = new ScsSettings;
        ScsInfo *info = new ScsInfo;
        scs_float *b = SCS_NULL;
        scs_float *c = SCS_NULL;
        scs_int s[2];
};

int cython_wrapper(scs_float*, scs_int, scs_int, scs_float, scs_int, scs_int, std::string, scs_int);


#endif /* cpp_version.h */
