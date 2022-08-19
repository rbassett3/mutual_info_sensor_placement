// #include <queue>
// #include <list>
// #include <cstring>
// #include <math.h>
// #include <iostream>
// #include "../scs/include/scs.h"
#include "cpp_version.h"
//load the shared library containing the Fortran functions (when using c++ in w/ jupyter)
//#pragma cling load("../fortran_funcs/fortran_funcs")
//#pragma cling load("../scs/out/libscsdir")

//define the struct that is put into the queue
//upper bound, free sensors, u \in {-1,0,1}, warm start info,
//and sensors placed
// struct queue_item {double ub;
//                     double* u; ScsSolution* sol; int sp;};

// //define a function so that we can sort the queue_inputs
// struct compare_queue_item{
//     bool operator()(const queue_item &item1, const queue_item &item2){
//     return item1.ub < item2.ub;
//     }
//     //< is reversed of > in previous line bc want to process max first
// };

// //load the external functions from the Fortran library
// //extern"C"{ void spm_greedy_(const double[], const int*, const int*, double[], double*);}
// extern"C"{ void spm_greedy_(double[], int*, int*, double[], double*);}
// extern"C"{ void find_best_in_subtree_(double[], double[], int*, int*, int*, double[], double*, int*);}
// extern"C"{ void choosebranchingvar_greedy_(double[], int*, double[], int*);}
// extern"C"{ void choosebranchingvar_max_var_(double[], int*, double[], int*);}

// extern"C"{ void fast_logdet_(double[], int*, double*);}
// extern"C"{ void best_from_many_outer_rot_roundings_(double[], int*, double[], int*, int*, double*, double*);}

// class BranchAndBound{
//     public:
//         //load the external functions from the Fortran library
//         //extern"C"{
//         //declarations of internal functions
//         //Arguments for default constructor is:
//         //Sigma, n, k, tol, depth, rounds
//         BranchAndBound(double*, int, int, double, int, int);
//         void flat2mat_nodiag(const int&, int&, int&);
//         void flat2mat(const int&, int&, int&);
//         void mat2flat(int&, const int&, const int&);
//         void bigmat2flat(int&, const int&, const int&);
//         void bigflat2mat(const int&, int&, int&);
//         void DeclareProblem();
//         void BeginBranch();
//         void TerminateBranch(double[], int, int);
//         void Branch(queue_item&, queue_item&, queue_item&);
//         void UpdateLbUb(queue_item&);
//         void Updateb(double*);
//         //declarations of branch and bound variables
//         int depth;
//         scs_int n;
//         scs_int k;
//         scs_int ind1;
//         scs_int ind2;
//         scs_int bigi;
//         scs_float *Sigma;
//         scs_float best_lb;
//         scs_float new_lb;
//         scs_float *best_u;
//         scs_float *new_u;
//         scs_float tol;
//         scs_float *A_vals;
//         scs_int A_rows;
//         scs_int A_cols;
//         scs_int *p;
//         scs_int *row_inds;
//         int rounds;
//         scs_float *X;
//         scs_float logdet_sigma;
//         std::priority_queue<queue_item, std::vector<queue_item>, compare_queue_item> Queue;
//         int sdps_solved = 0;
//         int leaf_nodes = 0;
//         bool warmstart = 1;

//         //declarations of SCS variables
//         ScsMatrix A;
//         ScsData *D = new ScsData;
//         ScsCone *K = new ScsCone;
//         ScsWork *W;
//         ScsSettings *stgs = new ScsSettings;
//         ScsInfo *info = new ScsInfo;
//         scs_float *b;
//         scs_float *c;

// };

void BranchAndBound::flat2mat_nodiag(const int &i, int &ind1, int &ind2){
    ind2 = floor((2*n-1 - sqrt(pow(2*n-1, 2) - 4*2*i))/2);
    ind1 = (i - (n*(n-1)/2 - (n-ind2)*(n-ind2-1)/2 - (ind2+1)));
}

void BranchAndBound::flat2mat(const int &i, int &ind1, int &ind2){
    ind2 = floor((2*n+1 - sqrt(pow(2*n+1,2) - 4*2*i))/2);
    ind1 = i - (n*(n+1)/2 - (n-ind2)*(n-ind2+1)/2 - ind2);
}

void BranchAndBound::mat2flat(int &i, const int &ind1, const int &ind2){
    if(ind1 >= ind2){
        i = n*(n+1)/2 - (n-ind2)*(n-ind2+1)/2 + ind1 - ind2;
        }
    else{
        i = n*(n+1)/2 - (n-ind1)*(n-ind1+1)/2 + ind2 - ind1;
        }
}

void BranchAndBound::bigmat2flat(int &i, const int &ind1, const int &ind2){
    if(ind1 >= ind2){
        i = 2*n*(2*n+1)/2 - (2*n-ind2)*(2*n-ind2+1)/2 + ind1 - ind2;
        }
    else{
        i = 2*n*(2*n+1)/2 - (2*n-ind1)*(2*n-ind1+1)/2 + ind2 - ind1;
        }
}

void BranchAndBound::bigflat2mat(const int &i, int &ind1, int &ind2){
    ind2 = floor((4*n+1 - sqrt(pow(4*n+1,2) - 4*2*i))/2);
    ind1 = i - (2*n*(2*n+1)/2 - (2*n-ind2)*(2*n-ind2+1)/2 - ind2);
}

void BranchAndBound::DeclareProblem(){
    //std::cout << "Entered DeclareProblem" << std::endl;;
    scs_int cur_elt;
    //declare the ScsMatrix A
    A_rows = 2*(n*(n+1)/2 - n) + 2*n + n*(n+1)/2 + n*(2*n+1) + 3*n;
    A_cols = n*(n+1)/2 - n + n*(n+1)/2 + n;
    int A_nnz = 8*((n*(n+1))/2 - n) + (n*(n+1))/2 + 2*n + n;
    //The entries of the matrix. Its length is
    //8*number of X + 2*number of entries of Z 
    //+ 2*number of diagonal entries of Z + n entries for the t variables
    A_vals = new scs_float[A_nnz];
    p = new scs_int[A_cols+1]; //pointers to where each column starts in M
    row_inds = new scs_int[A_nnz];
    //fill in the values of A
    cur_elt = 0;
    //std::cout << "Constructing X portion of A";
    //X variables
    for (scs_int i=0; i < (n*(n+1))/2 - n; ++i){
        p[i] = cur_elt;
        flat2mat_nodiag(i, ind1, ind2);
        //eltwise upper bound
        A_vals[cur_elt] = 1.0;
        row_inds[cur_elt] = i;
        //eltwise lower bound
        A_vals[cur_elt+1] = -1.0;
        row_inds[cur_elt+1] = i + (n*(n+1))/2 - n;
        //In what follows, I put in the mins and max b/c apparently
        //the row inds must be increasing
        //within each column
        //sum upper bound (row entry)
        A_vals[cur_elt+2] = 1.0;
        row_inds[cur_elt+2] = 2*((n*(n+1))/2 - n) + std::min(ind1, ind2);
        //sum upper bound (col entry)
        A_vals[cur_elt+3] = 1.0;
        row_inds[cur_elt+3] = 2*((n*(n+1))/2 - n) + std::max(ind1, ind2);
        //sum lower bound (row entry)
        A_vals[cur_elt+4] = -1.0;
        row_inds[cur_elt+4] = 2*((n*(n+1))/2 - n) + n + std::min(ind1, ind2);
        //sum lower bound (col entry)
        A_vals[cur_elt+5] = -1.0;
        row_inds[cur_elt+5] = 2*((n*(n+1))/2 - n) + n + std::max(ind1, ind2);
        //X is psd.
        mat2flat(bigi, ind1, ind2);
        A_vals[cur_elt+6] = -sqrt(2);
        row_inds[cur_elt+6] = bigi + 2*((n*(n+1))/2 - n) + 2*n;
        //concatenated X and Z is psd
        bigmat2flat(bigi, ind1+n, ind2+n);
        //fancy array arithmetic b/c Sigma is read as stacked rows
        //but b/c its symmetic this will also work if its column stacked.
        A_vals[cur_elt+7] = -sqrt(2)*Sigma[ind1*n + ind2]/2.0;
        row_inds[cur_elt+7] = bigi + 2*((n*(n+1))/2 - n) + 2*n + (n*(n+1))/2;
        cur_elt += 8;
    }

    //std::cout << "Constructing Z portion of A";
    //The Z variables
    for (scs_int i=0; i<(n*(n+1))/2; ++i){
        p[i+ (n*(n+1))/2 - n] = cur_elt;
        //concatenated X and Z is psd
        flat2mat(i, ind1, ind2);
        if (ind1 == ind2){ //on the diagonal
            //diagonal piece of matrix in the concatenated X and Z
            bigmat2flat(bigi, ind1, ind2);
            A_vals[cur_elt] = -1.0;
            row_inds[cur_elt] = bigi + 2*((n*(n+1))/2 - n) +
                                2*n + (n*(n+1))/2;
            cur_elt += 1;
        }
        bigmat2flat(bigi, ind1+n, ind2);
        A_vals[cur_elt] = -sqrt(2);
        row_inds[cur_elt] = bigi + 2*((n*(n+1))/2 - n) + 2*n + (n*(n+1))/2;
        cur_elt += 1;

        if (ind1 == ind2){ //on the diagonal
            //exponential cone constraint
            A_vals[cur_elt] = -1.0;
            row_inds[cur_elt] = 2*(n*(n+1)/2 - n) + 2*n +
                                  n*(n+1)/2 + n*(2*n+1) +
                                  3*ind1+2;
            cur_elt += 1;
        }
    }
    //std::cout << "Constructing t portion of A" ;
    //the t variables
    for (scs_int i=0; i<n; i++){
          p[i + (n*(n+1))/2 - n + (n*(n+1))/2] = cur_elt;
          A_vals[cur_elt] = -1.0;
          row_inds[cur_elt] = 2*(n*(n+1)/2 - n) + 2*n +
                              n*(n+1)/2 + n*(2*n+1) + 3*i;
          cur_elt += 1;
        }
    p[A_cols] = cur_elt;
    A = (ScsMatrix){.x=A_vals, .i=row_inds, .p=p, .m=A_rows, .n=A_cols};

    //construct b
    b = new scs_float[2*(n*(n+1)/2 -n) + 2*n + n*(n+1)/2 + n*(2*n+1) + 3*n];
    //input the 2*(n*(n+1)//2-n) eltwise constraints.
    //first upper bound
    scs_int cur_row = 0;
    for (scs_int i=0; i<(n*(n+1))/2 - n; ++i){
        flat2mat_nodiag(i, ind1, ind2);
        //upper bound
        b[cur_row+i] = 1.0;
        //lower bound
        b[cur_row + (n*(n+1))/2 - n + i] = 1.0;
        }
    cur_row += 2*((n*(n+1))/2 - n);
    //input the 2n sum constraints
    //first ub
    int max_sum_bound = std::max(2*k-n, n-2*k);
    for (scs_int i=0; i<n; ++i){
        //upper bound
        b[cur_row + i] = max_sum_bound - 1.0;
        //lower bound
        b[cur_row + n + i] = max_sum_bound + 1.0;
        }
    cur_row += 2*n;
    //input the psd constraints
    //X is psd
    for (scs_int i=0; i<(n*(n+1))/2; ++i){
        flat2mat(i, ind1, ind2);
        if (ind1 == ind2){
            b[cur_row + i] = 1.0;
            }
        else{
            b[cur_row+i] = 0.0;
        }
    }
    cur_row += (n*(n+1))/2;
    //concatenated X and Z is psd
    for (scs_int i=0; i<n*(2*n+1); ++i){
        bigflat2mat(i, ind1, ind2);
        if (ind1 >= n and ind2 >= n){
            if (ind1 == ind2){
                b[cur_row + i] = Sigma[n*(ind1-n) + ind2-n];
            }
            else{
                b[cur_row + i] = sqrt(2)*Sigma[n*(ind1-n) + ind2-n]/2;
            }
        }
        else{
        b[cur_row+i] = 0.0;
        }
    }
    cur_row += n*(2*n+1);

    //add exponential cone constraint
    //https://docs.mosek.com/modeling-cookbook/expo.html
    for (scs_int i=0; i<3*n; ++i){
        if (i%3==1){
            b[cur_row + i] = 1.0;
        }
        else{
            b[cur_row + i] = 0.0;
        }
    }
        
    //construct c
    c = new scs_float[A_cols];
    for (scs_int i=0; i<A_cols; i++){
        if (i >= A_cols-n){
            c[i] = -1.0;
        }
        else{
            c[i] = 0.0;
        }
    }

    //assemble the problem data
    //D = &(ScsData){.m=A_rows, .n=A_cols, .A=&A, .P=SCS_NULL, .b=b, .c=c};
    D->m=A_rows;
    D->n=A_cols;
    D->A=&A;
    D->P=SCS_NULL;
    D->b=b;
    D->c=c;
    //print A
    /*
    std::cout << "The A matrix's nonzero values" << std::endl;
    for (int i=0; i<8*((n*(n+1))/2 - n) + (n*(n+1))/2 + 2*n + n; i++){
        std::cout << A_vals[i] << std::endl;
    }
    

    std::cout << "The A matrix's row indices" << std::endl;
    for (int i=0; i<8*((n*(n+1))/2 - n) + (n*(n+1))/2 + 2*n + n; i++){
        std::cout << row_inds[i] << std::endl;
    }

    std::cout << "The A matrix's indices where new col begins" << std::endl;
    for (int i=0; i<A_cols+1; i++){
        std::cout << p[i] << std::endl;
    }

    std::cout << "The b vector" << std::endl;
    for (int i=0; i<2*(n*(n+1)/2 - n) + 2*n + n*(n+1)/2 + n*(2*n+1) + 3*n; i++){
        std::cout << b[i] << std::endl;
    }

    std::cout << "The c vector" << std::endl;
    for (int i=0; i<n*(n+1)/2 - n + n*(n+1)/2 + n; i++){
        std::cout << c[i] << std::endl;
    }
    */
        
    //D = new ScsData {m=A_rows, n=A_cols, A=A, P=SCS_NULL, b=b, c=c};
    int s[2] = {n, 2*n};
    //K = &(ScsCone){.z=0, .l=2*(n*(n+1)/2 -n + n), .bu=SCS_NULL, .bl=SCS_NULL,
    //                    .bsize=0, .q=SCS_NULL, .qsize=0,
    //                    .s=s, .ssize=2, .ep=n, .ed=0, .p=SCS_NULL, .psize=0};
    K->z=0;
    K->l=2*(n*(n+1)/2 -n + n);
    K->bu=SCS_NULL;
    K->bl=SCS_NULL;
    K->bsize=0;
    K->q=SCS_NULL;
    K->qsize=0;
    K->s=s;
    K->ssize=2;
    K->ep=n;
    K->ed=0;
    K->p=SCS_NULL;
    K->psize=0;
    W = scs_init(D, K, stgs);
    scs_update(W, b, c);

}          

scs_int BranchAndBound::BeginBranch(){
    int nz, num_sensors_to_place;
    long leaf_nodes_cut;
    int status;
    //initialize with greedy (from Fortran)
    spm_greedy_(Sigma, &n, &k, best_u, &best_lb);
    std::cout << "Initial (Greedy) Lower Bound: " << best_lb << std::endl;
    //for (int i = 0; i < n; i++){
    //    std::cout << best_u[i];
    //}
    //std::cout << std::endl;
    //declare initial queue item
    //This should fill in the content of best_u and best_lb
    //std::cout << "Queue declared and first item pushed" << std::endl;
    DeclareProblem();

    queue_item cur_item, pos_item, neg_item;
    cur_item.ub = 1e6;
    cur_item.sp = 0;
    cur_item.u = new scs_float[n];
    memset(cur_item.u, 0.0, n*sizeof(scs_float));
    cur_item.sol = new ScsSolution;
    cur_item.sol->x = new scs_float[A_cols];
    cur_item.sol->y = new scs_float[A_rows];
    cur_item.sol->s = new scs_float[A_rows];
    Queue.push(cur_item); //doing this copies cur_item into the queue, so we can reuse it.
    //std::cout << "Problem declared" << std::endl;
    while (!Queue.empty()){
        cur_item = Queue.top();
        //std::cout << cur_item.ub << std::endl;
        Queue.pop(); //does this deallocate?
        //std::cout << "Popped item from queue" << std::endl;
        //check if ub/lb make this branch valid
        if (cur_item.ub >= best_lb){
            //count number of free_sensors remaining
            nz = 0;
            for (int i=0; i<n; i++){
                if ((cur_item.u)[i] == 0.0)
                    nz += 1;
            }
            //std::cout << "Computed size" << std::endl;
            //std::cout << cur_item.sp << std::endl;
            //either we've solved SDPs to required depth
            //or we are one step away from determining all placements
            if ((n - nz >= depth) || (cur_item.sp >= k-1) || (n - nz - cur_item.sp >= n-k-1)){
                std::cout << "Terminating branch \n" << std::endl;
                //Either sdp depth exceeded or all sensors placed.
                TerminateBranch(cur_item.u, nz, cur_item.sp);
                //deallocate 
                delete cur_item.u;
                delete cur_item.sol->x;
                delete cur_item.sol->y;
                delete cur_item.sol->s;
                delete cur_item.sol;
                //delete cur_item;
            }
            else{
                //std::cout << "Branching" << std::endl;
                pos_item = cur_item; //I hope this copies
                neg_item = cur_item; //I hope this copies
                status = Branch(cur_item, pos_item, neg_item, nz);
                if(status != 0){
                    std::cout<<"ERROR! SOLVER DID NOT TERMINATE CORRECTLY" << std::endl;
                    return -1;
                }
            }
        }
        else{
            std::cout << "Branch Cut!" << std::endl;
            num_sensors_to_place = k - cur_item.sp;
            choose_(&nz, &num_sensors_to_place, &leaf_nodes_cut);
            std::cout << "Cut " << leaf_nodes_cut << " leaf nodes" << '\n' << std::endl;
            delete cur_item.u;
            delete cur_item.sol->x;
            delete cur_item.sol->y;
            delete cur_item.sol->s;
            delete cur_item.sol;
            //delete cur_item;
        }
    }
    scs_finish(W);
    std::cout << "Solution is: " << std::endl;
    for(int i=0; i<n; i++){
    std::cout << best_u[i] << " , "; 
    }
    std::cout << std::endl;
    clear_memory();
    return 0;
}

scs_int BranchAndBound::UpdateLbUb(queue_item &cur_item){
    //update b with d_lb, d_ub, B_lb, and B_ub
    //std::cout << "Entered UpdateLbUb function" << std::endl;
    std::cout << "SDP solve " << sdps_solved << std::endl;
    Updateb(cur_item.u);

        //std::cout << "b has been updated" << std::endl;
    
    
    //update problem workspace
    //last argument is c but it doesn't get updated so it can be set to SCS_NULL
    scs_update(W, b, c);
    //std::cout << "Work has been updated" << std::endl;
    //solve problem
    //last argument is whether to warm start. 
    //Seems like it's b/c x,s,y are not allocated
    //std::cout << "Solving Problem" << std::endl;
    //cur_item.sol.x = new scs_float[A_cols];
    //cur_item.sol.y = new scs_float[A_rows];
    //cur_item.sol.s = new scs_float[A_rows];

    /*
    if(cur_item.sol->x = SCS_NULL){

            std::cout<< "Warning! Warm starting when sol not allocated" << std::endl;
        }
    if (initial_solve){
        cur_item.sol->x = new scs_float[A_cols];
        cur_item.sol->y = new scs_float[A_rows];
        cur_item.sol->s = new scs_float[A_rows];
        cur_item.sol->x[A_cols-1] = 0.0;
        scs_solve(W, cur_item.sol, info, 1);
        scs_solve(W, cur_item.sol, info, 1);
        initial_solve = 0;
    }
    else{
        cur_item.sol->x = new scs_float[A_cols];
        cur_item.sol->y = new scs_float[A_rows];
        cur_item.sol->s = new scs_float[A_rows];
        cur_item.sol->x[A_cols-1] = 0.0;
        if(cur_item.sol->x = SCS_NULL){

            std::cout<< "Warning! Warm starting when sol not allocated" << std::endl;
        }
        scs_solve(W, cur_item.sol, info, 1);
        scs_solve(W, cur_item.sol, info, 1);
    }
    */
    
    //if (!warmstart){ //allocate 
        //cur_item.sol->x = new scs_float[A_cols];
        //cur_item.sol->y = new scs_float[A_rows];
        //cur_item.sol->s = new scs_float[A_rows];
        //cur_item.sol->x[A_cols-1] = 0.0;
    //    }

    scs_solve(W, cur_item.sol, info, warmstart);
    warmstart = 1;
    //std::cout << "completed first solve" << std::endl;
    //scs_solve(W, cur_item.sol, info, 1);
    //std::cout << "completed second solve" << std::endl;
    //

   if (info->status_val != 1){
       std::cout << "SDP not solved to optimality" << std::endl;
       std::cout << "Status_val is " << info->status_val << std::endl;
       std::cout << "Resolving without warmstart" << std::endl;
       scs_solve(W, cur_item.sol, info, 0);
       std::cout << "New status_val is " << info->status_val << std::endl;
       if (info->status_val != 1){ //solver failed
           return -1;
       }
    }
    sdps_solved += 1;
    //std::cout << "Problem has been solved" << std::endl;
    //std::cout << "Here's the lower triangular part of X" << std::endl;
    //for (int i=0; i < (n*(n-1))/2; i++){
    //    std::cout << cur_item.sol.x[i] << std::endl;
    //}
    cur_item.ub = -(info->dobj) - logdet_sigma; 
    //std::cout<< "primal objective " << -(info->pobj) - logdet_sigma << std::endl; 
    //thought about adding tol here due to corner cases in the feasibility tolerance.
    //decided against it. This can make the dual objective not actually upper bound the 
    // true maximum. It's rare and when it does so it's only by a very small margin.
    std::cout << "(lower bound, upper bound): (" << best_lb << ", " << cur_item.ub << ")" << std::endl;
    //use the negative dual here b/c it bounds -(obj value) from above
    
    //TODO:
    if(rounds > 0 && best_lb < cur_item.ub){
       /* 
        std::cout << "X on entry " << std::endl;
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                std::cout<< X[n*i+j] << " ";
            }
            std::cout<< std::endl;
        }
        */
        
        for(int i=0; i<(n*(n+1))/2 - n; i++){
            flat2mat_nodiag(i, ind1, ind2);
            //not assigning diagonals to be ones
            //not using symmetry
            X[ind1*n + ind2] = cur_item.sol->x[i];
            ////////////////X[ind2*n + ind1] = cur_item.sol->x[i];
        }

        
        /* 
        std::cout << "X after assignment" << std::endl;
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                std::cout<< X[n*i+j] << " ";
            }
            std::cout<< std::endl;
        }
        */
        

        best_from_many_outer_rot_roundings_(Sigma, &n, X, &rounds, &k, new_u, &new_lb);
        std::cout << "Best lb from roundings: " << new_lb << std::endl;
        std::cout << "Best u from roundings: [";
        for(int i=0; i<n; i++){
            std::cout << new_u[i] << ", ";
        }
        std::cout << "]" << std::endl;
        if(new_lb > best_lb){
            best_lb = new_lb;
            memcpy(best_u, new_u, n*sizeof(scs_float));
        }
    }
    return 0;
}

void BranchAndBound::TerminateBranch(double* u, int nz, int sp){
    int num_to_place = k - sp;
    //std::cout << "Entered TerminateBranch function" << std::endl;
    //std::cout << "u is: " << std::endl;
    //for (int i=0; i<n; i++){ 
    //    std::cout << u[i] << " ";
    //}
    //std::cout << std::endl;
    //std::cout << nz << std::endl;
    //std::cout << num_to_place << std::endl;
    find_best_in_subtree_(Sigma, u, &n, &num_to_place, &nz, new_u, &new_lb, &leaf_nodes);
    //std::cout << "found best in subtree" << std::endl;
    if(new_lb - 1e-8 > best_lb){
        std::cout << "Improved lower bound: " << new_lb << '\n' << std::endl;
        //std::cout << "u is: [" << std::endl;
        //for (int i=0; i<n; i++){
        //    std::cout<< new_u[i] << ", ";
        //    }
        //std::cout << "]" << std::endl;
        memcpy(best_u, new_u, n*sizeof(scs_float));
        best_lb = new_lb;
    }
}

scs_int BranchAndBound::Branch(queue_item &item, queue_item &pos_item, queue_item &neg_item, int nz){
    int i, num_sensors_to_place;
    long leaf_nodes_cut;
    int status;
    //std::cout << "In Branch function" << std::endl;

    // choose branching variable based on given input
    switch(branch_var){

        // add in bad greedy method for analysis
        case 1:
            choosebranchingvar_greedy_bad_(Sigma, &n, item.u, &i);
	    break;

        case 2:
            choosebranchingvar_max_var_(Sigma, &n, item.u, &i);
            break;

        case 3:
            choosebranchingvar_random_(&n, item.u, &i);
            break;

        default:
            choosebranchingvar_greedy_good_(Sigma, &n, item.u, &i);

    }
    //std::cout << "branch_var is " << branch_var << std::endl;
    std::cout << "Chose branching variable as " << i << std::endl;
    //set ith sensor to 1
    //copy cur_item
    //std::cout << "setting branching variable to 1" << std::endl;
    pos_item.u = new scs_float[n];
    memcpy(pos_item.u, item.u, n*sizeof(scs_float));
    pos_item.sol = new ScsSolution;
    pos_item.sol->x = new scs_float[A_cols];
    pos_item.sol->y = new scs_float[A_rows];
    pos_item.sol->s = new scs_float[A_rows];
    (pos_item.u)[i] = 1.0;
    pos_item.sp += 1; //we've placed a sensor
    status = UpdateLbUb(pos_item);
    if(status != 0){return -1;}
    //push these to queue if the bounds are appropriate
    if (pos_item.ub >= best_lb){
        Queue.push(pos_item);
        //std::cout << "pushed positive branch to queue" << std::endl;
    }
    else{
        std::cout << "Branch Cut" << std::endl;
        num_sensors_to_place = k - pos_item.sp;
        choose_(&nz, &num_sensors_to_place, &leaf_nodes_cut);
        std::cout << "Cut " << leaf_nodes_cut << " leaf nodes\n" << std::endl;
        delete[] pos_item.u;
        delete[] pos_item.sol->x;
        delete[] pos_item.sol->y;
        delete[] pos_item.sol->s;
        delete pos_item.sol;
    }


       
    //set ith sensor to -1
    //std::cout << "setting branching variable to -1" << std::endl;
    //except that copying arrays only copies their pointers
    //so we reallocate and copy the original u
    //neg_item.sol.x = new scs_float[A_cols];
    //neg_item.sol.y = new scs_float[A_rows];
    //neg_item.sol.s = new scs_float[A_rows];

    neg_item.u = new scs_float[n];
    memcpy(neg_item.u, item.u, n*sizeof(scs_float));
    neg_item.sol = new ScsSolution;
    neg_item.sol->x = new scs_float[A_cols];
    neg_item.sol->y = new scs_float[A_rows];
    neg_item.sol->s = new scs_float[A_rows];
    (neg_item.u)[i] = -1.0;
    status = UpdateLbUb(neg_item);
    if(status != 0){return -1;}
    if (neg_item.ub >= best_lb){
        Queue.push(neg_item);
        //std::cout << "pushed negative branch to queue" << std::endl;
    }
    else{
        std::cout << "Branch Cut" << std::endl;
        num_sensors_to_place = k - neg_item.sp;
        choose_(&nz, &num_sensors_to_place, &leaf_nodes_cut);
        std::cout << "Cut " << leaf_nodes_cut << " leaf nodes\n" << std::endl;
        delete[] neg_item.u;
        delete[] neg_item.sol->x;
        delete[] neg_item.sol->y;
        delete[] neg_item.sol->s;
        delete neg_item.sol;
    }
    
    delete[] item.u;
    delete[] item.sol->x;
    delete[] item.sol->y;
    delete[] item.sol->s;
    delete item.sol;
    //delete item;
    return 0;
}

void BranchAndBound::Updateb(double* u){

    scs_int cur_row = 0;
    //the eltwise constraints
    for (scs_int i=0; i<(n*(n+1))/2 - n; ++i){
        flat2mat_nodiag(i, ind1, ind2);
        if (u[ind1]*u[ind2] != 0.0){
            //upper bound
            b[cur_row+i] = u[ind1]*u[ind2];
            //lower bound
            b[cur_row + (n*(n+1))/2 - n + i] = -u[ind1]*u[ind2];
        }
        else{
            //upper bound
            b[cur_row + i] = 1.0;
            //lower bound
            b[cur_row + (n*(n+1))/2 - n + i] = 1.0;
        }
    }
    double max_sum_bound = std::max(2*k-n, n-2*k);
    cur_row += 2*((n*(n+1))/2 - n);
    //input the 2n sum constraints
    //first ub
    for (scs_int i=0; i<n; ++i){
        switch(int(u[i])){
            case 1: 
                //upper bound
                b[cur_row + i] = 2.0*k - n - 1.0 + 1e-6;
                //lower bound
                b[cur_row + n + i] = -(2.0*k - n - 1.0);
            break;
            case -1:
                //upper bound
                b[cur_row + i] = -2.0*k + n - 1.0 + 1e-6;
                //lower bound
                b[cur_row + n + i] = -(-2.0*k + n - 1.0);
            break;
            case 0:
                //upper bound
                b[cur_row + i] = max_sum_bound - 1.0;
                //lower bound
                b[cur_row + n + i] = max_sum_bound + 1.0;
            break;
            default:
                std::cout << "You should never be here" << std::endl;
        }
    }
    //std::cout << "The u vector" << std::endl;
    //std::cout << "[";
    //for (int i=0; i<n; i++){
    //    std::cout << u[i] << ", ";
    //}
    //std::cout << "]" << std::endl;

    //std::cout << "The first elts of the b vector" << std::endl;
    //std::cout << "[";
    //for (int i=0; i<2*(n*(n+1)/2 - n) + 2*n; i++){
    //[    std::cout << b[i] << ", ";
    //}
    //std::cout << "]" << std::endl;


}

//Sigma, n, k, tol, depth, rounds
BranchAndBound::BranchAndBound(double *sigma, int N, int num_to_place, double Tol, int Depth, int Rounds, int BranchVar){
    n = N;
    //Sigma = new scs_float[int(pow(n,2))];
    Sigma = sigma;
    k = num_to_place;
    tol = Tol;
    depth = Depth;
    rounds = Rounds;
    branch_var = BranchVar;
    best_u = new scs_float[n];
    new_u = new scs_float[n];
    X = new scs_float[n*n];
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if(i==j){
                X[n*i+i] = 1.0;
            }
            else{
                X[n*i+j] = 0.0;
            }
        }
    }
    best_lb = -1e6; //just fix some initial value that's very small
    //set default settings
    scs_set_default_settings(stgs);
    stgs->verbose = 0;
    //this mirrors what cvxpy does
    //=https://github.com/cvxpy/cvxpy/blob/master/cvxpy/reductions/solvers/conic_solvers/scs_conif.py
    stgs->eps_abs = tol;
    stgs->eps_rel = 0;
    //stgs->alpha = 1.8; //scs defaults to 1.5, cvxpy to 1.8
    //stgs->scale = 5.0; //scs defaults to 0.1, cvxpy to 5.0
    stgs->max_iters = 100000;
    stgs->warm_start = 0;
    //stgs->eps_infeas = 1e-12;
    //stgs->acceleration_lookback = 0;
    //copy Sigma and compute its logdet
    scs_float* Sigma_cp;
    Sigma_cp = new scs_float[int(pow(n,2))];
    for (int i=0; i< pow(n,2); i++){
        Sigma_cp[i]=Sigma[i];
    }
    fast_logdet_(Sigma_cp, &n, &logdet_sigma);
    delete[] Sigma_cp;
}

//todo: write a function that deallocates everything
void BranchAndBound::clear_memory(){
    delete[] best_u;
    delete[] new_u;
    delete[] X;
    delete[] p;
    delete[] A_vals;
    delete[] row_inds;
}

int cython_wrapper(double* Sigma, int n, int k, double tol, int depth, int rounds, std::string filename, int branch_var){
    BranchAndBound BnB = BranchAndBound(Sigma, n, k, tol, depth, rounds, branch_var);
    scs_int status;
    //BranchAndBound *BnB = new BranchAndBound(Sigma, n, k, tol, depth, rounds);
    //BnB->BeginBranch();
    status = BnB.BeginBranch();

    std::ofstream myfile(filename.c_str(), std::ios::app);
    if(myfile.is_open())
    {
        if(status != 0){myfile << "ERROR! Solver did not terminate correctly" << std::endl;}

        switch(branch_var){
            case 1:
                myfile << "onion_" << "greedy_bad_" << n << std::endl;
                break;

            case 2:
                myfile << "onion_" << "max_var_" << n << std::endl;
                break;

            case 3:
                myfile << "onion_" << "random_" << n << std::endl;
                break;

            default:
                myfile << "onion_" << "greedy_good_" << n << std::endl;
                break;
        }

        myfile << BnB.leaf_nodes << std::endl;
        myfile << BnB.sdps_solved << std::endl;
        myfile.close();
    }
    return status;
}
