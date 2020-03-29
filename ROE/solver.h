#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

typedef double FTYPE;

struct solverStruct{

    int Ns;
    int N;
    FTYPE** U;
    FTYPE** newU;
    FTYPE** F;
    FTYPE R;
    FTYPE g;
    FILE* output;
    FTYPE* sqrtRho;
    FTYPE** flux;
    FTYPE cc;
    FTYPE saveStep;
    FTYPE eLim;
    FTYPE CFL;

};

void solverAlloc(struct solverStruct* solver);

void solverInit(struct solverStruct* solver, struct inputStruct* input);

void solverInitU(struct solverStruct* solver, int n, FTYPE p0, FTYPE T0, FTYPE p1, FTYPE T1);

void solverPrint(struct solverStruct* solver, int jj);

void calcUPH(struct solverStruct* solver, int ii,  FTYPE* u, FTYPE* p, FTYPE* h);

void solverCalcF(struct solverStruct* solver);

void solverCalcFlux(struct solverStruct* solver);

void solverPropagate(struct solverStruct* solver);

void solverSimulate(struct solverStruct* solver);

FTYPE solverEntropyFix(struct solverStruct* solver, FTYPE e);

void solverUpdateCFL(struct solverStruct* solver, FTYPE e0, FTYPE e2);

#endif // SOLVER_H_INCLUDED
