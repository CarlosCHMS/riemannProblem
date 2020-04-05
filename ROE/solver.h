#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

typedef double FTYPE;

struct solverStruct{

    int Ns;
    int N;
    int fluxType;

    FTYPE** U;
    FTYPE** newU;
    //FTYPE** F;
    FTYPE R;
    FTYPE g;
    FILE* output;
    //FTYPE* sqrtRho;
    FTYPE** flux;
    FTYPE cc;
    FTYPE saveStep;
    FTYPE eLim;
    FTYPE CFL;

    FTYPE AUSMa;
    FTYPE AUSMb;

};

void solverAlloc(struct solverStruct* solver);

void solverInit(struct solverStruct* solver, struct inputStruct* input);

void solverInitU(struct solverStruct* solver, int n, FTYPE p0, FTYPE T0, FTYPE p1, FTYPE T1);

void solverPrint(struct solverStruct* solver, int jj);

void solverCalcF(struct solverStruct* solver);

void solverPropagate(struct solverStruct* solver);

void solverSimulate(struct solverStruct* solver);

// ROE functions

void calcUPH(struct solverStruct* solver, int ii,  FTYPE* u, FTYPE* p, FTYPE* h);

void solverCalcFluxROE(struct solverStruct* solver);

FTYPE solverEntropyFixROE(struct solverStruct* solver, FTYPE e);

void solverUpdateCFLROE(struct solverStruct* solver, FTYPE e0, FTYPE e2);

// functions created for AUSMp application

void solverCalcFluxAUSM(struct solverStruct* solver);

void calcUPHA(struct solverStruct* solver, int ii,  FTYPE* u, FTYPE* p, FTYPE* h, FTYPE* a);

void solverCalcSplitP(struct solverStruct* solver, FTYPE* M, FTYPE* Mp, FTYPE* Pp);

void solverCalcSplitN(struct solverStruct* solver, FTYPE* M, FTYPE* Mn, FTYPE* Pn);

void solverCheckSplit(struct solverStruct* solver, int flag);

void solverUpdateCFL(struct solverStruct* solver, FTYPE e0);

#endif // SOLVER_H_INCLUDED
