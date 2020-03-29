#include<stdio.h>
#include<stdlib.h>
#include "input.h"
#include "solver.h"

int main(){

    struct inputStruct* input = malloc(sizeof(struct inputStruct));
    struct solverStruct* solver = malloc(sizeof(struct solverStruct));

    inputInit(input);

    inputPrintParameters(input);

    solverInit(solver, input);

    solverSimulate(solver);

    return 0;

}
