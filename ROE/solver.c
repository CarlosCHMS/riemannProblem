#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "input.h"
#include "solver.h"


void solverAlloc(struct solverStruct* solver){

    int ii;

    solver->U = (FTYPE**)malloc(solver->N*sizeof(FTYPE*));
    solver->newU = (FTYPE**)malloc(solver->N*sizeof(FTYPE*));
    solver->F = (FTYPE**)malloc(solver->N*sizeof(FTYPE*));
    solver->flux = (FTYPE**)malloc((solver->N-1)*sizeof(FTYPE*));

    solver->sqrtRho = (FTYPE*)malloc(solver->N*sizeof(FTYPE));

    for(ii=0; ii<solver->N; ii++){
        solver->U[ii] = (FTYPE*)malloc(3*sizeof(FTYPE));
        solver->newU[ii] = (FTYPE*)malloc(3*sizeof(FTYPE));
        solver->F[ii] = (FTYPE*)malloc(3*sizeof(FTYPE));

        if(ii<(solver->N-1)){
            solver->flux[ii] = (FTYPE*)malloc(3*sizeof(FTYPE));
        };

    };

}

void solverInit(struct solverStruct* solver, struct inputStruct* input){

    solver->Ns = input->Ns;
    solver->N = input->Np;
    solver->R = input->R;
    solver->g = input->gamma;
    solver->cc = input->cc;
    solver->saveStep = input->saveStep;
    solver->eLim = 1e-6;

    solver->output = fopen("./output.csv", "w");

    solverAlloc(solver);

    solverInitU(solver, input->n, input->p0, input->T0, input->p1, input->T1);

    solverPrint(solver, 0);

    solverSimulate(solver);

    fclose(solver->output);

}

void solverInitU(struct solverStruct* solver, int n, FTYPE p0, FTYPE T0, FTYPE p1, FTYPE T1){

    // rho: volumetric density
    // e: volumetric internal energy
    // U[0]: density
    // U[1]: mass flux
    // U[2]: volumetric total internal energy

    int ii;
    FTYPE u, rho0, e0, rho1, e1;

    u = 0.0;

    rho0 = p0/(solver->R*T0);
    e0 = rho0*(solver->R*T0/(solver->g-1) + u*u/2);

    rho1 = p1/(solver->R*T1);
    e1 = rho1*(solver->R*T1/(solver->g-1) + u*u/2);

    //printf("%f, %f,\n", rho0, e0);

    for(ii=0; ii<solver->N; ii++){
        if(ii < n){

            solver->U[ii][0] = rho0;
            solver->U[ii][1] = 0.0;
            solver->U[ii][2] = e0;

        }else{

            solver->U[ii][0] = rho1;
            solver->U[ii][1] = 0.0;
            solver->U[ii][2] = e1;

        };

    };

}

void solverPrint(struct solverStruct* solver, int jj){

    int ii;

    fprintf(solver->output, "%i, %i, %f\n", solver->N, jj, solver->CFL);

    for(ii=0;ii<solver->N;ii++){

        fprintf(solver->output, " %f,", solver->U[ii][0]);

    };
    fprintf(solver->output, "\n");

    for(ii=0;ii<solver->N;ii++){

        fprintf(solver->output, " %f,", solver->U[ii][1]);

    };
    fprintf(solver->output, "\n");

    for(ii=0;ii<solver->N;ii++){

        fprintf(solver->output, " %f,", solver->U[ii][2]);

    };
    fprintf(solver->output, "\n");

}

void calcUPH(struct solverStruct* solver, int ii,  FTYPE* u, FTYPE* p, FTYPE* h){

    //u: velocity
    //p: pressure
    //h: massic entalpy

    *u = solver->U[ii][1]/solver->U[ii][0];
    *p = (solver->U[ii][2] - solver->U[ii][1]*(*u)/2)*(solver->g - 1);
    *h = (solver->U[ii][2] + *p)/solver->U[ii][0];

}

void solverCalcF(struct solverStruct* solver){

    int ii;
    FTYPE u, p, h;

    for(ii=0; ii<solver->N; ii++){

        calcUPH(solver, ii, &u, &p, &h);
        solver->F[ii][0] = solver->U[ii][1];
        solver->F[ii][1] = solver->U[ii][1]*u + p;
        solver->F[ii][2] = solver->U[ii][1]*h;

        solver->sqrtRho[ii] = sqrt(solver->U[ii][0]);

    };

}

void solverCalcFlux(struct solverStruct* solver){

    FTYPE u0, p0, h0, u1, p1, h1, rhom, rd, um, hm, cm;
    FTYPE alp0, alp1, alp2, e0, e1, e2;
    FTYPE* vec0;
    FTYPE* vec1;
    FTYPE* vec2;

    vec0 = (FTYPE*)malloc(3*sizeof(FTYPE));
    vec1 = (FTYPE*)malloc(3*sizeof(FTYPE));
    vec2 = (FTYPE*)malloc(3*sizeof(FTYPE));

    int ii, jj;

    for(ii=0; ii<(solver->N-1); ii++){

        calcUPH(solver, ii, &u0, &p0, &h0);
        calcUPH(solver, ii+1, &u1, &p1, &h1);

        rhom = solver->sqrtRho[ii]*solver->sqrtRho[ii+1];

        rd = (solver->sqrtRho[ii] + solver->sqrtRho[ii+1]);
        um = (solver->sqrtRho[ii]*u0 + solver->sqrtRho[ii+1]*u1)/rd;
        //pm = (solver->sqrtRho[ii]*p0 + solver->sqrtRho[ii+1]*p1)/rd
        hm = (solver->sqrtRho[ii]*h0 + solver->sqrtRho[ii+1]*h1)/rd;

        cm = sqrt((hm - um*um/2)*(solver->g-1));

        alp0 = ((p1-p0) - cm*rhom*(u1-u0))/(2*(cm*cm));
        alp1 = (solver->U[ii+1][0] - solver->U[ii][0]) - (p1 - p0)/(cm*cm);
        alp2 = ((p1-p0) + cm*rhom*(u1-u0))/(2*(cm*cm));

        vec0[0] = 1;
        vec0[1] = um - cm;
        vec0[2] = hm - um*cm;

        vec1[0] = 1;
        vec1[1] = um;
        vec1[2] = um*um/2;

        vec2[0] = 1;
        vec2[1] = um + cm;
        vec2[2] = hm + um*cm;

        // run
        e0 = abs(um - cm);
        e0 = solverEntropyFix(solver, e0);

        e1 = abs(um);

        e2 = abs(um + cm);
        e2 = solverEntropyFix(solver, e2);

        solverUpdateCFL(solver, e0, e2);

        for(jj=0;jj<3;jj++){
            solver->flux[ii][jj] = (solver->F[ii][jj] + solver->F[ii+1][jj])/2;
            solver->flux[ii][jj] -= (vec0[jj]*e0*alp0 + vec1[jj]*e1*alp1 + vec2[jj]*e2*alp2)/2;
        };

        /*
        // check
        e0 = (um - cm);
        e1 = (um);
        e2 = (um + cm);

        for(jj=0;jj<3;jj++){
            solver->flux[ii][jj] = (solver->F[ii+1][jj] - solver->F[ii][jj])/2;
            solver->flux[ii][jj] -= (vec0[jj]*e0*alp0 + vec1[jj]*e1*alp1 + vec2[jj]*e2*alp2)/2;
        };

        printf("%f, %f, %f,\n", solver->flux[ii][0], solver->flux[ii][1], solver->flux[ii][2]);
        */
    };

}

void solverPropagate(struct solverStruct* solver){

    int ii, jj;
    FTYPE** aux;
    FTYPE u;

    solverCalcF(solver);

    solverCalcFlux(solver);

    for(jj=0; jj<3; jj++){
        solver->newU[0][jj] = solver->U[0][jj] - solver->cc*(solver->flux[0][jj]);
    };
    // Left boundary condition
    // Interpolate massa flux
    solver->newU[0][1] = solver->newU[1][1]/3;

    for(ii=1; ii<(solver->N-1); ii++){
        for(jj=0; jj<3; jj++){
            solver->newU[ii][jj] = solver->U[ii][jj] - solver->cc*(solver->flux[ii][jj] - solver->flux[ii-1][jj] );

        };

    };

    for(jj=0; jj<3; jj++){
        solver->newU[solver->N-1][jj] = solver->U[solver->N-1][jj] + solver->cc*(solver->flux[solver->N-2][jj]);

    };
    // Right boundary condition
    // Interpolate mass flux
    solver->newU[solver->N-1][1] = solver->newU[solver->N-2][1]/3;

    aux = solver->U;
    solver->U = solver->newU;
    solver->newU = aux;

}

void solverSimulate(struct solverStruct* solver){

    int ii;
    float save = solver->saveStep;

    for(ii=0; ii<solver->Ns; ii++){
        solver->CFL = 0.0;
        solverPropagate(solver);
        printf("\nCFL in the iteration %i: %f", ii, solver->CFL);

        if(ii > save){
            printf("\nCalculating solution %i.", ii);
            solverPrint(solver, ii);
            save += solver->saveStep;
        };

    };
    printf("\nCalculating solution %i.", ii);
    solverPrint(solver, ii);

}

FTYPE solverEntropyFix(struct solverStruct* solver, FTYPE e){

    if(e < solver->eLim){

        e = (solver->eLim + e*e/solver->eLim);

    };

    return e;

}

void solverUpdateCFL(struct solverStruct* solver, FTYPE e0, FTYPE e2){

    FTYPE newCFL;

    newCFL = e0*solver->cc;
    if(newCFL > solver->CFL){
        solver->CFL = newCFL;
    };

    newCFL = e2*solver->cc;
    if(newCFL > solver->CFL){
        solver->CFL = newCFL;
    };

}
