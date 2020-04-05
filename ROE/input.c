#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "solver.h"
#include "input.h"

void inputInit(struct inputStruct* input){

    /*
    Initialize the inputStruct reading the input file
    */

    FILE* f1;
    int ii, jj, kk;
    char c;
    char s[100];

    f1 = fopen("aux.csv", "r");

    ii = 0;
    jj = 0;
    kk = 0;

    //Read input file
    while(c!=EOF){
        c = getc(f1);

        if(c == ','){
            //When the string gets a coma:
            s[jj] = '\0';
            jj = 0;
            //printf("%s, ", s);
            //printf("\ni: %i, k: %i, s: %s", ii, kk, s);
            if(ii==0){

                if(kk==1){
                    input->Ns = strtod(s, NULL);
                };

            }else if(ii==1){

                if(kk==1){
                    input->Np = strtod(s, NULL);
                };

            }else if(ii==2){

                if(kk==1){
                    input->n = strtod(s, NULL);
                };

            }else if(ii==3){

                if(kk==1){
                    input->saveStep = strtod(s, NULL);
                };

            }else if(ii==4){

                if(kk==1){
                    input->p0 = strtod(s, NULL);
                };

            }else if(ii==5){

                if(kk==1){
                    input->T0 = strtod(s, NULL);
                };

            }else if(ii==6){

                if(kk==1){
                    input->p1 = strtod(s, NULL);
                };

            }else if(ii==7){

                if(kk==1){
                    input->T1 = strtod(s, NULL);
                };

            }else if(ii==8){

                if(kk==1){
                    input->cc = strtod(s, NULL);
                };

            }else if(ii==9){

                if(kk==1){
                    input->R = strtod(s, NULL);
                };

            }else if(ii==10){

                if(kk==1){
                    input->gamma = strtod(s, NULL);
                };

            }else if(ii==11){

                if(kk==1){
                    input->fluxType = strtod(s, NULL);
                };

            };

            kk ++;
        }else if(c == '\n'){
            ii++;
            kk = 0;
        }else{
            if(c != ' '){
                s[jj] = c;
                jj++;
            };
        };
        //printf("\n%i", ii);
    };
    //printf("\noiii");

    fclose(f1);

}

void inputPrintParameters(struct inputStruct* input){

    /*
    Print the inputs to the command shell
    */

    printf("\nParameters:");
    printf("\nNs: %i", input->Ns);
    printf("\nNp: %i", input->Np);
    printf("\nn: %f", input->n);
    printf("\nsaveStep: %f", input->saveStep);
    printf("\np0: %f", input->p0);
    printf("\nT0: %f", input->T0);
    printf("\np1: %f", input->p1);
    printf("\nT1: %f", input->T1);
    printf("\ncc: %f", input->cc);
    printf("\nR: %f", input->R);
    printf("\ngamma: %f", input->gamma);
    printf("\nfluxType: %i", input->fluxType);

};
