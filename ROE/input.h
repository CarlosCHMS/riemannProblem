#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

    typedef double FTYPE;

    struct inputStruct{

        //FTYPE c;
        //FTYPE k;
        int Ns;
        int Np;

        FTYPE n;
        FTYPE saveStep;
        FTYPE p0;
        FTYPE T0;
        FTYPE p1;
        FTYPE T1;
        FTYPE cc;
        FTYPE R;
        FTYPE gamma;

    };

    void inputInit(struct inputStruct* input);

    void inputPrintParameters(struct inputStruct* input);

#endif // INPUT_H_INCLUDED
