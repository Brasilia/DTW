#include <vector>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>

using namespace std;

float Min(float a, float b){
    if (a < b)
        return a;
    return b;
}

float Max(float a, float b){
    if (a > b)
        return a;
    return b;
}

struct VectorMult {
    int dim;
    int seriesClass;
    vector<float> **axes;
};

bool ReadFile(const char *filename, vector<VectorMult> *inVectors, int dimensions = 1){
    //Input file
    FILE *fin = NULL;
    fin = fopen(filename, "r");
    if (fin == NULL){
        fprintf(stderr, "File could not be opened. Check file path.");
        return false;
    }

    //Buffers
    const int LINE_MAX_LENGTH = 100000;
    char line[LINE_MAX_LENGTH];
    float f;
    int nlines = 0;
    int seriesClass;

    while (!feof(fin)){
        if(fgets(line, LINE_MAX_LENGTH, fin) == NULL) //não há mais nada para ler?
            break;
        //Variables to track str position
        int nread;
        int stroffset;
        //Reads from string (line) to integer seriesClass
        sscanf(line, "%d%n", &seriesClass, &stroffset);
        //Allocates memory for new entry
        (*inVectors).resize(nlines+1);
        (*inVectors)[nlines].seriesClass = seriesClass;
        (*inVectors)[nlines].axes = (vector<float>**)malloc(sizeof(vector<float>*)*dimensions);
        for (int i = 0; i < dimensions; i++){
            (*inVectors)[nlines].axes[i] = new vector<float>;
        }
        //Reads line until eof, storing float values
        int scan = 1;
        while (scan != EOF){
            //Cicles through all dimensions
            for (int i = 0; i < dimensions; i++){
                scan = sscanf(&line[stroffset], "%f%n", &f, &nread);
                stroffset += nread;
                (*inVectors)[nlines].axes[i]->push_back(f);
            }
        }
        nlines++; //next line
    }
    fclose(fin);
    return true;
}

void PrintDataset(vector<VectorMult> inVectors, int dimensions){
    for (unsigned int i = 0; i < inVectors.size(); i++){
        for (int j = 0; j < dimensions; j++){
            cout << "j " << j << endl;
            for (unsigned int k = 0; k < inVectors[i].axes[j]->size(); k++){
                cout << inVectors[i].axes[j]->at(k) << " ";
            }
        }
        cout << endl;
    }
}

int main (int argc, char** argv){
    int dimensions = 1;
    vector<VectorMult> inVectors;
    if(!ReadFile("treino.txt", &inVectors, dimensions)){
        return -1;
    };

    cout << "vsize = " << inVectors.size() << endl;

    getchar();

    PrintDataset(inVectors, dimensions);


    return 0;
}
