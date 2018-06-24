#include <vector>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <ctime>

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

//Reads file to data structure
bool ReadFile(const char *filename, vector<VectorMult> *inVectors, int *maxSize, int dimensions = 1){
    //Input file
    FILE *fin = NULL;
    fin = fopen(filename, "r");
    if (fin == NULL){
        fprintf(stderr, "File \"%s\" could not be opened. Check file path.", filename);
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
        *maxSize = max(*maxSize, (int)((*inVectors)[nlines].axes[0]->size()) );
        nlines++; //next line
    }
    fclose(fin);
    return true;
}

void FreeVectors(vector<VectorMult> *v, int dimension = 1){
    for (unsigned int i = 0; i < v->size(); i++){
        for(int j = 0; j < dimension; j++){
            delete( (*v)[i].axes[j] );
        }
        free( (*v)[i].axes );
    }
}

//Prints data from data structure
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

float PointDist(float m, float n){
    return abs(m-n);
    //return pow((m-n),2);
}

float DTWDistance(vector<float> m, vector<float> n, float** warp){
    //mSize+1 rows, nSize+1 cols
    m.insert(m.begin(), 0);
    n.insert(n.begin(), 0);
    int mSize = m.size();
    int nSize = n.size();

    for (int i = 1; i < mSize; i++){
        for (int j = 1; j <nSize; j++){
            float cost = PointDist(m[i], n[j]);
            warp[i][j] = cost + min(warp[i-1][j-1], min(warp[i][j-1], warp[i-1][j]));
        }
    }
    float distance = warp[mSize-1][nSize-1];
    return distance;
}

void AllocWarpM(float*** warpM, int maxSizeH, int maxSizeW){
    *warpM = (float**)malloc(sizeof(float*)*maxSizeH);
    for (int i = 0; i < maxSizeH; i++){
        (*warpM)[i] = (float*)malloc(sizeof(float)*maxSizeW);
    }
    for (int i = 1; i < maxSizeH; i++){
        (*warpM)[i][0] = HUGE_VAL;
    }
    for (int i = 1; i < maxSizeW; i++){
        (*warpM)[0][i] = HUGE_VAL;
    }
    (*warpM)[0][0] = 0;
}

void PrintWarpM(float** warpM, int maxSizeH, int maxSizeW){
    for (int i = 0; i < maxSizeH; i++){
        for (int j = 0; j < maxSizeW; j++){
            printf("%.2f ", warpM[i][j]);
        }
        cout << endl;
    }
}

void FreeWarpM(float** warpM, int maxSizeH){
    for (int i = 0; i < maxSizeH; i++){
        free(warpM[i]);
    }
    free(warpM);
}

int main (int argc, char** argv){

    int dimensions = 1;
    int maxSizeH = 0, maxSizeW = 0;
    vector<VectorMult> trainVectors;
    vector<VectorMult> testVectors;
    if(!ReadFile("treino.txt", &trainVectors, &maxSizeH, dimensions)){
        return -1;
    };
    if(!ReadFile("teste.txt", &testVectors, &maxSizeW, dimensions)){
        return -1;
    };
    maxSizeH ++;
    maxSizeW ++;

    cout << maxSizeH << " maxH/maxW " << maxSizeW << endl;
    float **warpM = NULL;
    AllocWarpM(&warpM, maxSizeH, maxSizeW);

    cout << "vsize = " << trainVectors.size() << endl;
    cout << "vsize = " << testVectors.size() << endl;

    //getchar();

    clock_t time = clock();
    //Compare the test dataset to the train dataset
    int accCounter = 0;
    for (unsigned int i = 0; i < testVectors.size(); i++){ //test cases
        float minDist = HUGE_VAL;
        int nearestNClass = 0;
        int closest = -1;
        for (unsigned int j = 0; j < trainVectors.size(); j++){ //templates
            float dist = DTWDistance(*trainVectors[j].axes[0], *testVectors[i].axes[0], warpM);
            if (dist < minDist){
                minDist = dist;
                nearestNClass = trainVectors[j].seriesClass;
                closest = j;
            }
        }
        //cout << "Test case: " << i << "\nClosest: " << closest << "\nNNClass: " << nearestNClass << "\nTrueClass: " << testVectors[i].seriesClass << endl;
        if (nearestNClass == testVectors[i].seriesClass){
            accCounter++;
        }
        //getchar();
    }
    time = (clock() - time);
    float accuracy = (float)accCounter/testVectors.size();
    cout << "Accuracy: " << accuracy << endl;
    cout << "Total Time: " << time/1000 << "s" << endl;
    cout << "Avg Time: " << ((float)time/testVectors.size()) << "ms" << endl;

    FreeWarpM(warpM, maxSizeH);
    FreeVectors(&testVectors, dimensions);
    FreeVectors(&trainVectors, dimensions);

    getchar();
    return 0;
}
