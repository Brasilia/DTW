#include <vector>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <ctime>

using namespace std;

int MAXSIZE_H;
int MAXSIZE_W;

struct VectorMult {
    int dim;
    int seriesClass;
    vector<float> **axes;
};

void AllocWarpM(float*** warpM, int maxSizeH, int maxSizeW){
    maxSizeH = max(maxSizeH, maxSizeW) + 1 +1; //formar matriz quadrada, com folga para casos base
    maxSizeW = maxSizeH; //formar matriz quadrada
    *warpM = (float**)malloc(sizeof(float*)*maxSizeH);
    for (int i = 0; i < maxSizeH; i++){
        (*warpM)[i] = (float*)malloc(sizeof(float)*maxSizeW);
    }
    for (int i = 1; i < maxSizeH; i++){
        (*warpM)[i][0] = HUGE_VALF;
    }
    for (int i = 1; i < maxSizeW; i++){
        (*warpM)[0][i] = HUGE_VALF;
    }
    (*warpM)[0][0] = 0;
}

void PrintWarpM(float** warpM, int maxSizeH, int maxSizeW){
    for (int i = 0; i < maxSizeH; i++){
        for (int j = 0; j < maxSizeW; j++){
            int val;
            if(warpM[i][j] == HUGE_VALF){
                val = 1;
            } else {
                val = 0;
            }
            printf("%d ", val);
        }
        cout << endl;
    }
    cout << endl;
}

void FreeWarpM(float** warpM, int maxSizeH){
    for (int i = 0; i < maxSizeH; i++){
        free(warpM[i]);
    }
    free(warpM);
}

//Prepares the matrix for Sakoe-Chiba, quadratic time
void MaxWarpM(float** warpM, int sizeH, int sizeW){
    sizeH = max(sizeH, sizeW);//FIXME tirar +1
    sizeW = sizeH;
    for(int i = 0; i < MAXSIZE_H+1-1; i++){ //+0??
        for(int j = 0; j < MAXSIZE_W+1-1; j++){
            warpM[i][j] = HUGE_VALF;
        }
    }
    warpM[0][0] = 0;
    //PrintWarpM(*warpM, MAXSIZE_H, MAXSIZE_W);
}

//Prepares the matrix for Sakoe-Chiba, linear time
void MaxWarpM(float** warpM, int sizeH, int sizeW, int bandW){
    int maxSize = max(sizeH, sizeW);
    bandW+=1;
    for(int i = bandW; i < maxSize+1; i++){
        warpM[i-bandW][i] = HUGE_VALF;
        warpM[i][i-bandW] = HUGE_VALF;
    }
    warpM[0][0] = 0;
    //PrintWarpM(*warpM, MAXSIZE_H, MAXSIZE_W);
}

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
                if (scan == EOF) //end of file, abort
                    break;
                stroffset += nread;
                (*inVectors)[nlines].axes[i]->push_back(f);
            }
        }
        *maxSize = max(*maxSize, (int)((*inVectors)[nlines].axes[0]->size()) );
        nlines++; //next line
//        for (int i = 0; i < (int)((*inVectors)[nlines-1].axes[0]->size()); i++){
//            cout << (*(*inVectors)[nlines-1].axes[0])[i] << endl;
//        }
//        cout << endl<<endl;
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
            cout << "j " << j << "  Size: " << inVectors[i].axes[j]->size() << endl;
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

//Calculates point distance between multidimensional vectors, either dependent or independent
float PointDist(VectorMult s1, VectorMult s2, int p1, int p2, int dim, bool dependent = false){
    float sum = 0;
    if (dependent){
            for (int i = 0; i < dim; i++){
            sum+= pow( (*s1.axes[i])[p1] - (*s2.axes[i])[p2], 2 );
        }
        return sqrt(sum);
    } else {
        for (int i = 0; i < dim; i++){
            sum+= abs((*s1.axes[i])[p1] - (*s2.axes[i])[p2]);
        }
        return sum;
    }
}

float DTWDistance(vector<float> m, vector<float> n, float** warp){
    //mSize+1 rows, nSize+1 cols
    int mSize = m.size()+1;
    int nSize = n.size()+1;
    //Makes path
    int i, j;
    for (i = 1; i < mSize; i++){
        for (j = 1; j < nSize; j++){
            float cost = PointDist(m[i-1], n[j-1]);
            warp[i][j] = cost + min(warp[i-1][j-1], min(warp[i][j-1], warp[i-1][j]));
        }
    }
    float distance = warp[i-1][j-1];
    return distance;
}

//With Sakoe-Chiba band
float DTWDistance(vector<float> m, vector<float> n, float** warp, int band){
    //mSize+1 rows, nSize+1 cols
    int mSize = m.size()+1;
    int nSize = n.size()+1;
    int bandW = (band * max(mSize, nSize))/100;
    //Crops bigger series so that the final corner can be reached
    int maxSize = mSize + bandW;
    mSize = min (maxSize, mSize);
    nSize = min (maxSize, nSize);
    //Sets out-of-range cells to infinity
    MaxWarpM(warp, mSize, nSize, bandW);
    //Makes path
    int y, x;
    for (int i = 1; i < mSize; i++){
        for (int j = max(1, i-bandW); j <= min(nSize-1, i+bandW); j++){
            float cost = PointDist(m[i-1], n[j-1]);
            warp[i][j] = cost + min(warp[i-1][j-1], min(warp[i][j-1], warp[i-1][j]));
            y = i;
            x = j;
        }
    }
    float distance = warp[y][x];
    return distance/(mSize);
}

//Multidimensional, no band
float DTWDistance(VectorMult s1, VectorMult s2, int dim, float** warp){
    //mSize+1 rows, nSize+1 cols
    int mSize = s1.axes[0]->size()+1;
    int nSize = s2.axes[0]->size()+1;
//    cout << "mSize: " << mSize << endl;
//    cout << "nSize: " << nSize << endl;
//    cout << s1.axes[0]->size() << " " << s1.axes[1]->size() << " " << s1.axes[2]->size() << " " << endl;
//    cout << s2.axes[0]->size() << " " << s2.axes[1]->size() << " " << s2.axes[2]->size() << " " << endl;
    //Makes path
    int i, j;
    for (i = 1; i < mSize; i++){
        for (j = 1; j < nSize; j++){
            float cost = PointDist(s1, s2, i-1, j-1, dim, false);
//            cout << "cost "<< cost<<endl;
//            getchar();
            warp[i][j] = cost + min(warp[i-1][j-1], min(warp[i][j-1], warp[i-1][j]));
        }
    }
    float distance = warp[i-1][j-1];
    //cout << distance <<endl;
    return distance;
}

//Sakoe-Chiba Band, Multidimensional
float DTWDistance(VectorMult s1, VectorMult s2, int dim, float** warp, int band){
    if (band < 0){ //No Band
        return DTWDistance(s1, s2, dim, warp);
    }
    //mSize+1 rows, nSize+1 cols
    int mSize = s1.axes[0]->size()+1;
    int nSize = s2.axes[0]->size()+1;
    int bandW = (band * max(mSize, nSize))/100;
    //Crops bigger series so that the final corner can be reached
    int maxSize = mSize + bandW;
    mSize = min (maxSize, mSize);
    nSize = min (maxSize, nSize);
    //Sets out-of-range cells to infinity
    MaxWarpM(warp, mSize, nSize, bandW);
    //Makes path
    int x, y;
    for (int i = 1; i < mSize; i++){
        for (int j = max(1, i-bandW); j <= min(nSize-1, i+bandW); j++){
            float cost = PointDist(s1, s2, i-1, j-1, dim, false);
            warp[i][j] = cost + min(warp[i-1][j-1], min(warp[i][j-1], warp[i-1][j]));
            y = i;
            x = j;
        }
    }
    float distance = warp[y][x];
    return distance/(mSize);
}

//Runs test dataset against train dataset and returns accuracy
float DTWTest(vector<VectorMult> trainVectors, vector<VectorMult> testVectors, float**warpM, int dimensions){
    int accCounter = 0;
    for (unsigned int i = 0; i < testVectors.size(); i++){ //test cases
        float minDist = HUGE_VALF;
        int nearestNClass = 0;
        //int closest = -1;
        for (unsigned int j = 0; j < trainVectors.size(); j++){ //templates
            //float dist = DTWDistance(*trainVectors[j].axes[0], *testVectors[i].axes[0], warpM);
            //float dist = DTWDistance(*trainVectors[j].axes[0], *testVectors[i].axes[0], warpM, 20);
            float dist = DTWDistance(trainVectors[j], testVectors[i], dimensions, warpM, 20);
            if (dist < minDist){
                minDist = dist;
                nearestNClass = trainVectors[j].seriesClass;
                //closest = j;
            }
        }
        //cout << "Test case: " << i << "\nClosest: " << closest << "\nNNClass: " << nearestNClass << "\nTrueClass: " << testVectors[i].seriesClass << endl;
        if (nearestNClass == testVectors[i].seriesClass){
            accCounter++;
        }
        //getchar();
    }
    //Returns accuracy
    return (float)accCounter/testVectors.size();
}



int main (int argc, char** argv){

    int dimensions = 3;
    int maxSizeH = 0, maxSizeW = 0;
    vector<VectorMult> trainVectors;
    vector<VectorMult> testVectors;
    if(!ReadFile("teste3D.txt", &trainVectors, &maxSizeH, dimensions)){
        return -1;
    };
    if(!ReadFile("teste3D.txt", &testVectors, &maxSizeW, dimensions)){
        return -1;
    };
    //PrintDataset(trainVectors, dimensions);
    //PrintDataset(testVectors, dimensions);

    maxSizeH ++;
    maxSizeW ++;
    MAXSIZE_H = maxSizeH;
    MAXSIZE_W = maxSizeW;

    cout << maxSizeH << " maxH/maxW " << maxSizeW << endl;
    float **warpM = NULL;
    AllocWarpM(&warpM, maxSizeH, maxSizeW);

    cout << "vsize1 = " << trainVectors.size() << endl;
    cout << "vsize2 = " << testVectors.size() << endl;

    //getchar();

    clock_t time = clock();
    //Compare the test dataset to the train dataset
    float accuracy = DTWTest(trainVectors, testVectors, warpM, dimensions);
    time = (clock() - time);

    cout << "Accuracy: " << accuracy << endl;
    cout << "Total Time: " << time/1000 << "s" << endl;
    cout << "Avg Time: " << ((float)time/testVectors.size()) << "ms" << endl;

    FreeWarpM(warpM, maxSizeH);
    FreeVectors(&testVectors, dimensions);
    FreeVectors(&trainVectors, dimensions);

    getchar();
    return 0;
}
