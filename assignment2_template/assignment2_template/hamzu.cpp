#include<iostream>
using namespace std;

int main() {
    int size=3;
    double mask[]={0.25, 0.5, 0.25};

    double Ixx[][4]={{1,2,1, 3},
    {4,1,5,1},
    {1,3,1,6},
    {7,1,2,1} };
    int width=4;
    int height=4;

    double smix[4][4];
    double sum;
    double productsum;
int n;
    for (int j=0; j<height; j++) {
        for (int i=0; i<width; i++) {
            if (i<(size-1) /2) {
                sum=0.0;
                productsum=0.0;
             n=0;
            for (int l=((size-1)/2)-i; l<size; l++) {
                sum += mask[l];
                productsum = productsum + (mask[l]*Ixx[j][n]);
                n++;
            }
            smix[j][i]=productsum/sum;
        }
        else if (i> width-1-((size-1)/2)) {
            sum=0.0;
            productsum=0.0;
            n=width-1;
            for (int m=(((size-1)/2)+(width-1-i)); m>=0; m--) {
                sum += mask[m];
                productsum += mask[m]*Ixx[j][n];
                n--;
            }
            smix[j][i]=productsum/sum;
        }
        else  {
            sum =0.0;
            productsum=0.0;
            n=0;
            for (int k=0; k<size; k++) {
                sum += mask[k];
                productsum = productsum+ (Ixx[j][(i-((size-1)/2))+n]* mask[k]);
                n++;
            }
            smix[j][i]=productsum/sum;

        }

    }

}


for (int j=0;j<height;j++) {
    for (int i=0; i<width; i++) {
        cout<<smix[j][i]<<" ";
    }
    cout<<endl;
}




}