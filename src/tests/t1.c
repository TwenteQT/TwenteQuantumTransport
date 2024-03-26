#include <stdio.h>
int main(){
    double a1=1.3,a2=1.9;
    int i1,i2,i3,i4;
    i1=a1;
    i2=a2;
    i3=ceil(a1);
    i4=ceil(a2);
    printf("%d, %d, %d, %d\n",i1,i2,i3,i4);
    printf("%f, %f, %f, %f\n",ceil(a1),ceil(a2),floor(a1),floor(a2));
    
}