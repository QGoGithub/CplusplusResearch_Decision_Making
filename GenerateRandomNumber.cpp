//pdf(x) = 1 if x>465
//       = 0 if x<0
//       = x/465 otherwise
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

//This is a sample program to generate a random numbers based on probability desity function of spiner
//pdf(x) = 1 if x>460
//       = 0 if x<0
//       = x/460 otherwise
int N = 10;
int main(int argc, char **argv)
{
    int p = 0;
    for (int i = 0; i < N; i++)
    {
        p = rand() % 500;
        if (p > 465)
            cout << 0 << " ";
        else if (p < 0)
            cout << 0 << " ";
        else
            cout << p * 0.1 / 465 << " ";
    }
    cout << "...";
}
