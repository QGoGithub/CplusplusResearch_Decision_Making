//a C++ Program to check if two matrices can be multiplied

#include<conio.h>
#include<iostream>
#include<math.h>

using namespace std;

int main(int argc, char **argv)
{
    cout<<"Enter the dimension of the matrix:\n ";
    int rowA;cin>>rowA;
    int colA;cin>>colA;

    cout<<"Enter the dimension of the other matrix:\n ";
    int rowB;cin>>rowB;
    int colB;cin>>colB;

    if(colA == rowB)
    {
        cout<<"Matrices are multipliable";
    }
    else
    {
        cout<<"Matrices are not multipliable";
    }
}
