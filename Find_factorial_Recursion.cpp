/*
 * C++ Program to Find Factorial of a Number using Recursion
 */
#include <cstring>
#include <iostream>
#include <cstdlib>
#define ll long long
using namespace std;
/*
 * Find Factorial of a Number using Recursion
 */
ll fact_recur(int n)
{
    if (n == 0 || n == 1)
        return 1;
    else
        return n * fact_recur(n - 1);
}
/*
 * Main
 */
int main()
{
    int n;
    while (1)
    {
        cout<<"Enter interger to compute factorial(0 to exit): ";
        cin>>n;
        if (n == 0)
            break;
        cout<<fact_recur(n)<<endl;
    }
    return 0;
}
