#include "mainwindow.h"
#include <QApplication>
#include <iostream>
#include <stdio.h>


using namespace std;
//test 1: write a function that reverse a string in place
void stdreverse(char* S)
{
    char temp;
    int l = strlen(S);
    char *end_str=S+l-1;
    char *str=S;

    while(end_str>str)
    {
        temp = *end_str;
        *end_str = *str;
        *str = temp;
        end_str--;
        str++;
    }

    std::cout<<"check the string after reverse:"<<std::endl;
    for(size_t i=0;i<l;++i)
    {
        std::cout<<"print char: "<<S[i]<<std::endl;
    }

}


//implement a binary search method on an infinite array
int binarySearch(int sortedArray[],int first, int last, int key)
{
    while(first<=last)
    {
        int mid = (first+last)/2;
        if (key > sortedArray[mid])
                   first = mid + 1;  // repeat search in top half.
               else if (key < sortedArray[mid])
                   last = mid - 1; // repeat search in bottom half.
               else
                   return mid;     // found it. return position /////
     }
     return -(first + 1);    // failed to find key
}

void stairCase(int n)
{
    for(int i=0;i<n;++i)
    {
        int sp =i+1;
        for(int k=1;k<sp;++k)
            cout<<"#";
        cout<<endl;//"\n");
    }
}



int binarySearchFirstOccur(int sortedArray[], int first, int last, int key)
{
    int retVal=-1;
    while(first<=last)
    {
        std::cout<<"does it even come here"<<std::endl;
        int mid = (first+last)/2;
        if(key==sortedArray[mid])
        {
            retVal = mid;
            last=mid-1;
        }
        else if(key<sortedArray[mid])
        {
            last=mid-1;
        }
        else
        {
            first=mid+1;
        }
    }
    return retVal;
}


class Someclass {
public:
    int x;

public:
    Someclass(int xx) : x(xx){}
    Someclass(const Someclass& a) {x=a.x;x++;}
    void operator =(const Someclass& a1) {x = a1.x ; x--;}

};


template <typename T>
class Foo{
    T tVar;
public:
    Foo(T t): tVar(t){ }
    Foo()
    {

    }
};

class FooDerived:public Foo<std::string>{};


void f() throw (float){throw 10.0f;
                      }


class A {
public:
    int data;
    virtual void fun()
    {
    }
    virtual void fun1(){
    }
};

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    A aaa;
    cout<<"size of aaa"<<sizeof(aaa)<<std::endl;

    double foo1=4.0;
    const float& bar = foo1;
    foo1 +=2.0;
    cout<<"check bar:"<<bar<<std::endl;



    int xxx=4,*y;
    y=&xxx;
    (*y)++;
    cout<<"check *y"<<*y<<std::endl;

    char StrA[]="hello";
    char *Str=&StrA[0];
    stdreverse(Str);

    //test the binary search
    int arr[]={1,2,2,3,3,3,4,5,6};
    int x=3;
    int n = sizeof(arr)/sizeof(arr[0]);
    int c = binarySearchFirstOccur(arr,0,n-1,x);
    try
    {
        printf("%d occurs %d times ",x,c);
        f();
        std::cout<<"check c "<<c<<std::endl;

    }catch(...)
    {
        cout<<"check ddd"<<std::endl;
    }

    stairCase(6);

    Someclass aa(4);
    Someclass b =aa;
    cout<<"check b: "<<b.x<<endl;


        FooDerived d;


    //getchar();


    return a.exec();

}


