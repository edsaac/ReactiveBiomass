#include "fvCFD.H"

class MyClass
{
    public:
        dimensionedScalar* m_x = nullptr;
        dimensionedScalar* m_y = nullptr;
  
    MyClass(dimensionedScalar* x, dimensionedScalar* y)
        : m_x(x),m_y(y)
    {
        this->print();
    }

    void print()
    {
        Foam::Info<< *m_x << endl;
        Foam::Info<< *m_y << endl;
    }
};

int main(int argc, char *argv[])
{
    // vector zeros   (0,0,0);
    // vector numbers (-10,-1,10);
    // Foam::Info << Foam::max(zeros,numbers) << endl;

    // dimensionedVector dimZeros   ("dimZeros",dimensionSet(0,1,-1,0,0),zeros);
    // dimensionedVector dimNumbers ("dimNumbers",dimensionSet(0,1,-1,0,0),numbers);
    // Foam::Info << Foam::max(dimZeros,dimNumbers) << endl;

    dimensionedScalar xOut ("xOut",dimensionSet(0,1,-1,0,0),15.0);
    dimensionedScalar yOut ("yOut",dimensionSet(0,1,-1,0,0),-10.0);
        
    Info<< "---------------------------" <<endl ;
    MyClass myinstance2(&xOut,&yOut);

    Info<< "---------------------------" <<endl ;
    yOut *= 123;
    myinstance2.print();

    return 0;
}
