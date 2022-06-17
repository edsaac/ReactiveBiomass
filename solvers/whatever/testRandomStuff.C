#include "fvCFD.H"

int main(int argc, char *argv[])
{
    vector zeros   (0,0,0);
    vector numbers (-10,-1,10);
    Foam::Info << Foam::max(zeros,numbers) << endl;

    dimensionedVector dimZeros   ("dimZeros",dimensionSet(0,1,-1,0,0),zeros);
    dimensionedVector dimNumbers ("dimNumbers",dimensionSet(0,1,-1,0,0),numbers);
    Foam::Info << Foam::max(dimZeros,dimNumbers) << endl;
}
