#include "constants.H"
#include "messageStream.H"

#define KBOLTZ Foam::constant::physicoChemical::k
#define PI Foam::constant::mathematical::pi
#define GRAVITY Foam::constant::universal::G

int main(int argc, char *argv[])
{
    Foam::Info<< "Boltzmann constant" <<Foam::endl ;
    Foam::Info<< KBOLTZ <<Foam::endl;

    Foam::Info<< "Pi" <<Foam::endl ;
    Foam::Info<< PI <<Foam::endl;

    Foam::Info<< "Gravity" <<Foam::endl ;
    Foam::Info<< GRAVITY <<Foam::endl;

    return 0;
}
