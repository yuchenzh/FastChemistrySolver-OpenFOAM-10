#include "fvCFD.H"

#include "RREntry.H"
namespace Foam
{


    Ostream& operator << (Ostream& os, const RREntry& re)
    {
        re.write(os);
        return os;
    }

    Istream& operator >> (Istream& is, RREntry& re)
    {
        re.read(is);
        return is;
    }

    bool operator == (const RREntry& re1, const RREntry& re2)
    {
        return (&re1 == &re2)? true: false;
    }

    bool operator != (const RREntry& re1, const RREntry& re2)
    {
        return !(re1 == re2);
    }
}