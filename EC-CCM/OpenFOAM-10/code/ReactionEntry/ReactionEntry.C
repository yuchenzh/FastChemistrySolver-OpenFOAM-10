#include "fvCFD.H"

#include "ReactionEntry.H"
namespace Foam
{
    Foam::word num2Word(const Foam::label& num, const Foam::label& digits)
    {
        word tempWord = std::to_string(num);
        tempWord.insert(0, digits - tempWord.size(), '0');
        return tempWord;
    }


    Ostream& operator << (Ostream& os, const ReactionEntry& re)
    {
        re.write(os);
        return os;
    }

    Istream& operator >> (Istream& is, ReactionEntry& re)
    {
        re.read(is);
        return is;
    }

    bool operator == (const ReactionEntry& re1, const ReactionEntry& re2)
    {
        return (&re1 == &re2)? true: false;
    }

    bool operator != (const ReactionEntry& re1, const ReactionEntry& re2)
    {
        return !(re1 == re2);
    }

    ReactionEntry operator+ (const ReactionEntry& re1, const ReactionEntry& re2)
    {
        scalar totalCount = re1.count + re2.count;
        scalar re1Ratio = scalar(re1.count)/totalCount;
        scalar re2Ratio = scalar(re2.count)/totalCount;

        if (re1.debug && re2.debug)
        {
            scalarField re1Ystd = (re1.count == 1) ? scalarField(re1.Y.size(),0.0) : re1.Ystd;
            scalarField re2Ystd = (re2.count == 1) ? scalarField(re2.Y.size(),0.0) : re2.Ystd;
            scalarField M1 = sqr(re1Ystd) * (re1.count - 1);
            scalarField M2 = sqr(re2Ystd) * (re2.count - 1);
            scalarField avg = re1.Y * re1Ratio + re2.Y * re2Ratio;
            label newCount = re1.count + re2.count;
            M1 = M1 + M2 + sqr(re1.Y - re2.Y)/(1/scalar(re1.count) + 1/scalar(re2.count));
            return ReactionEntry(avg, 
                            re1.T*re1Ratio + re2.T*re2Ratio, 
                            re1.p*re1Ratio + re2.p*re2Ratio, 
                            re1.dtChem*re1Ratio + re2.dtChem*re2Ratio,
                            re1.rho0*re1Ratio + re2.rho0*re2Ratio,
                            re1.rho*re1Ratio + re2.rho*re2Ratio, 
                            newCount,
                            true,
                            sqrt(M1/(newCount - 1))
                            );
        }


        return ReactionEntry(re1.Y*re1Ratio + re2.Y*re2Ratio, 
                            re1.T*re1Ratio + re2.T*re2Ratio, 
                            re1.p*re1Ratio + re2.p*re2Ratio, 
                            re1.dtChem*re1Ratio + re2.dtChem*re2Ratio,
                            re1.rho0*re1Ratio + re2.rho0*re2Ratio,
                            re1.rho*re1Ratio + re2.rho*re2Ratio, 
                            totalCount);
    }

    void mergeReactionEntries(reactionEntries& dest, const reactionEntries& target)
    {
        for (auto it = target.begin(); it != target.end(); it++)
        {
            if (dest.found(it.key()))
            {
                dest.set(it.key(), dest.find(it.key())() + it());
            }
            else
            {
                dest.set(it.key(), it());
            }
        }
    }
}