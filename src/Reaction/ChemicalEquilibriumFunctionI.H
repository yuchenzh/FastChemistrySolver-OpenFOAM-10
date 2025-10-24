inline void OptReaction::update_Pow_pByRT_SumVki
(
    const double T_
) const
{
    for(auto it = Pow_pByRT_SumVki_I.begin(); it!=Pow_pByRT_SumVki_I.end();it++)
    {
        const int SumVki = it->first;
        const double val = std::pow(this->Pstd/(this->Ru*T_),SumVki);
        it->second = val;
    }
}

inline void OptReaction::update_Pow_pByRT_SumVki2
(
    const double T_
) const
{
    const double PByRT = this->Pstd/(this->Ru*T_);
    Pow_pByRT_SumVki[0] = 1/(PByRT*PByRT);
    Pow_pByRT_SumVki[1] = 1/(PByRT);      
    Pow_pByRT_SumVki[2] = 1;   
    Pow_pByRT_SumVki[3] = PByRT;    
    Pow_pByRT_SumVki[4] = PByRT*PByRT;   
}
