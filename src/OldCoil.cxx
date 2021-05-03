#include "../include/OldCoil.h"

#include <vector>

std::vector<OldCoil> backup_vec;
std::vector<thread_param> ar;
thread_pool tp;

Type OldCoil::ECalc(Type r, Type Theta){

    return Calc(r, Theta, EAsyncThread);
}

Type OldCoil::BhCalc(Type r, Type Theta){

    return Calc(r, Theta, BhAsyncThread);
}

Type OldCoil::BzCalc(Type r, Type Theta){

    return Calc(r, Theta, BzAsyncThread);
}

Type OldCoil::Calc(Type r, Type Theta, std::function<long unsigned int(int, void*)> func)
{
    int stop = 0, remainder;
    Type ret = 0;
    Type tmp = 0;
    if(ar.size() != tp.size() || backup_vec[0] != *this)
    {
        backup_vec.resize(tp.size());
        for(int it = 0; it < backup_vec.size(); ++it)
        {
            backup_vec[it] = *this;
            for(int it2 = 0; it2 < 4; it2++)
                for(int it3 = 0; it3 < 4; it3++)
                    backup_vec[it].DWeights[it2][it3] = this->DWeights[it2][it3];
            backup_vec[it].FiWeights = this->FiWeights;
            backup_vec[it].EdgeList1 = this->EdgeList1;
            backup_vec[it].EdgeList2 = this->EdgeList2;
        }

        ar.resize(tp.size());
    }
    ZeroMemory(&ar.front(), ar.size() * sizeof(thread_param));
    remainder = int(Type(a) / Type(IncA)) % tp.size();
    for(int it = 0; it < tp.size(); ++it)
    {
        ar[it].r = r;
        ar[it].Theta = Theta;
        ar[it].c = &backup_vec[it];
        ar[it].completed = false;
        if(it == tp.size() - 1)
        {
            ar[it].dR_start = tmp;
            ar[it].dR_stop = a;
            break;
        }
        ar[it].dR_start = tmp;
        if(remainder)
        {
            tmp += (int(Type(a) / Type(IncA) / tp.size()) + 1) * IncA;
            --remainder;
        }
        else
            tmp += int(Type(a) / Type(IncA) / tp.size()) * IncA;
        ar[it].dR_stop = tmp;
        ar[it].dR_stop -= IncA / 4;
    }


    for(int it = 0; it < tp.size(); ++it)
        tp.push(func, &(ar[it]));

    while(stop < tp.size())
    {
        stop = 0;
        for(auto it : ar)
            if(it.completed)
                ++stop;
    }

    for(auto it : ar)
    {
        ret += it.res;
    }

    return ret;
}