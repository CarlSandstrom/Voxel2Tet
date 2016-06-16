#include "TimeStamp.h"
namespace voxel2tet
{

TimeStamp::TimeStamp()
{
    inittime = clock();
}

void TimeStamp::Initialize(std::string Message)
{
    Table.clear();
    StartTimer(Message);
}

void TimeStamp::StartTimer(std::string Message)
{
    this->CurrentMessage = Message;
    t0 = clock();
}

void TimeStamp::StopTimer()
{
    tend = clock();
    std::pair<double, std::string> Item;
    Item.first = (double)(tend-t0) / CLOCKS_PER_SEC;
    Item.second = this->CurrentMessage;
    this->Table.push_back(Item);
}

std::vector<std::pair<double, std::string>> TimeStamp::GetTable()
{
    return Table;
}

void TimeStamp::PrintTable()
{
    for (std::pair<double, std::string> Item: this->Table) {
        printf("%s\t%f\n", Item.second.c_str(), Item.first);
    }
}

int TimeStamp::GiveTotalTime()
{
    return clock()-inittime;
}

}
