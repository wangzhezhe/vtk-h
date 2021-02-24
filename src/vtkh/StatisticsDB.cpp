#include <vtkh/vtkh.hpp>
#include <vtkh/StatisticsDB.hpp>

namespace vtkh
{
std::map<std::string, StatisticsDB*> StatisticsDB::Stats;

StatisticsDB*
StatisticsDB::GetStats(const std::string& nm)
{
  if (Stats.find(nm) == Stats.end())
  {
    StatisticsDB* stats = new StatisticsDB;
    //stats->SetCommunicator(MPI_Comm_f2c(vtkh::GetMPICommHandle()));
    Stats[nm] = stats;
  }

  return Stats[nm];
}

void
StatisticsDB::DumpStats(const std::string &fname, const std::string &preamble, bool append)
{
  this->calcStats();

    int numRanks = 1;
#ifdef VTKH_PARALLEL
    int rank = vtkh::GetMPIRank();
    numRanks = vtkh::GetMPISize();
    if (rank != 0)
        return;
#endif

    if (!append || !outputStream.is_open())
        outputStream.open(fname, std::ofstream::out);

    if (!preamble.empty())
        outputStream<<preamble;

    if (!this->timers.empty())
    {
        outputStream<<"TIMERS:"<<std::endl;
        for (auto &ti : this->timers)
            outputStream<<ti.first<<": "<<ti.second.GetTime()<<std::endl;
        outputStream<<std::endl;
        outputStream<<"TIMER_STATS"<<std::endl;
        for (auto &ti : this->timers)
            outputStream<<ti.first<<" "<<this->timerStat(ti.first)<<std::endl;

        outputStream<<"TIMER_RANK"<<std::endl;
        for (const auto& it : this->TimerRankValues)
        {
          for (int i = 0; i < numRanks; i++)
            outputStream<<it.first<<"_"<<i<<" "<<it.second[i]<<std::endl;
        }

        outputStream<<"EPOCH_TIMER_RANK"<<std::endl;
        for (const auto& it : this->EpochRankValues)
        {
          for (int i = 0; i < numRanks; i++)
            outputStream<<it.first<<"_"<<i<<" "<<it.second[i]<<std::endl;
        }
    }
    if (!this->counters.empty())
    {
        outputStream<<std::endl;
        outputStream<<"COUNTERS:"<<std::endl;
        for (auto &ci : this->counters)
            outputStream<<ci.first<<" "<<this->totalVal(ci.first)<<std::endl;
        outputStream<<std::endl;
        outputStream<<"COUNTER_STATS"<<std::endl;
        for (auto &ci : this->counters)
            outputStream<<ci.first<<" "<<this->counterStat(ci.first)<<std::endl;

        outputStream<<"COUNTER_RANK"<<std::endl;
        for (const auto& it : this->CounterRankValues)
        {
          for (int i = 0; i < numRanks; i++)
            outputStream<<it.first<<"_"<<i<<" "<<it.second[i]<<std::endl;
        }
    }
    if (!this->eventStats.empty())
    {
        std::string eventFname = fname + ".events";
        std::ofstream eventStream;
        eventStream.open(eventFname, std::ofstream::out);

        for (int i = 0; i < numRanks; i++)
        {
            auto events = eventStats[i];
            for (auto ei : events)
            {
                auto eventNm = ei.first;
                auto eventHistory = ei.second;
                eventStream<<eventNm<<"_"<<i<<",";
                for (auto &hi : eventHistory.history)
                    eventStream<<hi.first<<","<<hi.second<<",";
                eventStream<<eventNm<<std::endl;
            }
        }
        eventStream.close();
    }
}

};
