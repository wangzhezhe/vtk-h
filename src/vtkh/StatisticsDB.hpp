#ifndef __STATISTICS_DB_H
#define __STATISTICS_DB_H

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <unistd.h>
#define TIMEINFO timeval

#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <ratio>

#include <vtkh/vtkh_exports.h>
#include <vtkh/utils/StreamUtil.hpp>

namespace vtkh
{

class VTKH_API StopWatch
{
public:
  enum StopWatchType
  {
    REGULAR,
    ALL_HISTORY,
    EPOCH_HISTORY
  };

  StopWatch(StopWatchType swType= StopWatch::REGULAR)
  : Time(0), Running(false), Mode(swType) {}
  ~StopWatch() {}

  void Start()
  {
    this->StartTime = std::chrono::high_resolution_clock::now();
    this->Running = true;
  }

  void EpochStart()
  {
    if (this->Mode != EPOCH_HISTORY)
      throw "Time is not an epoch time";

    this->EpochTime = 0;
  }

  double EpochStop()
  {
    if (this->Mode != EPOCH_HISTORY)
      throw "Time is not an epoch time";

    this->History.push_back(this->EpochTime);
    double tmp = this->EpochTime;
    this->EpochTime = 0;
    return tmp;
  }

  double Stop()
  {
    this->EndTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeDiff = this->EndTime-this->StartTime;
    double dt = timeDiff.count();
    this->Time += dt;

    if (this->Mode == EPOCH_HISTORY)
      this->EpochTime += dt;

    this->Running = false;
    return this->Time;
  }

  void Reset()
  {
    this->Time=0.0;
    this->EpochTime = 0.0;
    this->History.resize(0);
    this->Running = false;
  }

  bool IsRunning() const {return this->Running;}
  double GetTime() const {return this->Time;}
  StopWatchType GetMode() const {return this->Mode;}
  const std::vector<double>& GetHistory() const { return this->History;}

private:
  StopWatchType Mode;
  double EpochTime, Time;
  bool Running, KeepHistory;
  std::chrono::high_resolution_clock::time_point StartTime, EndTime;
  std::vector<double> History;
};


class VTKH_API EventHistory
{
public:
    EventHistory() : t0(-1000) {}
    void SetT0(double t) {t0=t;}

    double getTime()
    {
#ifdef VTKH_PARALLEL
        return MPI_Wtime();
#else
        struct TIMEINFO ti;
        gettimeofday(&ti, 0);
        return (double)(ti.tv_sec + ti.tv_usec/1e6f);
#endif
    }

    void begin() {ts = getTime()-t0;}
    void end() {history.push_back(std::make_pair(ts, getTime()-t0));}

    void Normalize(float tmin, float tmax)
    {
        float dt = tmax-tmin;
        for (long unsigned int i = 0; i < history.size(); i++)
        {
            float t0 = history[i].first;
            float t1 = history[i].second;
            t0 = (t0-tmin)/dt;
            t1 = (t1-tmin)/dt;
            history[i].first = t0;
            history[i].second = t1;
        }
    }

    double ts = 0 , t0 = 0;
    std::vector<std::pair<double,double>> history;
};

class VTKH_API StatisticsDB
{
public:
  static std::map<std::string, StatisticsDB*> Stats;
  static StatisticsDB* GetStats(const std::string& nm);

  template <typename T>
  struct statValue
  {
    statValue() {}
    statValue(std::vector<T> vals)
    {
      if (vals.empty())
      {
        minI = maxI = -1;
        min = max = med = sum = T(-1);
        mean = std_dev = -1.f;
      }
      else
      {
        sum = accumulate(vals.begin(), vals.end(), T(0));

        //Get min/max info.
        auto res = minmax_element(vals.begin(), vals.end());
        minI = res.first - vals.begin();
        maxI = res.second - vals.end();
        min = (*res.first);
        max = (*res.second);

        //compute mean/median
        int n = vals.size();
        sort(vals.begin(), vals.end());
        mean = (float)sum / (float)n;
        med = vals[vals.size()/2];

        //compute standard deviation
        float x = 0;
        for (int i = 0; i < n; i++)
          x += (vals[i]-mean)*(vals[i]-mean);

        std_dev = sqrt(x/(float)n);
        values = vals;
      }
    }

    int minI, maxI;
    T min,max, med, sum;
    float mean, std_dev;
    std::vector<T> values;

    friend std::ostream &
    operator<<(std::ostream &os, const statValue<T> &s)
    {
      return os<<"AVG: "<<s.mean<<" MED: "<<s.med<<" ("<<s.min<<","<<s.max<<":"<<s.std_dev<<")";
    }
  };

  StatisticsDB() : statsComputed(false) {}
  StatisticsDB(const StatisticsDB &s)
  {
    statsComputed = s.statsComputed;
    timers = s.timers;
    counters = s.counters;
    timerStats = s.timerStats;
    counterStats = s.counterStats;
  }

  ~StatisticsDB() {}

  void insert(const std::vector<StatisticsDB> &v)
  {
    statsComputed = false;
    for (long unsigned int i = 0; i < v.size(); i++)
    {
      const StatisticsDB &s = v[i];
      for (auto ti = s.timers.begin(); ti != s.timers.end(); ti++)
      {
        if (timers.find(ti->first) != timers.end())
          throw std::runtime_error("Duplicate timer: "+ti->first);
        timers[ti->first] = ti->second;
      }

      for (auto ci = s.counters.begin(); ci != s.counters.end(); ci++)
      {
        if (counters.find(ci->first) != counters.end())
          throw std::runtime_error("Duplicate counter: "+ci->first);
        counters[ci->first] = ci->second;
      }

      for (auto ei = s.events.begin(); ei != s.events.end(); ei++)
      {
        if (events.find(ei->first) != events.end())
          throw std::runtime_error("Duplicate event: "+ei->first);
        events[ei->first] = ei->second;
      }
    }
  }



  //Timers.
  void AddTimer(const std::string &nm, StopWatch::StopWatchType swType)
  {
    if (timers.find(nm) != timers.end())
      throw nm + " timer already exists!";
    timers[nm] = StopWatch(swType);
    timers[nm].Reset();
  }
  void EpochStart(const std::string &nm) {vt(nm); timers[nm].EpochStart();}
  void EpochStop(const std::string &nm) {vt(nm); timers[nm].EpochStop();}
  void Start(const std::string &nm) {vt(nm); timers[nm].Start();}
  float Stop(const std::string &nm) {vt(nm); return timers[nm].Stop();}
  float Time(const std::string &nm) {vt(nm); return timers[nm].GetTime();}
  void Reset(const std::string &nm) {vt(nm); timers[nm].Reset();}

  //Counters.
  void AddCounter(const std::string &nm)
  {
    if (counters.find(nm) != counters.end())
      throw nm + " counter already exists!";
    counters[nm] = 0;
  }
  void Increment(const std::string &nm) {vc(nm); counters[nm]++;}
  void Increment(const std::string &nm, unsigned long val) {vc(nm); counters[nm]+=val;}
  unsigned long val(const std::string &nm) {vc(nm); return counters[nm];}

  //Events.
  void AddEvent(const std::string &nm)
  {
    if (events.find(nm) != events.end())
      throw nm + " event already exists!";
    events[nm] = EventHistory();
  }
  void Begin(const std::string &nm) {ve(nm); events[nm].begin();}
  void End(const std::string &nm) {ve(nm); events[nm].end();}
void SetEventT0()
  {
    double t0 = 0.0;
#ifdef VTKH_PARALLEL
    t0 = MPI_Wtime();
#endif
    for (auto it = events.begin(); it != events.end(); it++)
      it->second.SetT0(t0);
  }

  //Output to file
  void DumpStats(const std::string &fname, const std::string &preamble="", bool append=false);

  statValue<float> timerStat(const std::string &nm) {cs(); return timerStats[nm];}
  statValue<unsigned long> counterStat(const std::string &nm) {cs(); return counterStats[nm];}
  unsigned long totalVal(const std::string &nm) {cs(); return counterStats[nm].sum;}

  void cs() {calcStats();}
  void calcStats()
  {
    if (this->statsComputed)
      return;

    int rank = 0, nProcs = 1;

#ifdef VTKH_PARALLEL
    MPI_Comm mpiComm;
    mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
    MPI_Comm_size(mpiComm, &nProcs);
    MPI_Comm_rank(mpiComm, &rank);
#endif

    std::vector<float> regularTimers, epochTimers;
    std::vector<unsigned long> counterVals;
    int numRegTimers = 0, numEpochTimers = 0, numCounters = 0;

    if (!this->timers.empty())
    {
      //Get the regular timers first.
      std::vector<float> vals;
      std::vector<std::string> names;
      for (const auto& t : this->timers)
        if (true) //t.second.GetMode() == StopWatch::REGULAR)
        {
          names.push_back(t.first);
          vals.push_back(t.second.GetTime());
        }

      numRegTimers = vals.size();
      if (numRegTimers > 0)
      {
        regularTimers.resize(numRegTimers*nProcs, 0.0f);
        for (int i = 0; i < numRegTimers; i++)
          regularTimers[rank*numRegTimers + i] = vals[i];

#ifdef VTKH_PARALLEL
        //Send all this goulash to rank 0.
        if (nProcs > 1)
        {
          std::vector<float> tmp(regularTimers.size(), 0.0f);
          MPI_Reduce(&(regularTimers[0]), &(tmp[0]), regularTimers.size(), MPI_FLOAT, MPI_SUM, 0, mpiComm);
          regularTimers = tmp;
        }
#endif
        //create the stats.
        if (rank == 0)
        {
          for (int i = 0; i < numRegTimers; i++)
          {
            for (int j = 0; j < nProcs; j++)
              this->TimerRankValues[names[i]].push_back(regularTimers[i*nProcs + j]);

            this->timerStats[names[i]] = statValue<float>(this->TimerRankValues[names[i]]);
          }
        }
      }
    }

    //epoch timers
    std::vector<int> numEpochs;
    std::vector<std::string> epochNames;
    for (const auto& t : this->timers)
    {
      if (t.second.GetMode() == StopWatch::EPOCH_HISTORY)
      {
        const auto& history = t.second.GetHistory();
        numEpochs.push_back(history.size());
        epochNames.push_back(t.first);
        numEpochTimers++;
      }
    }

    if (numEpochTimers > 0)
    {
      std::vector<int> historySz(nProcs, 0);
      for (int i = 0; i < numEpochTimers; i++)
      {
        int maxVal = 0;
#ifdef VTKH_PARALLEL
        MPI_Allreduce(&(numEpochs[i]), &maxVal, 1, MPI_INT, MPI_MAX, mpiComm);
#endif
        std::vector<float> hvals(nProcs*maxVal, 0.0f);
        const auto& history = this->timers[epochNames[i]].GetHistory();
        for (long unsigned int j = 0; j < history.size(); j++)
          hvals[rank*maxVal + j] = history[j];
#ifdef VTKH_PARALLEL
        if (nProcs > 1)
        {
          std::vector<float> tmp(hvals.size(), 0.0f);
          MPI_Reduce(&(hvals[0]), &(tmp[0]), hvals.size(), MPI_FLOAT, MPI_SUM, 0, mpiComm);
          hvals = tmp;
        }
#endif
        if (rank == 0)
        {
          this->EpochRankValues[epochNames[i]].resize(nProcs);
          for (int r = 0; r < nProcs; r++)
          {
            for (int j = 0; j < maxVal; j++)
            {
              float v = hvals[r*maxVal + j];
              if (v > 0.0f)
                this->EpochRankValues[epochNames[i]][r].push_back(v);
            }
          }
        }
      }
    }

    if (!this->counters.empty())
    {
      numCounters = this->counters.size();
      counterVals.resize(numCounters*nProcs, 0);

      auto it = this->counters.begin();
      std::vector<std::string> names;
      for (int i = 0; i < numCounters; i++, it++)
      {
        names.push_back(it->first);
        counterVals[rank*numCounters + i] = it->second;
      }

#ifdef VTKH_PARALLEL
      if (nProcs > 1)
      {
          std::vector<unsigned long> tmp(counterVals.size(), 0);
          MPI_Reduce(&(counterVals[0]), &(tmp[0]), counterVals.size(), MPI_UNSIGNED_LONG, MPI_SUM, 0, mpiComm);
          counterVals = tmp;
      }
#endif

      if (rank == 0)
      {
        for (int i = 0; i < numCounters; i++)
        {
          for (int j = 0; j < nProcs; j++)
            this->CounterRankValues[names[i]].push_back(counterVals[j*numCounters + i]);
          this->counterStats[names[i]] = statValue<unsigned long>(this->CounterRankValues[names[i]]);
        }
      }
    }

    if (!this->events.empty())
    {
      int sz = events.size();

      //Normalize all the values.
      std::vector<float> vals0(nProcs,0.0f), valsN(nProcs, 0.0f);

      //find min/max per rank.
      float myMin = std::numeric_limits<float>::max();
      float myMax = std::numeric_limits<float>::min();
      int myMaxSz = -1;

      auto it = events.begin();
      for (int i = 0; i < sz; i++, it++)
      {
        int n = it->second.history.size();
        if (n == 0)
          continue;
        float v0 = it->second.history[0].first;
        float vn = it->second.history[n-1].second;

        if (v0 < myMin) myMin = v0;
        if (vn > myMax) myMax = vn;
        if (n > myMaxSz) myMaxSz = n;
      }

      float allMin, allMax;
      int allMaxSz;
#ifdef VTKH_PARALLEL
      MPI_Allreduce(&myMin, &allMin, 1, MPI_FLOAT, MPI_MIN, mpiComm);
      MPI_Allreduce(&myMax, &allMax, 1, MPI_FLOAT, MPI_MAX, mpiComm);
      MPI_Allreduce(&myMaxSz, &allMaxSz, 1, MPI_INT, MPI_MAX, mpiComm);

      int buffSz = allMaxSz*2 + 1, tag = 0;
      for (it = events.begin(); it != events.end(); it++)
      {
        //Normalize timings.
        //it->second.Normalize(allMin, allMax);

        //Rank 0 will recv everything.
        if (rank == 0)
        {
          eventStats.resize(nProcs);
          for (int i = 0; i < nProcs; i++)
            eventStats[i][it->first] = EventHistory();

          std::vector<float> buff(buffSz);
          for (int i = 1; i < nProcs; i++)
          {
            MPI_Status status;
            MPI_Recv(&buff[0], buffSz, MPI_FLOAT, i, tag, mpiComm, &status);
            int n = int(buff[0]);
            for (int j = 0; j < n; j+=2)
              eventStats[i][it->first].history.push_back(std::make_pair(buff[1+j], buff[1+j+1]));
          }
          //Put rank 0 data into global stats.
          eventStats[0][it->first] = it->second;
        }
        else
        {
          std::vector<float> evData(buffSz, 0.0f);
          int sz = it->second.history.size();

          evData[0] = sz*2;
          for (int j = 0; j < sz; j++)
          {
            evData[1+j*2+0] = it->second.history[j].first;
            evData[1+j*2+1] = it->second.history[j].second;
          }
          MPI_Send(&evData[0], buffSz, MPI_FLOAT, 0, tag, mpiComm);
        }
      }
#endif
    }

    this->statsComputed = true;

    }

private:
  void vt(const std::string &nm) const
  {
    if (timers.find(nm) == timers.end())
    {
      std::string msg = nm + " timer not found.";
      std::cerr<<"Error: "<<msg<<std::endl;
      throw msg;
    }
  }

  void vc(const std::string &nm) const
  {
    if (counters.find(nm) == counters.end())
    {
      std::string msg = nm + " counter not found.";
      std::cerr<<"Error: "<<msg<<std::endl;
      throw msg;
    }
  }

  void ve(const std::string &nm) const
  {
    if (events.find(nm) == events.end())
    {
      std::string msg = nm + " event not found.";
      std::cerr<<"Error: "<<msg<<std::endl;
      throw msg;
    }
  }

  std::map<std::string, StopWatch> timers;
  std::map<std::string, EventHistory> events;
  std::map<std::string, unsigned long> counters;

  bool statsComputed;
  std::map<std::string, statValue<float> > timerStats, epochTimerStats;
  std::map<std::string, statValue<unsigned long> > counterStats;
  std::vector<std::map<std::string,EventHistory>> eventStats;
  std::map<std::string, std::vector<float>> TimerRankValues;
  std::map<std::string, std::vector<std::vector<float>>> EpochRankValues;
  std::map<std::string, std::vector<unsigned long>> CounterRankValues;

  std::ofstream outputStream;
};

#ifdef VTKH_ENABLE_STATISTICS
#define ADD_COUNTER(nm) StatisticsDB::GetStats("stats")->AddCounter(nm)
#define COUNTER_INC(nm, val) StatisticsDB::GetStats("stats")->Increment(nm, val)

#define ADD_TIMER(nm) StatisticsDB::GetStats("stats")->AddTimer(nm, vtkh::StopWatch::REGULAR)
#define ADD_EPOCH_TIMER(nm) StatisticsDB::GetStats("stats")->AddTimer(nm, vtkh::StopWatch::EPOCH_HISTORY)

#define EPOCH_TIMER_START(nm) StatisticsDB::GetStats("stats")->EpochStart(nm)
#define EPOCH_TIMER_STOP(nm) StatisticsDB::GetStats("stats")->EpochStop(nm)
#define TIMER_START(nm) StatisticsDB::GetStats("stats")->Start(nm)
#define TIMER_STOP(nm) StatisticsDB::GetStats("stats")->Stop(nm)
#define DUMP_STATS(fname) StatisticsDB::GetStats("stats")->DumpStats(fname)

#define ADD_EVENT(nm) StatisticsDB::GetStats("stats")->AddEvent(nm)
#define SET_EVENT_T0() StatisticsDB::GetStats("stats")->SetEventT0()
#define EVENT_BEGIN(nm) StatisticsDB::GetStats("stats")->Begin(nm)
#define EVENT_END(nm) StatisticsDB::GetStats("stats")->End(nm)
#else
#define ADD_COUNTER(nm)
#define COUNTER_INC(nm, val)

#define ADD_TIMER(nm)
#define TIMER_START(nm)
#define TIMER_STOP(nm)
#define DUMP_STATS(fname)

#define ADD_EVENT(nm)
#define SET_EVENT_T0()
#define EVENT_BEGIN(nm)
#define EVENT_END(nm)
#endif

}

#endif //__STATISTICS_DB_H
