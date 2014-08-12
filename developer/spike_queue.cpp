#include "spike_queue.h"
#include "config.h"
                
nest::SpikeQueue::SpikeQueue()
{
}

nest::SpikeQueue::SpikeInfo::SpikeInfo(long_t stamp, double_t ps_offset, double_t weight) :
  stamp_(stamp), 
  ps_offset_(ps_offset), 
  weight_(weight)
{}
