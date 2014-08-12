#ifndef SPIKE_QUEUE_H
#define SPIKE_QUEUE_H

#include <queue>

#include "nest.h"
#include "scheduler.h"

namespace nest
{
  /**
   * Priority queue for all spikes, first spike first.
   */
  class SpikeQueue {
  public:
    
    SpikeQueue();
    
    /**
     * Add spike to queue.
     * @param  stamp      Delivery time 
     * @param  ps_offset  Precise timing offset of spike time
     * @param  weight     Weight of spike.
     */
    void add_spike(const long_t stamp, const double_t ps_offset, 
		   const double_t weight);
    
    /**
     * Return next spike.
     * @param req_stamp  Request spike with this stamp.  Queue
     *                   should never contain spikes with smaller
     *                   stamps.  Spikes with larger stamps are
     *                   left in queue.
     * @param ps_offset  PS-sense offset of spike time
     * @param weight     Spike weight
     * @returns          true if spike available, false otherwise
     */
    bool get_next_spike(const long_t req_stamp, 
			double_t& ps_offset, double_t& weight);

    /**
     * Clear buffer
     */
    void clear();

  private:        

    /**
     * Information about spike.
     */
    struct SpikeInfo {
      SpikeInfo(long_t stamp, double_t ps_offset, double_t weight);
      
      bool operator<(const SpikeInfo& b) const;
      
      // data elements must not be const, since heap implementation
      // in DEC STL uses operator=().
      long_t   stamp_;     //<! spike's time stamp
      double_t ps_offset_; //<! spike offset is PS sense
      double_t weight_;    //<! spike weight
    };

    std::priority_queue<SpikeInfo> queue_;  //!< queue for spikes

  };
  

  inline
  void SpikeQueue::add_spike(const long_t stamp, const double_t ps_offset, 
			     const double_t weight)
  {
    queue_.push(SpikeInfo(stamp, ps_offset, weight));
  }

  inline
  bool SpikeQueue::get_next_spike(const long_t req_stamp, 
				  double_t& ps_offset, double_t& weight)
  {
    if ( queue_.empty() )
      return false;

    // NOTE: popping the queue will invalidate this reference!
    const SpikeInfo& next = queue_.top();

    if ( next.stamp_ > req_stamp )
      return false;  

    // no old spikes must be left in queue
    assert(next.stamp_ == req_stamp);

    ps_offset = next.ps_offset_;
    weight    = next.weight_;

    queue_.pop();

    return true;
  }
      
  inline
  bool SpikeQueue::SpikeInfo::operator<(const SpikeInfo& b) const
  {
    return stamp_ == b.stamp_ ? ps_offset_ < b.ps_offset_ : stamp_ > b.stamp_;
  }
  
}
#endif
