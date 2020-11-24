
// Restraint-handling code for Privateer
// (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL (https://www.gnu.org/licenses/lgpl.html)
//
// 2013-2019 Jon Agirre
// York Structural Biology Laboratory
// The University of York
// mailto: jon.agirre@york.ac.uk
//

#include <functional>
#include <thread>
#include <atomic>
#include <vector>
#include <memory>
#include <future>
#include <mutex>
#include <queue>

#ifndef PRIVATEER_PARALLELISM_H_INCLUDED
#define PRIVATEER_PARALLELISM_H_INCLUDED


namespace privateer {

	namespace control
	{
        template <typename T>
        class Queue 
        {
            public:
                bool push(T const & value);
                // deletes the retrieved element, do not use for non integral types
                bool pop(T & v); 
                bool empty(); 

            private:
                std::queue<T> m_queue;
                std::mutex    m_mutex;
        };
	}

	class thread_pool 
    {

    public:

        thread_pool() { this->init(); }
        thread_pool(size_t nThreads) { this->init(); this->resize(nThreads); }

        // the destructor waits for all the functions in the queue to be finished
        ~thread_pool() { this->interrupt(false); }

        // get the number of running threads in the pool
        inline size_t size() const { return m_threads.size(); }

        // number of idle threads
        inline size_t n_idle() const { return ma_n_idle; }

        // get a specific thread
        inline std::thread & get_thread( size_t i ) { return *m_threads.at(i); }


        // restart the pool
        void restart()
        {
            this->interrupt(false); // finish all existing tasks but prevent new ones
            this->init(); // reset atomic flags
        }

        // change the number of threads in the pool
        // should be called from one thread, otherwise be careful to not interleave, also with this->interrupt()
        void resize(size_t nThreads);

        // wait for all computing threads to finish and stop all threads
        // may be called asynchronously to not pause the calling thread while waiting
        // if kill == true, all the functions in the queue are run, otherwise the queue is cleared without running the functions
        void interrupt( bool kill = false );
        
        template<typename F, typename... Rest>
        auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))>;


        // run the user's function that accepts argument size_t - id of the running thread. returned value is templatized
        // operator returns std::future, where the user can get the result and rethrow the catched exceptions
        template<typename F>
        auto push(F && f) ->std::future<decltype(f(0))>;


        void thread_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn, bool par = true);

        // void async_par_for(unsigned start, unsigned end, std::function<void(unsigned i)> fn, bool par = true);



    private:

        // deleted
        thread_pool(const thread_pool &);// = delete;
        thread_pool(thread_pool &&);// = delete;
        thread_pool & operator=(const thread_pool &);// = delete;
        thread_pool & operator=(thread_pool &&);// = delete;

        // clear all tasks
        void clear_queue()
        {
            std::function<void(size_t id)> * _f = nullptr;
            while (m_queue.pop(_f))
                delete _f; // empty the queue
        }

        // reset all flags
        void init() { ma_n_idle = 0; ma_kill = false; ma_interrupt = false; }

        // each thread pops jobs from the queue until:
        //  - the queue is empty, then it waits (idle)
        //  - its abort flag is set (terminate without emptying the queue)
        //  - a global interrupt is set, then only idle threads terminate
        void setup_thread( size_t i );
        

        // ----------  =====  ----------

        std::vector<std::unique_ptr<std::thread>>        m_threads;
        std::vector<std::shared_ptr<std::atomic<bool>>>  m_abort;
        control::Queue<std::function<void(size_t id)> *>  m_queue;

        std::atomic<bool>    ma_interrupt, ma_kill;
        std::atomic<size_t>  ma_n_idle;

        std::mutex               m_mutex;
        std::condition_variable  m_cond;
    };
}

#endif
